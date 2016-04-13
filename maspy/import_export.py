from __future__ import print_function, division, unicode_literals
from future.utils import viewkeys, viewvalues, viewitems, listvalues, listitems

try: # python 2.7
    from itertools import izip as zip
except ImportError: # python 3.x series
    pass
###############################################################################
import io
import numpy
import os

from lxml import etree as ETREE

import maspy.auxiliary as aux
import maspy.xml
import maspy.core
import maspy.peptidemethods


##########################################################
### mzml import and export methods #######################
##########################################################
def importMzml(filepath, msrunContainer=None, siAttrFromSmi=None, specfilename=None):
    """Performs a complete import of a mzml file into a maspy MsrunContainer.

    siAttrFromSmi: allow here to specify a custom function that extracts params a from spectrumMetadataItem
    """
    siAttrFromSmi = defaultFetchSiAttrFromSmi if siAttrFromSmi is None else siAttrFromSmi
    if msrunContainer is None:
        msrunContainer = maspy.core.MsrunContainer()

    basename = os.path.basename(filepath)
    dirname = os.path.dirname(filepath)
    filename, extension = os.path.splitext(basename)

    #Check if the specified file is valid for an import
    if not os.path.isfile(filepath):
        raise IOError('File does not exist: %s' % filepath)
    elif extension.lower() != '.mzml':
        raise IOError('Filetype is not "mzml": %s' % filepath)
    elif filename in msrunContainer.info:
        print(filename, 'already present in the msrunContainer, aborting import.')
        return None

    specfilename = basename if specfilename is None else specfilename
    mzmlReader = maspy.xml.MzmlReader(filepath)
    masterContainer = {'rm': str(), 'ci': {}, 'si': {}, 'sai': {}, 'smi': {}}
    for xmlSpectrum in mzmlReader.parseSpectra():
        smi, binaryDataArrayList = smiFromXmlSpectrum(xmlSpectrum, basename)
        #Generate SpectrumItem
        si = maspy.core.Si(smi.id, smi.specfile)
        si.isValid = True
        siAttrFromSmi(smi, si)
        if si.msLevel > 1:
            si.precursorId = si.precursorId.split('scan=')[1] #TODO: change to use regex to extract from known vendor format
        #Generate SpectrumArrayItem
        sai = maspy.core.Sai(smi.id, smi.specfile)
        sai.arrays, sai.arrayInfo = maspy.xml.extractBinaries(binaryDataArrayList,
                                                         smi.attributes['defaultArrayLength'])
        #Store all items in the appropriate containers
        masterContainer['smi'][smi.id] = smi
        masterContainer['si'][smi.id] = si
        masterContainer['sai'][smi.id] = sai

    for xmlChromatogram in mzmlReader.chromatogramList:
        ci = ciFromXml(xmlChromatogram)
        masterContainer['ci'][ci.id] = ci
    masterContainer['rm'] = mzmlReader.metadataNode

    msrunContainer._addSpecfile(filename, dirname)
    msrunContainer.rmc[filename] = masterContainer['rm']
    msrunContainer.info[filename]['status']['rm'] = True
    msrunContainer.smic[filename] = masterContainer['smi']
    msrunContainer.info[filename]['status']['smi'] = True
    msrunContainer.sic[filename] = masterContainer['si']
    msrunContainer.info[filename]['status']['si'] = True
    msrunContainer.saic[filename] = masterContainer['sai']
    msrunContainer.info[filename]['status']['sai'] = True
    msrunContainer.cic[filename] = masterContainer['ci']
    msrunContainer.info[filename]['status']['ci'] = True

    return msrunContainer


def writeMzml(specfile, msrunContainer, outputdir, spectrumIds=None, chromatogramIds=None):
    #TODO: maybe change to use aux.openSafeReplace
    outputFile = io.BytesIO()

    #TODO: perform check that specfile is present in msrunContainer and at least the metadatanode.
    metadataTree = msrunContainer.rmc[specfile]
    #Generate a list of spectrum ids that should be written to mzML
    if spectrumIds is None and specfile in msrunContainer.smic:
        keyTuple = [(int(key), key) for key in viewkeys(msrunContainer.smic[specfile])]
        spectrumIds = [key for _, key in sorted(keyTuple)]
    spectrumCounts = len(spectrumIds)
    #Generate a list of chromatogram ids that should be written to mzML
    if chromatogramIds is None and specfile in msrunContainer.cic:
        chromatogramIds = [cId for cId in viewkeys(msrunContainer.cic[specfile])]
    chromatogramCounts = len(chromatogramIds)

    xmlFile = ETREE.xmlfile(outputFile, encoding='ISO-8859-1', buffered=False)
    xmlWriter = xmlFile.__enter__()
    xmlWriter.write_declaration()

    nsmap = {None: 'http://psi.hupo.org/ms/mzml', 'xsi': 'http://www.w3.org/2001/XMLSchema-instance'}
    xmlRoot = xmlWriter.element(metadataTree.tag, metadataTree.attrib, nsmap=nsmap)
    xmlRoot.__enter__()
    xmlWriter.write('\n')

    for metadataNode in metadataTree.getchildren():
        if metadataNode.tag == 'run':
            xmlRun = xmlWriter.element(metadataNode.tag, metadataNode.attrib)
            xmlRun.__enter__()
            xmlWriter.write('\n')
            for runChild in metadataNode.getchildren():
                if runChild.tag == 'spectrumList':
                    specDefaultProcRef = runChild.attrib['defaultDataProcessingRef']
                elif runChild.tag == 'chromatogramList':
                    chromDefaultProcRef = runChild.attrib['defaultDataProcessingRef']
                else:
                    #TODO: maybe recCopy?
                    xmlRun.append(runChild)

            #If any spectra should be written, generate the spectrumList Node
            if spectrumCounts > 0:
                specListAttribs = {'count': str(spectrumCounts),
                                   'defaultDataProcessingRef': specDefaultProcRef}
                xmlSpectrumList = xmlWriter.element('spectrumList', specListAttribs)
                xmlSpectrumList.__enter__()
                xmlWriter.write('\n')

                for index, key in enumerate(spectrumIds):
                    #TODO: proper container format instead of msrunContainer
                    smi = msrunContainer.smic[specfile][key]
                    sai = msrunContainer.saic[specfile][key]
                    xmlSpectrum = xmlSpectrumFromSmi(index, smi, sai)
                    xmlWriter.write(xmlSpectrum, pretty_print=True)
                xmlSpectrumList.__exit__(None, None, None)
                xmlWriter.write('\n')

            #If any chromatograms should be written, generate the chromatogramList Node
            if chromatogramCounts > 0:
                chromListAttribs = {'count': str(chromatogramCounts),
                                    'defaultDataProcessingRef': chromDefaultProcRef}
                xmlChromatogramList = xmlWriter.element('chromatogramList', chromListAttribs)
                xmlChromatogramList.__enter__()
                xmlWriter.write('\n')
                for index, key in enumerate(chromatogramIds):
                    #TODO: proper container format instead of msrunContainer
                    ci = msrunContainer.cic[specfile][key]
                    xmlChromatogram = xmlChromatogramFromCi(index, ci)
                    xmlWriter.write(xmlChromatogram, pretty_print=True)
                xmlChromatogramList.__exit__(None, None, None)
                xmlWriter.write('\n')

            xmlRun.__exit__(None, None, None)
            xmlWriter.write('\n')
        else:
            xmlWriter.write(maspy.xml.recCopyElement(metadataNode), pretty_print=True)
    # close
    xmlRoot.__exit__(None, None, None)
    xmlFile.__exit__(None, None, None)

    filepath = aux.joinpath(outputdir, specfile+'.mzML')
    with open(filepath, 'wb') as openfile:
        openfile.write(outputFile.getvalue())


##########################################################
# --- generate python objects from mzML xml elements --- #
##########################################################
def ciFromXml(xmlelement):
    ci = maspy.core.Ci()
    for key in ['id', 'dataProcessingRef']:
        if key in xmlelement.attrib:
            setattr(ci, key, xmlelement.attrib[key])
    ci.params, children = maspy.xml.extractParams(xmlelement)
    arrayLength = xmlelement.attrib['defaultArrayLength']

    for child in children:
        childElements, childParams = maspy.xml.sublistReader(child)
        childTag = maspy.xml.clearTag(child.tag)
        if childTag == 'precursor':
            #TODO: THIS HAS TO BE ADAPTED
            newElement = MzmlPrecursor(**element)
        elif childTag == 'product':
            #TODO: THIS HAS TO BE ADAPTED
            newElement = MzmlProduct(**element)
        elif childTag == 'binaryDataArrayList':
            for childElement in childElements:
                dataType, dataTypeParam = maspy.xml.findBinaryDataType(childElement['params'])
                ci.arrayInfo[dataType] = {'dataProcessingRef': None, 'params': childElement['params']}
                if 'dataProcessingRef' in childElement:
                    ci.arrayInfo[dataType]['dataProcessingRef'] = childElement['dataProcessingRef']
            ci.arrays, ci.arrayInfo = maspy.xml.extractBinaries(childElements, arrayLength)
    return ci


def smiFromXmlSpectrum(xmlelement, specfile):
    scanId = xmlelement.attrib['id'].split('scan=')[1] #TODO: change to use regex to extract from known vendor format
    smi = maspy.core.Smi(scanId, specfile)
    smi.params, spectrumChildren = maspy.xml.extractParams(xmlelement)
    smi.attributes.update(xmlelement.attrib)
    binaryDataArrayList = list()

    for spectrumChild in spectrumChildren:
        spectrumChildTag = maspy.xml.clearTag(spectrumChild.tag)
        elements, params = maspy.xml.sublistReader(spectrumChild)
        if params:
            #Define scanListParams here
            setattr(smi, spectrumChildTag+'Params', params)
        for element in elements:
            if spectrumChildTag == 'scanList':
                newElement = maspy.core.MzmlScan(**element)
            elif spectrumChildTag == 'precursorList':
                newElement = maspy.core.MzmlPrecursor(**element)
            elif spectrumChildTag == 'productList':
                newElement = maspy.core.MzmlProduct(**element)
            elif spectrumChildTag == 'binaryDataArrayList':
                binaryDataArrayList.append(element)
                continue
            else:
                raise Exception('unexpected spectrum xml child tag!')
            getattr(smi, spectrumChildTag).append(newElement)
    return smi, binaryDataArrayList


####################################################################################
# --- extract spectrum attributes from the mzml like structure of a Smi object --- #
####################################################################################
def fetchSpectrumInfo(smi):
    attributes = dict()
    for key, accession, dtype in [('msLevel', 'MS:1000511', int),
                                  ('tic', 'MS:1000285', float),
                                  ('basepeakMz', 'MS:1000504', float),
                                  ('basepeakI', 'MS:1000505', float)
                                  ]:
        param = maspy.xml.findParam(smi.params, accession)
        if param is not None:
            attributes[key] = dtype(param[1])
        else:
            #Or don't add key to attributes?
            attributes[key] = None
    return attributes


def fetchScanInfo(smi):
    if smi.scanList:
        attributes = dict()
        for key, accession, dtype in [('rt', 'MS:1000016', float),
                                      ('iit', 'MS:1000927', float),
                                      ]:
            param = maspy.xml.findParam(smi.scanList[0]['params'], accession)
            if param is not None:
                attributes[key] = dtype(param[1])
                #TODO: raise a warning if time is neither in minutes nor in seconds
                if key == 'rt' and param[2] == 'UO:0000031':
                    #convert retention time from minutes to seconds
                    attributes['rt'] *= 60
            else:
                attributes[key] = None
    else:
        attributes = None
    return attributes


def fetchParentIon(smi):
    if smi.precursorList:
        attributes = dict()
        attributes['precursorId'] = smi.precursorList[0].spectrumRef
        selectedIon = smi.precursorList[0]['selectedIonList'][0]
        for key, accession, dtype in [('mz', 'MS:1000744', float),
                                      ('i', 'MS:1000042', float),
                                      ('charge', 'MS:1000041', int)
                                      ]:
            param = maspy.xml.findParam(selectedIon, accession)
            if param is not None:
                attributes[key] = dtype(param[1])
            else:
                attributes[key] = None
    else:
        attributes = None
    return attributes


def defaultFetchSiAttrFromSmi(smi, si):
    """Default method to extract attributes from a spectrum metadata item (sai) and
    add them to a spectrum item (si)."""
    for key, value in viewitems(fetchSpectrumInfo(smi)):
        setattr(si, key, value)
    for key, value in viewitems(fetchScanInfo(smi)):
        setattr(si, key, value)
    if si.msLevel > 1:
        for key, value in viewitems(fetchParentIon(smi)):
            setattr(si, key, value)


#####################################################
# --- generate mzml elements from maspy objects --- #
#####################################################
def xmlGenScanList(scanList, scanListParams):
    numEntries = len(scanList)
    xmlScanList = ETREE.Element('scanList', {'count': str(numEntries)})
    maspy.xml.xmlAddParams(xmlScanList, scanListParams)
    for scan in scanList:
        #Note: no attributes supported
        xmlScan = ETREE.Element('scan', {})
        maspy.xml.xmlAddParams(xmlScan, scan['params'])

        #Generate the scanWindowList entry
        numScanWindows = len(scan['scanWindowList'])
        if numScanWindows > 0:
            xmlScanWindowList = ETREE.Element('scanWindowList', {'count': str(numScanWindows)})
            for scanWindow in scan['scanWindowList']:
                xmlScanWindow = ETREE.Element('scanWindow')
                maspy.xml.xmlAddParams(xmlScanWindow, scanWindow)
                xmlScanWindowList.append(xmlScanWindow)
            xmlScan.append(xmlScanWindowList)

        xmlScanList.append(xmlScan)
    return xmlScanList


def xmlGenPrecursorList(precursorList):
    numEntries = len(precursorList)
    xmlPrecursorList = ETREE.Element('precursorList', {'count': str(numEntries)})
    for precursor in precursorList:
        #Note: no attributes supported
        precursorAttrib = {}
        if precursor['spectrumRef'] is not None:
            precursorAttrib.update({'spectrumRef': precursor['spectrumRef']})
        xmlPrecursor = ETREE.Element('precursor', precursorAttrib)

        xmlActivation = ETREE.Element('activation')
        maspy.xml.xmlAddParams(xmlActivation, precursor['activation'])
        xmlPrecursor.append(xmlActivation)

        if precursor['isolationWindow'] is not None:
            xmlIsolationWindow = ETREE.Element('isolationWindow')
            maspy.xml.xmlAddParams(xmlIsolationWindow, precursor['isolationWindow'])
            xmlPrecursor.append(xmlIsolationWindow)

        #Generate the selectedIonList entry
        numSelectedIons = len(precursor['selectedIonList'])
        if numSelectedIons > 0:
            xmlSelectedIonList = ETREE.Element('selectedIonList', {'count': str(numSelectedIons)})
            for selectedIon in precursor['selectedIonList']:
                xmlSelectedIon = ETREE.Element('selectedIon')
                maspy.xml.xmlAddParams(xmlSelectedIon, selectedIon)
                xmlSelectedIonList.append(xmlSelectedIon)
            xmlPrecursor.append(xmlSelectedIonList)

        xmlPrecursorList.append(xmlPrecursor)
    return xmlPrecursorList


def xmlGenProductList(productList):
    #TODO
    pass


def xmlGenBinaryDataArrayList(binaryDataInfo, binaryDataDict, compression='zlib', arrayTypes=None):
    #compression: any other value than 'zlib' results in no compression
    #Use datatypes parameter to specify the order of the arrays
    arrayTypes = [_ for _ in viewkeys(binaryDataInfo)] if arrayTypes is None else arrayTypes
    numEntries = len(binaryDataInfo)
    xmlBinaryDataArrayList = ETREE.Element('binaryDataArrayList', {'count': str(numEntries)})
    for arrayType in arrayTypes:
        _, dataTypeParam = maspy.xml.findBinaryDataType(binaryDataInfo[arrayType]['params'])
        binaryData = binaryDataDict[arrayType]
        bitEncoding = '64' if binaryData.dtype.str == '<f8' else '32'
        binaryData, arrayLength = maspy.xml.encodeBinaryData(binaryData, bitEncoding, compression)

        # --- define binaryDataArray parameters --- #
        params = list()
        if bitEncoding == '64':
            params.append(('MS:1000523', None, None))
        else:
            params.append(('MS:1000521', None, None))
        if compression == 'zlib':
            params.append(('MS:1000574', None, None))
        else:
            params.append(('MS:1000576', None, None))
        mandatoryAccessions = ['MS:1000523', 'MS:1000521', 'MS:1000574', 'MS:1000576']
        for param in binaryDataInfo[arrayType]['params']:
            if param[0] not in mandatoryAccessions:
                params.append(param)

        #Note: not all attributes supported
        binaryDataArrayAttrib = {'encodedLength': str(len(binaryData))}
        for attr in ['dataProcessingRef']:
            if binaryDataInfo[arrayType][attr] is not None:
                binaryDataArrayAttrib[attr] = binaryDataInfo[arrayType][attr]
        xmlBinaryDataArray = ETREE.Element('binaryDataArray', binaryDataArrayAttrib)
        maspy.xml.xmlAddParams(xmlBinaryDataArray, params)

        xmlBinary = ETREE.Element('binary')
        xmlBinary.text = binaryData
        xmlBinaryDataArray.append(xmlBinary)
        xmlBinaryDataArrayList.append(xmlBinaryDataArray)
    return xmlBinaryDataArrayList


def xmlSpectrumFromSmi(index, smi, sai=None, compression='zlib'):
    """
    index: The zero-based, consecutive index of the spectrum in the SpectrumList. (mzML specification)
    smi: a SpectrumMetadataItem instance
    sai: a SpectrumArrayItem instance, if None no binaryDataArrayList is written
    """
    if sai is not None:
        arrayLength = [array.size for array in viewvalues(sai.arrays)]
        if len(set(arrayLength)) != 1:
            raise Exception('Unequal array size for different entries in sai.arrays')
        else:
            arrayLength = arrayLength[0]
    else:
        arrayLength = 0

    spectrumAttrib = {'index': str(index), 'id': smi.attributes['id'],
                      'defaultArrayLength': str(arrayLength)}

    xmlSpectrum = ETREE.Element('spectrum', **spectrumAttrib)
    maspy.xml.xmlAddParams(xmlSpectrum, smi.params)
    #Add the scanList
    if len(smi.scanList) > 0:
        xmlSpectrum.append(xmlGenScanList(smi.scanList, smi.scanListParams))
    if len(smi.precursorList) > 0:
        xmlSpectrum.append(xmlGenPrecursorList(smi.precursorList))
    if len(smi.productList) > 0:
        xmlSpectrum.append(xmlGenProductList(smi.productList))
    if sai is not None:
        xmlSpectrum.append(xmlGenBinaryDataArrayList(sai.arrayInfo,
                                                     sai.arrays,
                                                     compression=compression
                                                     ))
    return xmlSpectrum


def xmlChromatogramFromCi(index, ci, compression='zlib'):
    arrayLength = [array.size for array in viewvalues(ci.arrays)]
    if len(set(arrayLength)) != 1:
        raise Exception('Unequal array size for different entries in sai.arrays')
    else:
        arrayLength = arrayLength[0]

    chromatogramAttrib = {'index': str(index), 'id': ci.id,
                          'defaultArrayLength': str(arrayLength)}
    if 'dataProcessingRef' in ci.attrib:
        chromatogramAttrib.update({'dataProcessingRef': dataProcessingRef})

    xmlChromatogram = ETREE.Element('chromatogram', **chromatogramAttrib)
    maspy.xml.xmlAddParams(xmlChromatogram, ci.params)
    #TODO: add appropriate functions for precursor and product
    if ci.product is not None:
        pass
    if ci.precursor is not None:
        pass

    #Sort the array keys, that 'rt' is always the first, necessary for the software "SeeMS" to properly display chromatograms
    arrayTypes = set(ci.arrayInfo)
    if 'rt' in arrayTypes:
        arrayTypes.remove('rt')
        arrayTypes = ['rt'] + list(arrayTypes)
    else:
        arrayTypes = list(arrayTypes)

    xmlChromatogram.append(xmlGenBinaryDataArrayList(ci.arrayInfo,
                                                     ci.arrays,
                                                     compression=compression,
                                                     arrayTypes=arrayTypes
                                                     )
                           )
    return xmlChromatogram


# --- import for SiiContainer class --- #
def importPsmResults(siiContainer, fileLocation, specfile, psmType='percolator', psmEngine='comet', qValue=0.01):
    """Function to control the import of PSM results into :class:`SiiContainer`.

    :ivar siiContainer: Add PSMs to this instance
    :ivar fileLocation: Actual path to file
    :ivar specfile: Keyword (filename) to represent file in the :class:`SiiContainer`.
    :ivar psmType: can be used to specify post processing tools like percolator, used to choose import function
    :ivar psmEngine: specify peptide spectrum matching engine, used to choose import function
    :ivar qValue: define a qValue cut off for valid items, could be changed to (:var:`scoreCutOff` and :var:`scoreKey`)

    See also :func:`_importPercolatorResults` and :func:`_importFromPercolatorArray`
    """
    if specfile not in siiContainer.info:
        siiContainer.addSpecfile(specfile, os.path.dirname(fileLocation))

    if psmType == 'percolator':
        siiContainer.info[specfile]['scoreAttr'] = 'score'
        siiContainer.info[specfile]['largerBetter'] = True
        _psmArrays = _importPercolatorResults(fileLocation, psmEngine=psmEngine)
        _importFromPercolatorArrays(siiContainer, _psmArrays, specfile, qValueCutOff=qValue)


def _importPercolatorResults(fileLocation, psmEngine=None):
    """Reads percolator PSM results from a txt file.

    :ivar fileLocation: File path
    :ivar psmEngine: Specifies the used peptide spectrum matching engine ('comet', 'msgf', 'xtandem')

    :return: {attribute:numpy.array(), attribute:numpy.array()}

    See also :func:`importPsmResults` and :func:`_importFromPercolatorArray`
    """
    #HEADERLINE: xtandem seperates proteins with ';', msgf separates proteins by a tab
    with io.open(fileLocation, 'r', encoding='utf-8') as openfile:
        lines = openfile.readlines()
        headerDict = dict([[y,x] for (x,y) in enumerate(lines[0].strip().split('\t'))])
        scanEntryList = list()
        for line in lines[1:]:
            if len(line.strip()) == 0:
                continue
            fields = line.strip().split('\t')
            entryDict = dict()
            for headerName, headerPos in viewitems(headerDict):
                entryDict[headerName] = fields[headerPos]
            if psmEngine == 'msgf':
                entryDict['proteinIds'] = list(fields[headerDict['proteinIds']:])
            elif psmEngine == 'xtandem':
                entryDict['proteinIds'] = entryDict['proteinIds'].split(';')
            scanEntryList.append(entryDict)

    scanArrDict = dict()
    for headerName in viewkeys(headerDict):
        scanArrDict[headerName] = list()

    # Define list of headers #
    for scanEntryDict in scanEntryList:
        for headerName,entry in viewitems(scanEntryDict):
            if headerName in ['score','q-value','posterior_error_prob']:
                scanArrDict[headerName].append(float(entry))
            else:
                scanArrDict[headerName].append(entry)

    if psmEngine in ['comet','msgf','xtandem']:
        scanNrList = list()
        for entry in scanArrDict['PSMId']:
            if psmEngine in ['comet','msgf']:
                scanNr = entry.split('_')[-3]
            elif psmEngine in ['xtandem']:
                scanNr = entry.split('_')[-2]
            scanNrList.append( scanNr )
        scanArrDict['scanNr'] = scanNrList
    else:
        print('No valid psm engine specified, can\'t import percolator results!')

    for headerName in list(viewkeys(scanArrDict)):
        scanArrDict[headerName] = numpy.array(scanArrDict[headerName])
    return scanArrDict


def _importFromPercolatorArrays(siiContainer, psmArrays, specfile, qValueCutOff=None):
    """Writes :class:`SpectruIdentificationItem` into :class:`SiContainer`.

    :ivar siiContainer: Add PSMs to this instance
    :ivar psmArrays: contains PSM information, generated by :func:`_importFromPercolatorArray`
    :ivar specfile: Keyword (filename) to represent file in the :class:`SiContainer`
    :ivar qValueCutOff: define a qValue cut off for valid items

    See also :func:`importPsmResults`
    """
    scoreAttr = siiContainer.info[specfile]['scoreAttr']
    if siiContainer.info[specfile]['largerBetter']:
        sortMask = psmArrays[scoreAttr].argsort()[::-1]
    else:
        sortMask = psmArrays[scoreAttr].argsort()
    for key in psmArrays:
        psmArrays[key] = psmArrays[key][sortMask]

    for currPosition in range(0, len(psmArrays['scanNr'])):
        peptide = psmArrays['peptide'][currPosition]
        if peptide.find('.') != -1:
            peptide = peptide.split('.')[1]

        sequence = maspy.peptidemethods.removeModifications(peptide)
        scanNr = psmArrays['scanNr'][currPosition]
        qValue = psmArrays['q-value'][currPosition]
        score = psmArrays['score'][currPosition]
        psmId = psmArrays['PSMId'][currPosition]
        pep = psmArrays['posterior_error_prob'][currPosition]

        sii = maspy.core.Sii(scanNr, specfile)
        sii.peptide = peptide
        sii.sequence = sequence
        sii.qValue = qValue
        sii.score = score
        sii.pep = pep
        sii.isValid = False

        if sii.id in siiContainer.container[specfile]:
            siiList = siiContainer.container[specfile][sii.id]
            sii.rank = len(siiList) + 1
        else:
            sii.rank = 1
            siiContainer.container[specfile][sii.id] = list()

        if sii.rank == 1:
            if qValueCutOff is not None:
                if sii.qValue <= qValueCutOff:
                    sii.isValid = True
            else:
                sii.isValid = True

        siiContainer.container[specfile][sii.id].append(sii)


# --- import for FeatureContainer class --- #
def importPeptideFeatures(fiContainer, filelocation, specfile):
    """ Import peptide features from a featureXml file (eg. generated by OPENMS featureFinderCentroided).

    :param fiContainer: Spectra are added to to this instance of :class:`FeatureContainer`
    :param filelocation: Actual file path
    :param specfile: Keyword (filename) to represent file in the :class:`FeatureContainer`. Each filename
    can only occure once, therefore importing the same filename again is prevented.
    """
    if not os.path.isfile(filelocation):
        print('File does not exits:', filelocation)
    elif not filelocation.lower().endswith('.featurexml'):
        print('File is not a featurexml file:', filelocation)
    else:
        if specfile in fiContainer.info:
            print(specfile, 'is already present in the SiContainer, import interrupted.')
        else:
            fiContainer.addSpecfile(specfile, os.path.dirname(filelocation))
            featureDict = _importFeatureXml(filelocation)

            for featureId, featureEntryDict in viewitems(featureDict):
                rtArea = set()
                for convexHullEntry in featureEntryDict['convexHullDict']['0']:
                    rtArea.update([convexHullEntry[0]])

                fi = maspy.core.Fi(featureId, specfile)
                fi.rt = featureEntryDict['rt']
                fi.rtArea = max(rtArea) - min(rtArea)
                fi.rtLow = min(rtArea)
                fi.rtHigh = max(rtArea)
                fi.charge = featureEntryDict['charge']
                fi.mz = featureEntryDict['mz']
                fi.mh = maspy.peptidemethods.calcMhFromMz(featureEntryDict['mz'],
                                                          featureEntryDict['charge'])
                fi.intensity = featureEntryDict['intensity']
                fi.quality = featureEntryDict['overallquality']
                fi.isMatched = False
                fi.isAnnotated = False
                fi.isValid = True

                fiContainer.container[specfile][featureId] = fi


def _importFeatureXml(fileLocation):
    """Reads a featureXml file.

    :return: {featureKey1: {attribute1:value1, attribute2:value2, ...}, ...}

    See also :func:`importPeptideFeatures`
    """
    with io.open(fileLocation, 'r', encoding='utf-8') as openFile:
        readingFeature = False
        readingHull = False
        featureDict = dict()

        for i, line in enumerate(openFile):
            line = line.strip()
            if readingFeature == True:
                if line.find('<convexhull') != -1:
                    readingHull = True
                    hullNr = line.split('<convexhull nr=\"')[1].split('\">')[0]
                    hullList = list()
                elif readingHull == True:
                    if line.find('<pt') != -1:
                        x = float(line.split('x=\"')[1].split('\"')[0])
                        y = float(line.split('y=\"')[1].split('\"')[0])
                        # x = retentiontime, y = m/z
                        #retentionTimeList.append(x)
                        hullList.append([x,y])
                    elif line.find('</convexhull>') != -1:
                        featureDict[featureKey]['convexHullDict'][hullNr] = hullList
                        readingHull = False

                elif line.find('<position dim=\"0\">') != -1:
                    featureDict[featureKey]['dim0'] = float(line.split('<position dim=\"0\">')[1].split('</position>')[0])
                elif line.find('<position dim=\"1\">') != -1:
                    featureDict[featureKey]['dim1'] = float(line.split('<position dim=\"1\">')[1].split('</position>')[0])
                elif line.find('<intensity>') != -1:
                    featureDict[featureKey]['intensity'] = float(line.split('<intensity>')[1].split('</intensity>')[0])
                elif line.find('<overallquality>') != -1:
                    featureDict[featureKey]['overallquality'] = float(line.split('<overallquality>')[1].split('</overallquality>')[0])
                elif line.find('<charge>') != -1:
                    featureDict[featureKey]['charge'] = int( line.split('<charge>')[1].split('</charge>')[0] )

                elif line.find('<userParam') != -1:
                    if line.find('name=\"label\"') != -1:
                        featureDict[featureKey]['label'] = line.split('value=\"')[1].split('\"/>')[0]
                    elif line.find('name=\"score_fit\"') != -1:
                        featureDict[featureKey]['score_fit'] = float(line.split('value=\"')[1].split('\"/>')[0])
                    elif line.find('name=\"score_correlation\"') != -1:
                        featureDict[featureKey]['score_correlation'] = float(line.split('value=\"')[1].split('\"/>')[0])
                    elif line.find('name=\"FWHM\"') != -1:
                        featureDict[featureKey]['FWHM'] = float(line.split('value=\"')[1].split('\"/>')[0])
                    elif line.find('name=\"spectrum_index\"') != -1:
                        featureDict[featureKey]['spectrum_index'] = line.split('value=\"')[1].split('\"/>')[0]
                    elif line.find('name=\"spectrum_native_id\"') != -1:
                        featureDict[featureKey]['spectrum_native_id'] = line.split('value=\"')[1].split('\"/>')[0]

                elif line.find('</feature>') != -1:
                    #mzList = list()
                    #for retentionTime,mz in featureDict[featureKey]['convexHullDict']['0']:
                    #    mzList.append(mz)
                    featureDict[featureKey]['rt'] = featureDict[featureKey]['dim0']#numpy.median(retentionTimeList)
                    featureDict[featureKey]['mz'] = featureDict[featureKey]['dim1']#numpy.median(mzList)

                    readingFeature == False

            if line.find('<feature id') != -1:
                readingFeature = True
                featureKey = line.split('<feature id=\"')[1].split('\">')[0]
                featureDict[featureKey] = dict()
                featureDict[featureKey]['convexHullDict'] = dict()
                #retentionTimeList = list()
    return featureDict
