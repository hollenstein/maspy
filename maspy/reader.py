from __future__ import print_function, division, unicode_literals
from future.utils import viewkeys, viewvalues, viewitems, listvalues, listitems

try: # python 2.7
    from itertools import izip as zip
except ImportError: # python 3.x series
    pass
###############################################################################
import io
from lxml import etree as ETREE
import numpy
import os

import maspy.auxiliary as aux
import maspy.core
import maspy.peptidemethods
import maspy.xml


#####################################################################
### mzml import methods                       #######################
#####################################################################
def importMzml(filepath, msrunContainer=None, siAttrFromSmi=None, specfilename=None):
    """Performs a complete import of a mzml file into a maspy MsrunContainer.

    siAttrFromSmi: allow here to specify a custom function that extracts params a from spectrumMetadataItem
    """
    #TODO: docstring
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


# --- generate python objects from mzML xml elements --- #
def ciFromXml(xmlelement):
    #TODO: docstring
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
    #TODO: docstring
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


# --- extract spectrum attributes from the mzml like structure of a Smi object --- #
def fetchSpectrumInfo(smi):
    #TODO: docstring
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
    #TODO: docstring
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
    #TODO: docstring
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
    adding them to a spectrum item (si)."""
    for key, value in viewitems(fetchSpectrumInfo(smi)):
        setattr(si, key, value)
    for key, value in viewitems(fetchScanInfo(smi)):
        setattr(si, key, value)
    if si.msLevel > 1:
        for key, value in viewitems(fetchParentIon(smi)):
            setattr(si, key, value)


#####################################################################
### import for SiiContainer class             #######################
#####################################################################
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
    #Note: regarding headerline, xtandem seperates proteins with ';', msgf separates proteins by a tab
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


#####################################################################
### import for FeatureContainer class         #######################
#####################################################################
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
