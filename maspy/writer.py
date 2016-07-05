from __future__ import print_function, division, unicode_literals
from future.utils import viewkeys, viewvalues, viewitems, listvalues, listitems

try: # python 2.7
    from itertools import izip as zip
except ImportError: # python 3.x series
    pass
###############################################################################
import hashlib
import io
from lxml import etree as ETREE
import numpy
import os

import maspy.auxiliary as aux
import maspy.xml


##########################################################
### mzml import and export methods #######################
##########################################################
def writeMzml(specfile, msrunContainer, outputdir, spectrumIds=None,
              chromatogramIds=None, writeIndex=True):
    """ #TODO: docstring

    :param specfile: #TODO docstring
    :param msrunContainer: #TODO docstring
    :param outputdir: #TODO docstring
    :param spectrumIds: #TODO docstring
    :param chromatogramIds: #TODO docstring
    """
    #TODO: maybe change to use aux.openSafeReplace
    outputFile = io.BytesIO()

    #TODO: perform check that specfile is present in msrunContainer and at least
    #   the metadatanode.
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

    spectrumIndexList = list()
    chromatogramIndexList = list()

    xmlFile = ETREE.xmlfile(outputFile, encoding='ISO-8859-1', buffered=False)
    xmlWriter = xmlFile.__enter__()
    xmlWriter.write_declaration()

    nsmap = {None: 'http://psi.hupo.org/ms/mzml',
             'xsi': 'http://www.w3.org/2001/XMLSchema-instance'
             }
    mzmlAttrib = {'{http://www.w3.org/2001/XMLSchema-instance}schemaLocation': \
                    'http://psi.hupo.org/ms/mzml http://psidev.info/files/ms/mzML/xsd/mzML1.1.0.xsd',
                  'version': '1.1.0', 'id': metadataTree.attrib['id']
                  }

    if writeIndex:
        xmlIndexedMzml = xmlWriter.element('indexedmzML', nsmap=nsmap)
        xmlIndexedMzml.__enter__()
        xmlWriter.write('\n')
    xmlMzml = xmlWriter.element('mzML', mzmlAttrib, nsmap=nsmap)
    xmlMzml.__enter__()
    xmlWriter.write('\n')

    for metadataNode in metadataTree.getchildren():
        if metadataNode.tag != 'run':
            xmlWriter.write(maspy.xml.recCopyElement(metadataNode),
                            pretty_print=True
                            )
        else:
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

            #If any spectra should be written, generate the spectrumList Node.
            if spectrumCounts > 0:
                specListAttribs = {'count': str(spectrumCounts),
                                   'defaultDataProcessingRef': specDefaultProcRef
                                   }
                xmlSpectrumList = xmlWriter.element('spectrumList',
                                                    specListAttribs
                                                    )
                xmlSpectrumList.__enter__()
                xmlWriter.write('\n')

                for index, key in enumerate(spectrumIds):
                    smi = msrunContainer.smic[specfile][key]
                    sai = msrunContainer.saic[specfile][key]
                    #Store the spectrum element offset here
                    spectrumIndexList.append((outputFile.tell(),
                                              smi.attributes['id']
                                              ))

                    xmlSpectrum = xmlSpectrumFromSmi(index, smi, sai)
                    xmlWriter.write(xmlSpectrum, pretty_print=True)

                xmlSpectrumList.__exit__(None, None, None)
                xmlWriter.write('\n')

            #If any chromatograms should be written, generate the
            #chromatogramList Node.
            if chromatogramCounts > 0:
                chromListAttribs = {'count': str(chromatogramCounts),
                                    'defaultDataProcessingRef': chromDefaultProcRef
                                    }
                xmlChromatogramList = xmlWriter.element('chromatogramList',
                                                        chromListAttribs
                                                        )
                xmlChromatogramList.__enter__()
                xmlWriter.write('\n')
                for index, key in enumerate(chromatogramIds):
                    ci = msrunContainer.cic[specfile][key]
                    #Store the chromatogram element offset here
                    chromatogramIndexList.append((outputFile.tell(), ci.id))

                    xmlChromatogram = xmlChromatogramFromCi(index, ci)
                    xmlWriter.write(xmlChromatogram, pretty_print=True)
                xmlChromatogramList.__exit__(None, None, None)
                xmlWriter.write('\n')

            xmlRun.__exit__(None, None, None)
            xmlWriter.write('\n')

    #Close the mzml node
    xmlMzml.__exit__(None, None, None)
    #Optional: write the indexedMzml nodes and close the indexedMzml node
    if writeIndex:
        xmlWriter.write('\n')
        indexListOffset = outputFile.tell()
        _writeMzmlIndexList(xmlWriter, spectrumIndexList, chromatogramIndexList)
        _writeIndexListOffset(xmlWriter, indexListOffset)
        _writeMzmlChecksum(xmlWriter, outputFile)
        xmlIndexedMzml.__exit__(None, None, None)
    #Close the xml file
    xmlFile.__exit__(None, None, None)
    #Write the output mzML file
    filepath = aux.joinpath(outputdir, specfile+'.mzML')
    with open(filepath, 'wb') as openfile:
        openfile.write(outputFile.getvalue())


def _writeMzmlIndexList(xmlWriter, spectrumIndexList, chromatogramIndexList):
    """ #TODO: docstring

    :param xmlWriter: #TODO: docstring
    :param spectrumIndexList: #TODO: docstring
    :param chromatogramIndexList: #TODO: docstring
    """
    counts = 0
    if spectrumIndexList:
        counts += 1
    if chromatogramIndexList:
        counts += 1
    if counts == 0:
        return None
    #Create indexList node
    xmlIndexList = xmlWriter.element('indexList', {'count': str(counts)})
    xmlIndexList.__enter__()
    xmlWriter.write('\n')

    _writeIndexListElement(xmlWriter, 'spectrum', spectrumIndexList)
    _writeIndexListElement(xmlWriter, 'chromatogram', chromatogramIndexList)

    #Close indexList node
    xmlIndexList.__exit__(None, None, None)
    xmlWriter.write('\n')


def _writeIndexListElement(xmlWriter, elementName, indexList):
    """ #TODO: docstring

    :param xmlWriter: #TODO: docstring
    :param elementName: #TODO: docstring
    :param indexList: #TODO: docstring
    """
    if indexList:
        xmlIndex = xmlWriter.element('index', {'name': elementName})
        xmlIndex.__enter__()
        xmlWriter.write('\n')
        for offset, indexId in indexList:
            offsetElement = ETREE.Element('offset', {'idRef': indexId})
            offsetElement.text = str(offset)
            xmlWriter.write(offsetElement, pretty_print=True)
        xmlIndex.__exit__(None, None, None)
        xmlWriter.write('\n')


def _writeMzmlChecksum(xmlWriter, outputFile):
    """ #TODO: docstring

    :param xmlWriter: #TODO: docstring
    :param outputFile: #TODO: docstring
    """
    sha = hashlib.sha1(outputFile.getvalue())
    sha.update('<fileChecksum>')

    xmlChecksumElement = ETREE.Element('fileChecksum')
    xmlChecksumElement.text = sha.hexdigest()

    xmlWriter.write(xmlChecksumElement, pretty_print=True)


def _writeIndexListOffset(xmlWriter, offset):
    """ #TODO: docstring

    :param xmlWriter: #TODO: docstring
    :param offset: #TODO: docstring
    """
    xmlIndexListOffset = ETREE.Element('indexListOffset')
    xmlIndexListOffset.text = str(offset)

    xmlWriter.write(xmlIndexListOffset, pretty_print=True)


# --- generate mzml elements from maspy objects --- #
def xmlGenScanList(scanList, scanListParams):
    """ #TODO: docstring

    :params scanList: #TODO: docstring
    :params scanListParams: #TODO: docstring

    :returns: #TODO: docstring
    """
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
            xmlScanWindowList = ETREE.Element('scanWindowList',
                                              {'count': str(numScanWindows)}
                                              )
            for scanWindow in scan['scanWindowList']:
                xmlScanWindow = ETREE.Element('scanWindow')
                maspy.xml.xmlAddParams(xmlScanWindow, scanWindow)
                xmlScanWindowList.append(xmlScanWindow)
            xmlScan.append(xmlScanWindowList)

        xmlScanList.append(xmlScan)
    return xmlScanList


def xmlGenPrecursorList(precursorList):
    """ #TODO: docstring

    :params precursorList: #TODO: docstring

    :returns: #TODO: docstring
    """
    numEntries = len(precursorList)
    xmlPrecursorList = ETREE.Element('precursorList',
                                     {'count': str(numEntries)}
                                     )
    for precursor in precursorList:
        #Note: no attributes for external referencing supported
        precursorAttrib = {}
        if precursor['spectrumRef'] is not None:
            precursorAttrib.update({'spectrumRef': precursor['spectrumRef']})
        xmlPrecursor = ETREE.Element('precursor', precursorAttrib)

        #Add isolationWindow element
        if precursor['isolationWindow'] is not None:
            xmlIsolationWindow = ETREE.Element('isolationWindow')
            maspy.xml.xmlAddParams(xmlIsolationWindow,
                                   precursor['isolationWindow']
                                   )
            xmlPrecursor.append(xmlIsolationWindow)

        #Add selectedIonList element
        numSelectedIons = len(precursor['selectedIonList'])
        if numSelectedIons > 0:
            xmlSelectedIonList = ETREE.Element('selectedIonList',
                                               {'count': str(numSelectedIons)}
                                               )
            for selectedIon in precursor['selectedIonList']:
                xmlSelectedIon = ETREE.Element('selectedIon')
                maspy.xml.xmlAddParams(xmlSelectedIon, selectedIon)
                xmlSelectedIonList.append(xmlSelectedIon)
            xmlPrecursor.append(xmlSelectedIonList)

        #Add activation element
        xmlActivation = ETREE.Element('activation')
        maspy.xml.xmlAddParams(xmlActivation, precursor['activation'])
        xmlPrecursor.append(xmlActivation)


        xmlPrecursorList.append(xmlPrecursor)
    return xmlPrecursorList


def xmlGenProductList(productList):
    """ #TODO: docstring

    :params productList: #TODO: docstring

    :returns: #TODO: docstring
    """
    raise NotImplementedError('xmlGenProductList() is not yet implemented')


def xmlGenBinaryDataArrayList(binaryDataInfo, binaryDataDict,
                              compression='zlib', arrayTypes=None):
    """ #TODO: docstring

    :params binaryDataInfo: #TODO: docstring
    :params binaryDataDict: #TODO: docstring
    :params compression: #TODO: docstring
    :params arrayTypes: #TODO: docstring

    :returns: #TODO: docstring
    """
    #Note: any other value for "compression" than "zlib" results in no
    #   compression
    #Note: Use arrayTypes parameter to specify the order of the arrays
    if arrayTypes is None:
        arrayTypes = [_ for _ in viewkeys(binaryDataInfo)]
    numEntries = len(binaryDataInfo)
    xmlBinaryDataArrayList = ETREE.Element('binaryDataArrayList',
                                           {'count': str(numEntries)}
                                           )
    for arrayType in arrayTypes:
        _, dataTypeParam = maspy.xml.findBinaryDataType(binaryDataInfo[arrayType]['params'])
        binaryData = binaryDataDict[arrayType]
        bitEncoding = '64' if binaryData.dtype.str == '<f8' else '32'
        binaryData, arrayLength = maspy.xml.encodeBinaryData(binaryData,
                                                             bitEncoding,
                                                             compression
                                                             )

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
        mandatoryAccessions = ['MS:1000523', 'MS:1000521', 'MS:1000574',
                               'MS:1000576'
                               ]
        for param in binaryDataInfo[arrayType]['params']:
            if param[0] not in mandatoryAccessions:
                params.append(param)

        #Note: not all attributes supported
        binaryDataArrayAttrib = {'encodedLength': str(len(binaryData))}
        for attr in ['dataProcessingRef']:
            if binaryDataInfo[arrayType][attr] is not None:
                binaryDataArrayAttrib[attr] = binaryDataInfo[arrayType][attr]
        xmlBinaryDataArray = ETREE.Element('binaryDataArray',
                                           binaryDataArrayAttrib
                                           )
        maspy.xml.xmlAddParams(xmlBinaryDataArray, params)

        xmlBinary = ETREE.Element('binary')
        xmlBinary.text = binaryData
        xmlBinaryDataArray.append(xmlBinary)
        xmlBinaryDataArrayList.append(xmlBinaryDataArray)
    return xmlBinaryDataArrayList


def xmlSpectrumFromSmi(index, smi, sai=None, compression='zlib'):
    """ #TODO: docstring

    :param index: The zero-based, consecutive index of the spectrum in the
        SpectrumList. (mzML specification)
    :param smi: a SpectrumMetadataItem instance
    :param sai: a SpectrumArrayItem instance, if none is specified no
        binaryDataArrayList is written
    :param compression: #TODO: docstring

    :returns: #TODO: docstring
    """
    if sai is not None:
        arrayLength = [array.size for array in viewvalues(sai.arrays)]
        if len(set(arrayLength)) != 1:
            raise Exception('Unequal size for different array in sai.arrays')
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
    """ #TODO: docstring
    :param index: #TODO: docstring
    :param ci: #TODO: docstring
    :param compression: #TODO: docstring

    :returns: #TODO: docstring
    """
    arrayLength = [array.size for array in viewvalues(ci.arrays)]
    if len(set(arrayLength)) != 1:
        raise Exception('Unequal size for different array in sai.arrays')
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
        raise NotImplementedError()
    if ci.precursor is not None:
        raise NotImplementedError()

    #Sort the array keys, that 'rt' is always the first, necessary for example
    #   for the software "SeeMS" to properly display chromatograms.
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
