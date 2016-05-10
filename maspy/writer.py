from __future__ import print_function, division, unicode_literals
from future.utils import viewkeys, viewvalues, viewitems, listvalues, listitems

try: # python 2.7
    from itertools import izip as zip
except ImportError: # python 3.x series
    pass
###############################################################################
from lxml import etree as ETREE
import io
import numpy
import os

import maspy.auxiliary as aux
import maspy.xml


##########################################################
### mzml import and export methods #######################
##########################################################
def writeMzml(specfile, msrunContainer, outputdir, spectrumIds=None, chromatogramIds=None):
    #TODO: docstring
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


# --- generate mzml elements from maspy objects --- #
def xmlGenScanList(scanList, scanListParams):
    #TODO: docstring

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
    #TODO: docstring

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
    #TODO: docstring
    raise NotImplementedError('xmlGenProductList is not yet implemented')


def xmlGenBinaryDataArrayList(binaryDataInfo, binaryDataDict, compression='zlib', arrayTypes=None):
    #TODO: docstring

    #Note: any other value for "compression" than "zlib" results in no compression
    #Note: Use arrayTypes parameter to specify the order of the arrays
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
    #TODO: docstring
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
        raise NotImplementedError()
    if ci.precursor is not None:
        raise NotImplementedError()

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
