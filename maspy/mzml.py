from __future__ import print_function, division
from future.utils import viewkeys, viewvalues, viewitems, listvalues, listitems

import zlib
import numpy
from base64 import b64decode as b64dec
from base64 import b64encode as b64enc
from struct import unpack as unpack
from struct import pack as pack

from lxml import etree as etree
import maspy.auxiliary as aux


class FastReader(object):
    def __init__(self, inputPath):
        self.info = dict()
        self.info['inputPath'] = inputPath

        if self.info['inputPath'].endswith('.gz'):
            import gzip
            import codecs
            self.info['fileObject'] = codecs.getreader("utf-8")(
                gzip.open(self.info['inputPath'])
            )
        else:
            self.info['fileObject'] = open(self.info['inputPath'], 'rb')

        self.iter = iter(etree.iterparse(self.info['fileObject'], events=(b'start', b'end')))
        self.event, self.element = next(self)

    def __iter__(self):
        return self

    def __next__(self):
        """ The python 2.6+ iterator """
        return next(self)

    def next(self):
        event, element = next(self.iter, ('END', 'END'))

        # stop iteration and clear parsed elements when parsing is done
        if event == 'END':
            self.clearParsedElements()
            raise StopIteration()
        else:
            self.event, self.element = event, element

        return self.event, self.element

    def clearParsedElements(self):
        self.element.clear()
        while self.element.getprevious() is not None:
            del self.element.getparent()[0]


class PymzmlWriter(object):
    def __init__(self, inputPath, outputPath, debug=False):
        self.debug = debug

        self.info = dict()
        self.info['inputPath'] = inputPath
        self.info['outputPath'] = outputPath

        self.info['schemeLocation'] = dict()
        self.info['schemeLocation']['1.1.0'] = 'http://psidev.cvs.sourceforge.net/viewvc/psidev/psi/psi-ms/mzML/schema/mzML1.1.0.xsd'
        self.info['schemeLocation']['1.1.1'] = 'http://psidev.cvs.sourceforge.net/viewvc/psidev/psi/psi-ms/mzML/schema/mzML1.1.1_idx.xsd'

        self.info['nameSpaceMap'] = dict()
        self.info['nameSpaceMap'][None] = 'http://psi.hupo.org/ms/mzml'
        self.info['nameSpaceMap']['xsi'] = 'http://www.w3.org/2001/XMLSchema-instance'

        # Initialize mzML reader, if the input file is an indexed mzML file, go to mzML node
        self.reader = FastReader(self.info['inputPath'])
        if self.reader.event == b'start' and self.reader.element.tag.endswith('indexedmzML'):
            next(self.reader)

        # Open output file, write declaration and the XML root node, "mzML"
        if self.debug:
            import StringIO
            self.outputFile = StringIO.StringIO()
        else:
            self.outputFile = open(self.info['outputPath'], 'wb')

        self.xmlFile = etree.xmlfile(self.outputFile, encoding='ISO-8859-1', buffered=False)
        self.xmlWriter = self.xmlFile.__enter__()
        self.xmlWriter.write_declaration()

        self.xmlRoot = self.xmlWriter.element('mzML', nsmap=self.info['nameSpaceMap'],
                                              attrib={'id':self.reader.element.attrib['id'],
                                                      'version':'1.1.0',
                                                      '{http://www.w3.org/2001/XMLSchema-instance}schemaLocation':self.info['schemeLocation']['1.1.0']
                                                      }
                                              )
        self.xmlRoot.__enter__()
        self.xmlWriter.write('\n')

    def writeMetaData(self, newMetaNodeElements=dict()):
        next(self.reader)
        currSubChildTag = self.reader.element.tag
        # possible child nodes: cvList, fileDescription, referenceableParamGroupList, sampleList, softwareList, instrumentConfigurationList, dataProcessingList, run

        while not currSubChildTag.endswith('run'):
            while not(self.reader.element.tag == currSubChildTag and self.reader.event == b'end'):
                next(self.reader)

            # new child nodes can be added via the dictionary newMetaNodeElements #
            for newNodeTag, newNodeList in viewitems(newMetaNodeElements):
                if self.reader.element.tag.endswith(newNodeTag):
                    for newNode in newNodeList:
                        self.reader.element.append(newNode)
                    self.reader.element.attrib['count'] = str(len(self.reader.element.getchildren()))

            removeTreeFormating(self.reader.element)
            self.xmlWriter.write(self.reader.element, pretty_print=True)
            self.reader.clearParsedElements()

            # proceed to next child node #
            next(self.reader)
            currSubChildTag = self.reader.element.tag

    def writeRunNode(self, specWriter=None):
        xmlRunNode = self.xmlWriter.element('run',
                                            nsmap=self.info['nameSpaceMap'],
                                            attrib=self.reader.element.attrib
                                            )
        xmlRunNode.__enter__()
        self.xmlWriter.write('\n')

        next(self.reader)
        currSubChildTag = self.reader.element.tag
        while not (self.reader.element.tag.endswith('run') and self.reader.event == b'end'):
            if currSubChildTag.endswith('spectrumList'):
                ## WRITE SPECTRA HERE ##
                if specWriter is None:
                    self._writeAllSpectra()
                else:
                    specWriter()
                while not(self.reader.element.tag == currSubChildTag and self.reader.event == b'end'):
                    next(self.reader)

            else:
                while not(self.reader.element.tag == currSubChildTag and self.reader.event == b'end'):
                    next(self.reader)
                removeTreeFormating(self.reader.element)
                self.xmlWriter.write(self.reader.element, pretty_print=True)
                self.reader.clearParsedElements()

            # proceed to next child node #
            next(self.reader
            currSubChildTag = self.reader.element.tag

        xmlRunNode.__exit__(None,None,None)
        self.xmlWriter.write('\n')

    def _writeAllSpectra(self):
        xmlSpectrumListNode = self.xmlWriter.element('spectrumList',
                                                     nsmap=self.info['nameSpaceMap'],
                                                     attrib=self.reader.element.attrib
                                                     )
        xmlSpectrumListNode.__enter__()
        self.xmlWriter.write('\n')
        while not(self.reader.element.tag.endswith('spectrumList') and self.reader.event == b'end'):
            if self.reader.element.tag.endswith('spectrum') and self.reader.event == b'start':
                while not(self.reader.element.tag.endswith('spectrum') and self.reader.event == b'end'):
                    next(self.reader)
                removeTreeFormating(self.reader.element)
                self.xmlWriter.write(self.reader.element, pretty_print=True)
                self.reader.clearParsedElements()
            next(self.reader)
        xmlSpectrumListNode.__exit__(None, None, None)

    def cleanUp(self):
        self.xmlRoot.__exit__(None, None, None)
        self.xmlFile.__exit__(None, None, None)
        if not self.debug:
            self.outputFile.close()


def removeTreeFormating(element):
    for children in element.iter():
        if children.text is not None:
            if len(children.text.strip()) == 0:
                children.text = None
            else:
                children.text = children.text.strip()
        if children.tail is not None:
            if len(children.tail.strip()) == 0:
                children.tail = None
            else:
                children.tail = children.tail.strip()


def decodeBinaryData(binaryData, arrayLength, bitEncoding, compression):
    bitEncodedData = binaryData.encode("utf-8")
    bitDecodedData = b64dec(bitEncodedData)

    if bitEncoding == '64':
        floattype = 'd' # 64-bit
    else:
        floattype = 'f' # 32-bit

    if compression == 'zlib':
        decompressedData = zlib.decompress(bitDecodedData)
    else:
        decompressedData = bitDecodedData

    fmt = '{endian}{arraylength}{floattype}'.format(endian='<', arraylength=arrayLength, floattype=floattype)
    dataArray = numpy.array(unpack(fmt, decompressedData))
    return dataArray


def encodeBinaryData(dataArray, bitEncoding, compression):
    arrayLength = len(dataArray)
    if bitEncoding == '64':
        floattype = 'd' # 64-bit
    else:
        floattype = 'f' # 32-bit
    fmt = '{endian}{arraylength}{floattype}'.format(endian='<', arraylength=arrayLength, floattype=floattype)
    packedData = pack(fmt, *dataArray)

    if compression == 'zlib':
        compressedData = zlib.compress(packedData)
    else:
        compressedData = packedData

    encodedData = b64enc(compressedData)
    return encodedData, arrayLength


def writeXmlSpectrumList(writer, siContainer, selectedMsLevel=[1, 2], CorrFactor=None, CorrFactorMz=None,
                         ms2CorrFactorRt=None, ms2CorrFactorMz=None, ms2CorrectionType=None, scanSelector=None
                         ):
    # Adapted to pyms
    rtCorrection = True if CorrFactor is not None else False
    mzCorrection = True if CorrFactorMz is not None else False
    ms2RtCorrection = True if ms2CorrFactorRt is not None else False
    ms2MzCorrection = True if ms2CorrFactorMz is not None else False

    # Define which and how many spectra to write #
    # all spectra, only ms1, only ms2; #
    # ms2 can either be native or add mixed spectra #
    siArray = siContainer.getArrays(['rt'], filterAttribute='msLevel', selector=lambda x: x in aux.toList(selectedMsLevel))
    specfile = siContainer.specfiles[0]

    calibrateMs1 = True if (rtCorrection or mzCorrection) else False
    calibrateMs2 = True if (ms2RtCorrection or ms2MzCorrection) else False
    addMixedMs2 = False

    if scanSelector is None:
        scanSelector = lambda x: True

    spectraCounts = int()
    for scanNr in siArray['id']:
        if scanSelector(scanNr):
            spectraCounts += 1

    if addMixedMs2:
        pass

    while not(writer.reader.element.tag.endswith('spectrumList') and writer.reader.event == b'start'):
        next(writer.reader)
    writer.reader.element.attrib['count'] = str(spectraCounts)

    # Write and open spectrum list node
    xmlSpectrumListNode = writer.xmlWriter.element('spectrumList',
                                                   nsmap=writer.info['nameSpaceMap'],
                                                   attrib=writer.reader.element.attrib
                                                   )
    xmlSpectrumListNode.__enter__()
    writer.xmlWriter.write('\n')

    # Iterate over spectrum list and write spectra, if necessary correct ion m/z values
    spectrumIndexCounter = int()
    while not(writer.reader.element.tag.endswith('spectrumList') and writer.reader.event == b'end'):
        if writer.reader.element.tag.endswith('spectrum') and writer.reader.event == b'start':
            while not(writer.reader.element.tag.endswith('spectrum') and writer.reader.event == b'end'):
                next(writer.reader)

            nativeSpecId = writer.reader.element.attrib['id']
            scanNr = nativeSpecId.split('scan=')[1]
            msLevel = siContainer[(specfile, scanNr)].msLevel
            rt = siContainer[(specfile, scanNr)].rt
            precursorMzNode = None
            defaultArrayLength = writer.reader.element.attrib['defaultArrayLength']

            if msLevel not in selectedMsLevel:
                continue

            if not scanSelector(scanNr):
                continue

            removeTreeFormating(writer.reader.element)

            for spectrumChild in writer.reader.element.getchildren():
                if spectrumChild.tag.endswith('precursorList'):
                    for precursorListChild in spectrumChild.getchildren():
                        for precursorChild in precursorListChild.getchildren():
                            if precursorChild.tag.endswith('selectedIonList'):
                                for selectedIonListChild in precursorChild.getchildren():
                                    if selectedIonListChild.tag.endswith('selectedIon'):
                                        for selectedIonChild in selectedIonListChild.getchildren():
                                            if selectedIonChild.tag.endswith('cvParam'):
                                                if selectedIonChild.attrib['accession'] == 'MS:1000744':
                                                    precursorMzNode = selectedIonChild


                elif spectrumChild.tag.endswith('binaryDataArrayList'):
                    if (calibrateMs1 and msLevel == 1) or (calibrateMs2 and msLevel == 2):
                        for dataArrayListChild in spectrumChild.getchildren():
                            if dataArrayListChild.tag.endswith('binaryDataArray'):
                                spectrumCvParamList = list()
                                for dataArrayChild in dataArrayListChild.getchildren():
                                    if dataArrayChild.tag.endswith('cvParam'):
                                        spectrumCvParamList.append(dataArrayChild.attrib['accession'])
                                    elif dataArrayChild.tag.endswith('binary'):
                                        binaryData = dataArrayChild.text
                                if 'MS:1000514' in spectrumCvParamList: # is m/z array
                                    bitEncoding = '64' if 'MS:1000523' in spectrumCvParamList else '32'
                                    if 'MS:1000574' in spectrumCvParamList:
                                        compression = 'zlib'
                                    else:
                                        compression = None
                                    spectrumArray = decodeBinaryData(binaryData, defaultArrayLength, bitEncoding, compression)

                                    #########################
                                    # correct spectrum ions #
                                    #########################
                                    if msLevel == 1:
                                        if rtCorrection:
                                            spectrumArray = spectrumArray * (1 + (CorrFactor[rt] * 1e-6))
                                        if mzCorrection:
                                            corrArray = CorrFactorMz.corrArray(spectrumArray)
                                            spectrumArray = spectrumArray * (1 + corrArray * 1e-6)
                                    elif msLevel ==2:
                                        if ms2RtCorrection:
                                            if ms2CorrectionType == 'abs' or ms2CorrectionType == 'absolute':
                                                spectrumArray = spectrumArray  - ms2CorrFactorRt[rt]
                                            else:
                                                spectrumArray = spectrumArray * (1 + (ms2CorrFactorRt[rt] * 1e-6))

                                        if ms2MzCorrection:
                                            corrArray = ms2CorrFactorMz.corrArray(spectrumArray)
                                            if ms2CorrectionType == 'abs' or ms2CorrectionType == 'absolute':
                                                spectrumArray = spectrumArray - corrArray
                                            else:
                                                spectrumArray = spectrumArray * (1 + corrArray * 1e-6)

                                    binaryData, arrayLength = encodeBinaryData(spectrumArray, bitEncoding, compression)
                                    dataArrayChild.text = binaryData
                                    dataArrayListChild.attrib['encodedLength'] = str(len(binaryData))
                                    writer.reader.element.attrib['defaultArrayLength'] = str(arrayLength)
                                    # also add: base peak m/z 'MS:1000504', lowest observed m/z 'MS:1000528', highest observed m/z 'MS:1000527'
            if addMixedMs2:
                for mixedMs2 in mixedMs2Spectra[nativeSpecId]:
                    newSpecID = '1' # TODO: GENERATE THIS #
                    writer.reader.element.attrib['id'] = newSpecID
                    writer.reader.element.attrib['index'] = str(spectrumIndexCounter)
                    spectrumIndexCounter += 1

                    writer.xmlWriter.write(writer.reader.element, pretty_print=True)
                    writer.xmlWriter.flush()
                writer.reader.clearParsedElements()
            else:
                writer.reader.element.attrib['index'] = str(spectrumIndexCounter)
                spectrumIndexCounter += 1
                if msLevel == 2:
                    containerIdMs1 = siContainer[(specfile, scanNr)].ms1Id
                    rtMs1 = siContainer[containerIdMs1].rt
                    precursorMz = float(precursorMzNode.attrib['value'])
                    ## correct ms2 precursor here! ##
                    if rtCorrection:
                        precursorMz = precursorMz * (1 + (CorrFactor[rtMs1] * 1e-6))
                    if mzCorrection:
                        precursorMz = precursorMz * (1 + CorrFactorMz[precursorMz] * 1e-6)
                    precursorMzNode.attrib['value'] = str(precursorMz)

                writer.xmlWriter.write(writer.reader.element, pretty_print=True)
                writer.xmlWriter.flush()
                writer.reader.clearParsedElements()
        else:
            next(writer.reader)
    xmlSpectrumListNode.__exit__(None, None, None)
    writer.reader.clearParsedElements()
    writer.xmlWriter.write('\n')
