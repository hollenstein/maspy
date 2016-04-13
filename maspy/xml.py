from __future__ import print_function, division, unicode_literals
from future.utils import viewkeys, viewvalues, viewitems, listvalues, listitems

try: # python 2.7
    from itertools import izip as zip
except ImportError: # python 3.x series
    pass
###############################################################################
import os

import numpy
from base64 import b64decode as B64DEC
from base64 import b64encode as B64ENC
from struct import unpack as UNPACK
from struct import pack as PACK
from lxml import etree as ETREE
import io
import zlib

import maspy.auxiliary as aux
import maspy.ontology

##############################################################
### maspy unrelated xml (mzml) content #######################
##############################################################
#TODO: docstring
binaryDataArrayTypes = {'MS:1000514': 'mz', 'MS:1000515': 'i', 'MS:1000516': 'z', 'MS:1000517': 'sn',
                        'MS:1000595': 'rt', 'MS:1000617': 'lambda', 'MS:1000786': 'non-standard',
                        'MS:1000820': 'flow', 'MS:1000821': 'pressure', 'MS:1000822': 'temperature'
                        }


#TODO: docstring
mspath = aux.joinpath(os.path.dirname(aux.__file__), 'ontologies', 'psi-ms.obo')
unitpath = aux.joinpath(os.path.dirname(aux.__file__), 'ontologies', 'unit.obo')
oboTranslator = maspy.ontology.OBOOntology()
with io.open(mspath, 'r', encoding='utf-8') as openfile:
    oboTranslator.load(openfile)
with io.open(unitpath, 'r', encoding='utf-8') as openfile:
    oboTranslator.load(openfile)


###############################
# --- general xml methods --- #
###############################
def clearParsedElements(element):
    #TODO: docstring
    element.clear()
    while element.getprevious() is not None:
        del element.getparent()[0]


def clearTag(tag):
    #TODO: docstring
    return tag.split('}')[-1]


def recClearTag(element):
    #TODO: docstring
    children = element.getchildren()
    if len(children) > 0:
        for child in children:
            recClearTag(child)
    element.tag = clearTag(element.tag)


def recRemoveTreeFormating(element):
    #TODO: docstring
    children = element.getchildren()
    if len(children) > 0:
        for child in children:
            recRemoveTreeFormating(child)
    if element.text is not None:
        if len(element.text.strip()) == 0:
            element.text = None
        else:
            element.text = element.text.strip()
    if element.tail is not None:
        if len(element.tail.strip()) == 0:
            element.tail = None
        else:
            element.tail = element.tail.strip()


def recCopyElement(oldelement):
    #TODO: docstring
    #Note: doesn't copy "text" or "tail" of xml elements
    newelement = ETREE.Element(oldelement.tag, oldelement.attrib)
    if len(oldelement.getchildren()) > 0:
        for childelement in oldelement.getchildren():
            newelement.append(recCopyElement(childelement))
    return newelement


####################################
# --- working with param tuple --- #
####################################
def cvParamFromDict(attributes):
    """ Python representation of a mzML cvParam = tuple(accession, value, unitAccession) """
    keys = ['accession', 'value', 'unitAccession']
    return tuple(attributes[key] if key in attributes else None for key in keys)


def userParamFromDict(attributes):
    """ Python representation of a mzML userParam = tuple(name, value, unitAccession, type) """
    keys = ['name', 'value', 'unitAccession', 'type']
    return tuple(attributes[key] if key in attributes else None for key in keys)


def refParamGroupFromDict(attributes):
    """ Python representation of a mzML  referencableParamGroup = ('ref', ref)
    NOTE: the utilization of referencableParamGroups is currently not implemented in maspy. """
    return ('ref', attributes['ref'])


def findParam(params, targetValue):
    """ Returns a param entry (cvParam or userParam) in a list of params if its
    'accession' (cvParam) or 'name' (userParam) matches the targetValue.
    return: cvParam, userParam or None if no matching param was found """
    for param in params:
        if param[0] == targetValue:
            return param
    return None


def getParam(xmlelement):
    """ Extract parameters from an xmlelement
    param: tuple, either a 'userParam', 'cvParam' or 'referenceableParamGroupRef'
    return: param or False if xmlelement is not a parameter
    """
    elementTag = clearTag(xmlelement.tag)
    if elementTag in ['userParam', 'cvParam', 'referenceableParamGroupRef']:
        if elementTag == 'cvParam':
            param = cvParamFromDict(xmlelement.attrib)
        elif elementTag == 'userParam':
            param = userParamFromDict(xmlelement.attrib)
        else:
            param = refParamGroupFromDict(xmlelement.attrib)
    else:
        param = False
    return param


def extractParams(xmlelement):
    #TODO: docstring
    params = list()
    children = list()
    for child in xmlelement.getchildren():
        param = getParam(child)
        if param:
            params.append(param)
        else:
            children.append(child)
    return params, children


def xmlAddParams(xmlelement, params):
    #TODO: docstring
    if not params:
        return None
    for param in params:
        if len(param) == 3:
            cvAttrib = {'cvRef': param[0].split(':')[0], 'accession': param[0], 'name':oboTranslator[param[0]].name}
            if param[1]:
                cvAttrib.update({'value': param[1]})
            else:
                cvAttrib.update({'value': ''})
            if param[2]:
                cvAttrib.update({'unitAccession': param[2], 'unitCvRef': param[2].split(':')[0], 'unitName': oboTranslator[param[2]].name})
            paramElement = ETREE.Element('cvParam', **cvAttrib)
        elif len(param) == 4:
            userAttrib = {'name': param[0]}
            if param[1]:
                userAttrib.update({'value': param[1]})
            else:
                userAttrib.update({'value': ''})
            if param[2]:
                userAttrib.update({'unitAccession': param[2], 'unitCvRef': param[2].split(':')[0]})
            if param[3]:
                userAttrib.update({'type': param[3]})
            paramElement = ETREE.Element('userParam', **userAttrib)
        elif param[0] == 'ref':
            refAttrib = {'ref': param[1]}
            paramElement = ETREE.Element('referenceableParamGroupRef', **refAttrib)
        xmlelement.append(paramElement)


####################################################################
# --- decode and encode function for binary data of mzml files --- #
####################################################################
def decodeBinaryData(binaryData, arrayLength, bitEncoding, compression):
    #TODO: docstring
    bitEncodedData = binaryData.encode("utf-8")
    bitDecodedData = B64DEC(bitEncodedData)

    if bitEncoding == '64':
        floattype = 'd' # 64-bit
        numpyType = numpy.float64
    else:
        floattype = 'f' # 32-bit
        numpyType = numpy.float32

    if compression == 'zlib':
        decompressedData = zlib.decompress(bitDecodedData)
    else:
        decompressedData = bitDecodedData

    fmt = '{endian}{arraylength}{floattype}'.format(endian='<', arraylength=arrayLength, floattype=floattype)
    dataArray = numpy.array(UNPACK(fmt, decompressedData), dtype=numpyType)
    return dataArray


def encodeBinaryData(dataArray, bitEncoding, compression):
    #TODO: docstring
    arrayLength = len(dataArray)
    if bitEncoding == '64':
        floattype = 'd' # 64-bit
    else:
        floattype = 'f' # 32-bit
    fmt = '{endian}{arraylength}{floattype}'.format(endian='<', arraylength=arrayLength, floattype=floattype)
    packedData = PACK(fmt, *dataArray)

    if compression == 'zlib':
        compressedData = zlib.compress(packedData)
    else:
        compressedData = packedData

    encodedData = B64ENC(compressedData)
    return encodedData, arrayLength


def findBinaryDataType(params):
    """
    from: http://www.peptideatlas.org/tmp/mzML1.1.0.html#binaryDataArray
    a binaryDataArray "MUST supply a *child* term of MS:1000518 (binary data type) only once"
    """
    binaryDataType = None
    cvParam = None
    for param in params:
        if param[0] in binaryDataArrayTypes:
            binaryDataType = binaryDataArrayTypes[param[0]]
            cvParam = param
            break
    return binaryDataType, cvParam


def extractBinaries(binaryDataArrayList, arrayLength):
    #TODO: docstring
    extractedArrays = dict()
    arrayInfo = dict()
    for binaryData in binaryDataArrayList:
        if binaryData['binary']:
            bitEncoding = '64' if findParam(binaryData['params'], 'MS:1000523') is not None else '32'
            compression = 'zlib' if findParam(binaryData['params'], 'MS:1000574') is not None else None
            dataType, dataTypeParam = findBinaryDataType(binaryData['params'])
            extractedArrays[dataType] = decodeBinaryData(binaryData['binary'], arrayLength, bitEncoding, compression)
            binaryData['binary'] = None
            arrayInfo[dataType] = {'dataProcessingRef': None, 'params': binaryData['params']}
            if 'dataProcessingRef' in binaryData:
                arrayInfo[dataType]['dataProcessingRef'] = binaryData['dataProcessingRef']
    return extractedArrays, arrayInfo


#############################
# --- Parse a mzml file --- #
#############################
class MzmlReader(object):
    #TODO: docstring
    #TODO: change to work as a with method
    def __init__(self, mzmlPath):
        self.mzmlPath = mzmlPath
        self.metadataNode = None
        self.chromatogramList = list()
        self._parsed = False

        if self.mzmlPath.endswith('.gz'):
            #Only import modules if necessary
            import gzip; import codecs
            self.openfile = codecs.getreader('utf-8')(
                gzip.open(self.mzmlPath)
            )
        else:
            #TODO: necessary to open with 'rb'?
            self.openfile = io.open(self.mzmlPath, 'rb')
        self.iterator = ETREE.iterparse(self.openfile, events=('start', 'end'))

    def __iter__(self):
        return self

    def __next__(self):
        """ The python 2.6+ iterator """
        return self.next()

    def next(self):
        #TODO: docstring
        try:
            self.event, self.element = next(self.iterator)
            self.elementTag = clearTag(self.element.tag)
        except StopIteration:
            clearParsedElements(self.element)
            raise StopIteration
        return self.event, self.element, self.elementTag

    def loadMetadata(self):
        #TODO: docstring
        if self._parsed:
            raise TypeError('Mzml file already parsed.')
        [None for _ in self._parseMzml()]
        self._parsed = True

    def parseSpectra(self):
        #TODO: docstring
        #Note: the spectra need to be iterated completely to save the metadataNode
        if self._parsed:
            raise TypeError('Mzml file already parsed.')
        self._parsed = True
        return self._parseMzml()

    def _parseMzml(self):
        #TODO: docstring
        for event, element, elementTag in self:
            if elementTag == 'mzML':
                metadataNode = ETREE.Element(self.elementTag, self.element.attrib)
                _, _, targetTag = self.next()
                break

        while targetTag != 'mzML':
            if targetTag == 'run':
                runNode = ETREE.Element('run', self.element.attrib)
                next(self)
                while self.event != 'end' or self.elementTag != 'run':
                    if self.elementTag == 'spectrumList':
                        #Add spectrumListNode
                        specListAttrib = {'defaultDataProcessingRef': self.element.attrib['defaultDataProcessingRef']}
                        specListNode = ETREE.Element('spectrumList', specListAttrib)
                        runNode.append(specListNode)
                        #Parse and yield spectrum xml elements
                        while self.event != 'end' or self.elementTag != 'spectrumList':
                            if self.event == 'end' and self.elementTag == 'spectrum':
                                yield self.element
                                clearParsedElements(self.element)
                            next(self)
                    elif self.elementTag == 'chromatogramList':
                        #Add chromatogramListNode
                        chromListAttrib = {'defaultDataProcessingRef': self.element.attrib['defaultDataProcessingRef']}
                        chromListNode = ETREE.Element('chromatogramList', chromListAttrib)
                        runNode.append(chromListNode)
                        #Parse and store chromatogram xml elements
                        while self.event != 'end' or self.elementTag != 'chromatogramList':
                            if self.event == 'end' and self.elementTag == 'chromatogram':
                                self.chromatogramList.append(self.element)
                                #Alternatively also the chromatogram xml elements could be yielded:
                                #yield self.element
                                #clearParsedElements(self.element)
                            next(self)
                    else:
                        runNode.append(self.element)
                    next(self)
                metadataNode.append(runNode)
                break
            else:
                while self.event != 'end' or self.elementTag != targetTag:
                    next(self)
                metadataNode.append(self.element)
            _, _, targetTag = next(self)
        recClearTag(metadataNode)
        recRemoveTreeFormating(metadataNode)
        self.metadataNode = recCopyElement(metadataNode)
        self.openfile.close()


def sublistReader(xmlelement):
    #TODO: docstring
    #Note: actually I'm not 100% sure how this function behaves
    elements = list()
    params, children = extractParams(xmlelement)
    for child in children:
        currElement = dict()
        currElement.update(child.attrib)
        childparams, subchildren = extractParams(child)
        if childparams:
            currElement['params'] = childparams
        for subchild in subchildren:
            subchildTag = clearTag(subchild.tag)
            if 'List' in subchildTag:
                listelements, listparams = sublistReader(subchild)
                simplelist = [listelement['params'] for listelement in listelements]
                currElement[subchildTag] = simplelist
            else:
                subchildparams, _ = extractParams(subchild)
                currElement[subchildTag] = subchildparams
                if subchildTag == 'binary' and subchild.text:
                    currElement[subchildTag] = subchild.text.strip()
        elements.append(currElement)
    return elements, params
