"""
#TODO: module description
"""
######################## Python 2 and 3 compatibility #########################
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
from future.utils import viewitems, viewkeys, viewvalues, listitems, listvalues

try:
    #python 2.7
    from itertools import izip as zip
except ImportError:
    #python 3 series
    pass
################################################################################
from base64 import b64decode as B64DEC
from base64 import b64encode as B64ENC
import io
from struct import unpack as UNPACK
from struct import pack as PACK
import os
import zlib

import numpy
from lxml import etree as ETREE

import maspy.auxiliary as aux
import maspy.ontology as ONTOLOGY

##############################################################
### maspy unrelated xml (mzml) content #######################
##############################################################
""" A dictionary with the psi-ms.obo ids of various binary data array types
(derived from id: MS:1000513) used to describe a data array of either a mzML
spectrum element or mzML chromatogram element.
"""
#TODO: check exact term of spectrum or chromatogram in mzml context
binaryDataArrayTypes = {'MS:1000514': 'mz', 'MS:1000515': 'i',
                        'MS:1000516': 'z', 'MS:1000517': 'sn',
                        'MS:1000595': 'rt', 'MS:1000617': 'lambda',
                        'MS:1000786': 'non-standard', 'MS:1000820': 'flow',
                        'MS:1000821': 'pressure', 'MS:1000822': 'temperature'
                        }


"""#TODO: docstring """
oboTranslator = ONTOLOGY.OBOOntology()
#Note: duplicate values have to be removed from the ".obo" files, to prevent
#errors in ONTOLOGY.OBOOntology()
msOntologyPath = aux.joinpath(os.path.dirname(aux.__file__), 'ontologies',
                              'psi-ms.obo'
                              )
unitOntologyPath = aux.joinpath(os.path.dirname(aux.__file__), 'ontologies',
                                'unit.obo'
                                )
with io.open(msOntologyPath, 'r', encoding='utf-8') as openfile:
    oboTranslator.load(openfile)
with io.open(unitOntologyPath, 'r', encoding='utf-8') as openfile:
    oboTranslator.load(openfile)


###############################
# --- general xml methods --- #
###############################
def clearParsedElements(element):
    """Deletes an element and all linked parent elements.

    This function is used to save memory while iteratively parsing
    an xml file by removing already processed elements.

    :param element: #TODO docstring
    """
    element.clear()
    while element.getprevious() is not None:
        del element.getparent()[0]


def clearTag(tag):
    """ #TODO: docstring
    eg "{http://psi.hupo.org/ms/mzml}mzML" returns "mzML"

    :param tag: #TODO docstring
    :returns:
    """
    return tag.split('}')[-1]


def recClearTag(element):
    """Applies maspy.xml.clearTag() to the tag attribute of the "element" and
    recursively to all child elements.

    :param element: an :instance:`xml.etree.Element`
    """
    children = element.getchildren()
    if len(children) > 0:
        for child in children:
            recClearTag(child)
    element.tag = clearTag(element.tag)


def recRemoveTreeFormating(element):
    """Removes whitespace characters, which are leftovers from previous xml
    formatting.

    :param element: an instance of lxml.etree._Element

    str.strip() is applied to the "text" and the "tail" attribute of the
    element and recursively to all child elements.
    """
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
    """Generates a copy of an xml element and recursively of all
    child elements.

    :param oldelement: an instance of lxml.etree._Element

    :returns: a copy of the "oldelement"

    .. warning::
        doesn't copy ``.text`` or ``.tail`` of xml elements
    """
    newelement = ETREE.Element(oldelement.tag, oldelement.attrib)
    if len(oldelement.getchildren()) > 0:
        for childelement in oldelement.getchildren():
            newelement.append(recCopyElement(childelement))
    return newelement


####################################
# --- working with param tuple --- #
####################################
def cvParamFromDict(attributes):
    """Python representation of a mzML cvParam = tuple(accession, value,
    unitAccession).

    :param attributes: #TODO: docstring

    :returns: #TODO: docstring
    """
    keys = ['accession', 'value', 'unitAccession']
    return tuple(attributes[key] if key in attributes else None for key in keys)


def userParamFromDict(attributes):
    """Python representation of a mzML userParam = tuple(name, value,
    unitAccession, type)

    :param attributes: #TODO: docstring

    :returns: #TODO: docstring
    """
    keys = ['name', 'value', 'unitAccession', 'type']
    return tuple(attributes[key] if key in attributes else None for key in keys)


def refParamGroupFromDict(attributes):
    """Python representation of a mzML  referencableParamGroup = ('ref', ref)

    :param attributes: #TODO: docstring

    :returns: #TODO: docstring

    .. note::
        altough the mzML element referencableParamGroups is imported, its
        utilization is currently not implemented in MasPy.
    """
    return ('ref', attributes['ref'])


def findParam(params, targetValue):
    """Returns a param entry (cvParam or userParam) in a list of params if its
    'accession' (cvParam) or 'name' (userParam) matches the targetValue.
    return: cvParam, userParam or None if no matching param was found

    :param params: #TODO: docstring
    :param targetValue: #TODO: docstring

    :returns: #TODO: docstring
    """
    for param in params:
        if param[0] == targetValue:
            return param
    return None


def getParam(xmlelement):
    """Converts an mzML xml element to a param tuple.

    :param xmlelement: #TODO docstring

    :returns: a param tuple or False if the xmlelement is not a parameter
        ('userParam', 'cvParam' or 'referenceableParamGroupRef')
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
    """ #TODO docstring

    :param xmlelement: #TODO docstring

    :returns: #TODO docstring
    """
    params = list()
    children = list()
    for child in xmlelement.getchildren():
        param = getParam(child)
        if param:
            params.append(param)
        else:
            children.append(child)
    return params, children


def xmlAddParams(parentelement, params):
    """Generates new mzML parameter xml elements and adds them to the
    'parentelement' as xml children elements.

    :param parentelement: :class:`xml.etree.Element`, an mzML element
    :param params: a list of mzML parameter tuples ('cvParam', 'userParam' or
        'referencableParamGroup')
    """
    if not params:
        return None
    for param in params:
        if len(param) == 3:
            cvAttrib = {'cvRef': param[0].split(':')[0], 'accession': param[0],
                        'name':oboTranslator[param[0]].name
                        }
            if param[1]:
                cvAttrib.update({'value': param[1]})
            else:
                cvAttrib.update({'value': ''})
            if param[2]:
                cvAttrib.update({'unitAccession': param[2],
                                 'unitCvRef': param[2].split(':')[0],
                                 'unitName': oboTranslator[param[2]].name
                                 })
            paramElement = ETREE.Element('cvParam', **cvAttrib)
        elif len(param) == 4:
            userAttrib = {'name': param[0]}
            if param[1]:
                userAttrib.update({'value': param[1]})
            else:
                userAttrib.update({'value': ''})
            if param[2]:
                userAttrib.update({'unitAccession': param[2],
                                   'unitCvRef': param[2].split(':')[0]
                                   })
            if param[3]:
                userAttrib.update({'type': param[3]})
            paramElement = ETREE.Element('userParam', **userAttrib)
        elif param[0] == 'ref':
            refAttrib = {'ref': param[1]}
            paramElement = ETREE.Element('referenceableParamGroupRef',
                                         **refAttrib
                                         )
        parentelement.append(paramElement)


####################################################################
# --- decode and encode function for binary data of mzml files --- #
####################################################################
def interpretBitEncoding(bitEncoding):
    """Returns a floattype string and a numpy array type.

    :param bitEncoding: Must be either '64' or '32'

    :returns: (floattype, numpyType)
    """
    if bitEncoding == '64':
        floattype = 'd' # 64-bit
        numpyType = numpy.float64
    elif bitEncoding == '32':
        floattype = 'f' # 32-bit
        numpyType = numpy.float32
    else:
        errorText = ''.join(['bitEncoding \'', bitEncoding, '\' not defined. ',
                             'Must be \'64\' or \'32\''
                             ])
        raise TypeError(errorText)
    return (floattype, numpyType)


def decodeBinaryData(binaryData, arrayLength, bitEncoding, compression):
    """Function to decode a mzML byte array into a numpy array. This is the
    inverse function of :func:`encodeBinaryData`. Concept inherited from
    :func:`pymzml.spec.Spectrum._decode` of the python library `pymzML
    <https://pymzml.github.io/>`_.

    :param binaryData: #TODO: docstring
    :param arrayLength: #TODO: docstring
    :param binEncoding: #TODO: docstring
    :param compression: #TODO: docstring

    :returns: #TODO: docstring
    """
    #TODO: should raise an error if a wrong compression is specified
    bitEncodedData = binaryData.encode("utf-8")
    bitDecodedData = B64DEC(bitEncodedData)
    floattype, numpyType = interpretBitEncoding(bitEncoding)

    if compression == 'zlib':
        decompressedData = zlib.decompress(bitDecodedData)
    else:
        decompressedData = bitDecodedData

    fmt = '{endian}{arraylength}{floattype}'.format(endian='<',
                                                    arraylength=arrayLength,
                                                    floattype=floattype
                                                    )
    dataArray = numpy.array(UNPACK(fmt, decompressedData), dtype=numpyType)
    return dataArray


def encodeBinaryData(dataArray, bitEncoding, compression):
    """Function to encode a ``numpy.array`` into a mzML byte array. This is the
    inverse function of :func:`decodeBinaryData`.

    :param dataArray: #TODO: docstring
    :param bitEncoding: #TODO: docstring
    :param compression: #TODO: docstring

    :returns: #TODO: docstring
    """
    #TODO: should raise an error if a wrong compression is specified
    arrayLength = len(dataArray)
    floattype, __ = interpretBitEncoding(bitEncoding)
    fmt = '{endian}{arraylength}{floattype}'.format(endian='<',
                                                    arraylength=arrayLength,
                                                    floattype=floattype
                                                    )
    packedData = PACK(fmt, *dataArray)

    if compression == 'zlib':
        compressedData = zlib.compress(packedData)
    else:
        compressedData = packedData

    encodedData = B64ENC(compressedData)
    return encodedData, arrayLength


def findBinaryDataType(params):
    """ #TODO: docstring
    from: http://www.peptideatlas.org/tmp/mzML1.1.0.html#binaryDataArray
    a binaryDataArray "MUST supply a *child* term of MS:1000518
    (binary data type) only once"

    :param params: #TODO: docstring

    :returns: #TODO: docstring
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
    """ #TODO: docstring

    :param binaryDataArrayList: #TODO: docstring
    :param arrayLength: #TODO: docstring

    :returns: #TODO: docstring
    """
    extractedArrays = dict()
    arrayInfo = dict()
    for binaryData in binaryDataArrayList:
        if findParam(binaryData['params'], 'MS:1000523') is not None:
            bitEncoding = '64'
        else:
            bitEncoding = '32'
        if findParam(binaryData['params'], 'MS:1000574') is not None:
            compression = 'zlib'
        else:
            compression = None
        dataType, dataTypeParam = findBinaryDataType(binaryData['params'])
        if binaryData['binary']:
            extractedArrays[dataType] = decodeBinaryData(binaryData['binary'],
                                                         arrayLength,
                                                         bitEncoding,
                                                         compression
                                                         )
        else:
            __, numpyType = interpretBitEncoding(bitEncoding)
            extractedArrays[dataType] = numpy.array([], dtype=numpyType)

        binaryData['binary'] = None
        arrayInfo[dataType] = {'dataProcessingRef': None,
                               'params': binaryData['params']
                               }
        if 'dataProcessingRef' in binaryData:
            arrayInfo[dataType]['dataProcessingRef'] = \
                binaryData['dataProcessingRef']
    return extractedArrays, arrayInfo


#############################
# --- Parse a mzml file --- #
#############################
class MzmlReader(object):
    """ #TODO: docstring

    :ivar mzmlPath: #TODO: docstring
    :ivar metadataNode: #TODO: docstring
    :ivar chromatogramList: #TODO: docstring

    """
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
        """ #TODO: docstring

        :returns: #TODO: docstring
        """
        try:
            self.event, self.element = next(self.iterator)
            self.elementTag = clearTag(self.element.tag)
        except StopIteration:
            clearParsedElements(self.element)
            raise StopIteration
        return self.event, self.element, self.elementTag

    def loadMetadata(self):
        """ #TODO: docstring """
        #TODO: change that spectra dont have to be iterated to extract metadata
        #node
        if self._parsed:
            raise TypeError('Mzml file already parsed.')
        [None for _ in self._parseMzml()]
        self._parsed = True

    def parseSpectra(self):
        """ #TODO: docstring

        :returns: #TODO: docstring
        """
        #Note: the spectra need to be iterated completely to save the
        #metadataNode
        if self._parsed:
            raise TypeError('Mzml file already parsed.')
        self._parsed = True
        return self._parseMzml()

    def _parseMzml(self):
        """ #TODO: docstring """
        #TODO: this is already pretty nested, reduce that eg by using a function
        #   processRunNode
        for event, element, elementTag in self:
            if elementTag == 'mzML':
                metadataNode = ETREE.Element(self.elementTag,
                                             self.element.attrib
                                             )
                _, _, targetTag = self.next()
                break

        while targetTag != 'mzML':
            if targetTag == 'run':
                runNode = ETREE.Element('run', self.element.attrib)
                next(self)
                while self.event != 'end' or self.elementTag != 'run':
                    if self.elementTag == 'spectrumList':
                        #Add spectrumListNode
                        specListAttrib = {'defaultDataProcessingRef':
                                          self.element.attrib['defaultDataProcessingRef']
                                          }
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
                        chromListAttrib = {'defaultDataProcessingRef':
                                           self.element.attrib['defaultDataProcessingRef']
                                           }
                        chromListNode = ETREE.Element('chromatogramList',
                                                      chromListAttrib
                                                      )
                        runNode.append(chromListNode)
                        #Parse and store chromatogram xml elements
                        while self.event != 'end' or self.elementTag != 'chromatogramList':
                            if self.event == 'end' and self.elementTag == 'chromatogram':
                                self.chromatogramList.append(self.element)
                                #Alternatively also the chromatogram xml
                                #elements could be yielded:
                                #   yield self.element
                                #   clearParsedElements(self.element)
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
    """ #TODO: docstring """
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
