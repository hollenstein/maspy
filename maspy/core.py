""" Modul docstring """

from __future__ import print_function, division, unicode_literals
from future.utils import viewkeys, viewvalues, viewitems, listvalues, listitems

try: # python 2.7
    from itertools import izip as zip
except ImportError: # python 3.x series
    pass
###############################################################################
from collections import defaultdict as ddict
import io
import numpy
from operator import itemgetter as ITEMGETTER
import os
import warnings

import json
from lxml import etree as ETREE
import zipfile

import maspy.auxiliary as aux
import maspy.constants
import maspy.peptidemethods


############################################
### Common container functions #############
############################################
def _getArrays(items, attr, defaultValue):
    """Return arrays with equal size of item attributes from a list of sorted "items"
    for fast and convenient data processing.

    :param attr: list of item attributes that should be added to the returned array.
    :param defaultValue: if an item is missing an attribute, the "defaultValue" is
        added to the array instead.

    return {'attribute1': numpy.array(), 'attribute2': numpy.array(), ...}
    """
    arrays = dict([(key, []) for key in attr])
    for item in items:
        for key in attr:
            arrays[key].append(getattr(item, key, defaultValue))
    for key in [_ for _ in viewkeys(arrays)]:
        arrays[key] = numpy.array(arrays[key])

    return arrays


def _getItems(container, containerKeys=None, sort=False, reverse=False, selector=lambda item: True):
    """Generator that yields filtered and/or sorted items from the specified "container".

    :param containerKeys: valid keys of the "container", if None all keys are considered
    :type containerKeys: a single dictionary key or a list of keys
    :param sort: if "sort" is specified the returned list of items is sorted according to the item
        attribute specified by "sort", if the attribute is not present the item is skipped.
    :param reverse: boolean to reverse sort order
    :param selector: a function which is called with each item and returns
        True (include item) or False (discard item). If not specified all items are returned
    """
    containerKeys = [_ for _ in viewkeys(container)] if containerKeys is None else aux.toList(containerKeys)

    if sort:
        sortIdentifier = list()
        for containerKey in containerKeys:
            for identifier in [_ for _ in viewkeys(container[containerKey])]:
                item = container[containerKey][identifier]
                if selector(item):
                    try:
                        sortIdentifier.append((getattr(item, sort), containerKey, identifier))
                    except AttributeError:
                        pass
        sortIdentifier.sort(key=ITEMGETTER(0), reverse=reverse)
        for _, containerKey, identifier in sortIdentifier:
            yield container[containerKey][identifier]
    else:
        for containerKey in containerKeys:
            for identifier in [_ for _ in viewkeys(container[containerKey])]:
                item = container[containerKey][identifier]
                if selector(item):
                    yield item


def _getListItems(container, containerKeys=None, sort=False, reverse=False, selector=lambda item: True):
    """Generator that yields filtered and/or sorted items from the specified "container",
    Note: use this function if the value of the container is not the item itself but a list of items

    :param containerKeys: valid keys of the "container", if None all keys are considered
    :type containerKeys: a single dictionary key or a list of keys
    :param sort: if "sort" is specified the returned list of items is sorted according to the item
    attribute specified by "sort", if the attribute is not present the item is skipped.
    :param reverse: boolean to reverse sort order
    :param selector: a function which is called with each item and returns
    True (include item) or False (discard item). If not specified all items are returned
    """
    containerKeys = [_ for _ in viewkeys(container)] if containerKeys is None else aux.toList(containerKeys)
    if sort:
        sortIdentifier = list()
        for containerKey in containerKeys:
            for identifier in [_ for _ in viewkeys(container[containerKey])]:
                for itemPos, item in enumerate(container[containerKey][identifier]):
                    if selector(item):
                        try:
                            sortIdentifier.append((getattr(item, sort), containerKey, identifier, itemPos))
                        except AttributeError:
                            pass
        sortIdentifier.sort(key=ITEMGETTER(0), reverse=reverse)
        for _, containerKey, identifier, itemPos in sortIdentifier:
            yield container[containerKey][identifier][itemPos]
    else:
        for containerKey in containerKeys:
            for identifier in [_ for _ in viewkeys(container[containerKey])]:
                for itemPos, item in enumerate(container[containerKey][identifier]):
                    if selector(item):
                        yield item


def _containerSetPath(container, folderpath, specfiles):
    """Changse the folderpath of the specified specfiles in container.info """
    if not os.path.exists(folderpath):
        warnings.warn('The specified directory does not exist %s' %(folderpath, ))
    for specfile in aux.toList(specfiles):
        if specfile in container.info:
            container.info[specfile]['path'] = folderpath
        else:
            warnings.warn('Specfile not present in container %s' %(specfile, ))


##############################################################
### MsrunContainer related classes and functions #############
##############################################################
class MsrunContainer(object):
    """Container for mass spectrometry data (eg MS1 and MS2 spectra), provides full support for mzml files.

    :ivar rmc: "run metadata container", contains mzml metadata xml strings, as imported from the mzML file
    :ivar cic: "chromatogram item container", see :class:`Ci`
    :ivar smic: "spectrum metadata item container", see :class:`Smi`
    :ivar saic: "spectrum array item container", see :class:`Sai`
    :ivar sic: "spectrum item container", see :class:`Si`
    :ivar info: contains information about the imported specfiles. ::

            {specfilename: {'path': str,
                            'status': {u'ci': bool, u'rm': bool, u'sai': bool,
                                       u'si': bool, u'smi': bool}
                            },
             ...
             }

        ``path`` contains information about the filelocation used for saving and
        loading msrun files in the maspy dataformat. ``status`` describes which
        datatypes are curerently imported.

        code example::

            {u'JD_06232014_sample1_A': {u'path': u'C:/filedirectory',
                                        u'status': {u'ci': True,
                                                    u'rm': True,
                                                    u'sai': True,
                                                    u'si': True,
                                                    u'smi': True
                                                    }
                                        }
             }

    .. note::
        The structure for the containers ``rmc``, ``cic``, ``smic``, ``saic``
        and ``sic`` is::

            {"specfilename": {"itemId": item, ...}, ...}


    """
    def __init__(self):
        self.rmc = {}
        self.cic = {}
        self.smic = {}
        self.saic = {}
        self.sic = {}
        self.info = {}

    def getArrays(self, attr=None, specfiles=None, sort=False, reverse=False, selector=None, defaultValue=None):
        """Return a condensed array of data selected from :class:`Si` objects from ``self.sic``
        for fast and convenient data processing.

        :param attr: list of :class:`Si` item attributes that should be added to the returned array.
            The attributes "id" and "specfile" are always included, in combination they serve as a unique id.
        :param defaultValue: if an item is missing an attribute, the "defaultValue" is
            added to the array instead.
        :param specfiles: filenames of msrun files - if specified return only items from those files
        :type specfiles: str or [str, str, ...]
        :param sort: if "sort" is specified the returned list of items is sorted according to the :class:`Si`
            attribute specified by "sort", if the attribute is not present the item is skipped.
        :param reverse: boolean to reverse sort order
        :param selector: a function which is called with each :class:`Si` item and has to return
            True (include item) or False (discard item). Default function is: ``lambda si: True``

        :returns: {'attribute1': numpy.array(), 'attribute1': numpy.array(), ...}
        """
        selector = lambda si: True if selector is None else selector
        attr = attr if attr is not None else []
        attr = set(['id', 'specfile'] + aux.toList(attr))
        items = self.getItems(specfiles, sort, reverse, selector)
        return _getArrays(items, attr, defaultValue)

    def getItems(self, specfiles=None, sort=False, reverse=False, selector=None):
        """Generator that yields filtered and/or sorted :class:`Si` instances
        from :class:`self.sic <maspy.core.MsrunContainer>`.

        :param specfiles: filenames of msrun files - if specified return only items from those files
        :type specfiles: str or [str, str, ...]
        :param sort: if "sort" is specified the returned list of items is sorted according to the :class:`Si`
            attribute specified by "sort", if the attribute is not present the item is skipped.
        :param reverse: boolean to reverse sort order
        :param selector: a function which is called with each :class:`Si` item and returns
            True (include item) or False (discard item). Default function is: ``lambda si: True``
        """
        selector = lambda si: True if selector is None else selector
        specfiles = [_ for _ in viewkeys(self.info)] if specfiles is None else aux.toList(specfiles)
        return _getItems(self.sic, specfiles, sort, reverse, selector)

    def getItem(self, specfile, identifier):
        """Returns a :class:`Si` instance from
        :class:`self.sic <maspy.core.MsrunContainer>`."""
        return self.sic[specfile][identifier]

    def addSpecfile(self, specfiles, path):
        """Adds specfile entries to self.info, but doesn't import any data yet.
        To actually import the MsrunContainer files, use :func:`MsrunContainer.load()`.
        """
        for specfile in aux.toList(specfiles):
            if specfile not in self.info:
                self._addSpecfile(specfile, path)
            else:
                warnings.warn('Specfile is already present in the MsrunContainer.\nname: %s \npath: %s' %(specfile, self.info[specfile]['path']))

    def _addSpecfile(self, specfile, path):
        """Adds a new specfile entry to MsrunContainer.info. """
        datatypeStatus = {'rm':False, 'ci':False, 'smi':False, 'sai':False, 'si':False}
        self.info[specfile] = {'path': path, 'status': datatypeStatus}

    def setPath(self, folderpath, specfiles=None):
        """Change the folderpath of the specified specfiles. If save is called and no container file is
        present in the specified path, a new file is generated, otherwise it is replaced.
        """
        specfiles = [_ for _ in viewkeys(self.info)] if specfiles is None else specfiles
        _containerSetPath(self, folderpath, specfiles)

    def removeData(self, specfiles=None, rm=False, ci=False, smi=False, sai=False, si=False):
        """Removes the specified datatypes of the specfiles from the msrunContainer.
        To completely remove the specfile, also from info, use :func:`MsrunContainer.removeSpecfile`
        """
        specfiles = [_ for _ in viewkeys(self.info)] if specfiles is None else specfiles
        #TODO: add check if specfiles are present in the container
        datatypes = self._processDatatypes(rm, ci, smi, sai, si)
        for specfile in aux.toList(specfiles):
            for datatype in datatypes:
                datatypeContainer = datatype+'c'
                dataContainer = getattr(self, datatypeContainer)
                try:
                    del dataContainer[specfile]
                except KeyError:
                    pass
                finally:
                    self.info[specfile]['status'][datatype] = False

    def removeSpecfile(self, specfiles):
        """Completely removes the specified specfiles from the msrunContainer."""
        for specfile in aux.toList(specfiles):
            for datatypeContainer in ['rmc', 'cic', 'smic', 'saic', 'sic']:
                dataContainer = getattr(self, datatypeContainer)
                try:
                    del dataContainer[specfile]
                except KeyError:
                    pass
            del self.info[specfile]

    def _processDatatypes(self, rm, ci, smi, sai, si):
        #TODO: docstring
        datatypes = list()
        for datatype, value in [('rm', rm), ('ci', ci), ('smi', smi), ('sai', sai), ('si', si)]:
            if value:
                datatypes.append(datatype)
        return datatypes

    def save(self, specfiles=None, rm=False, ci=False, smi=False, sai=False, si=False, compress=True, path=None):
        """

        :param specfiles:
        :param rm:
        :param ci:
        :param smi:
        :param sai:
        :param si:
        :param compress:
        :param path:
        """
        #TODO: docstring
        specfiles = [_ for _ in viewkeys(self.info)] if specfiles is None else specfiles
        datatypes = self._processDatatypes(rm, ci, smi, sai, si)
        if len(datatypes) == 0:
            datatypes = ['rm', 'ci', 'smi', 'sai', 'si']

        for specfile in aux.toList(specfiles):
            if specfile not in self.info:
                print('Error while saving', specfile, ', not found in msrunCountainer!')
                continue
            else:
                msrunInfo = self.info[specfile]
                specfilePath = msrunInfo['path'] if path is None else path

            with aux.PartiallySafeReplace() as msr:
                for datatype in datatypes:
                    filename = specfile + '.mrc_' + datatype
                    filepath = aux.joinpath(specfilePath, filename)
                    with msr.open(filepath, 'w+b') as openfile:
                        if datatype == 'rm':
                           self._writeRmc(openfile, specfile, compress=compress)
                        elif datatype == 'ci':
                           self._writeCic(openfile, specfile, compress=compress)
                        elif datatype == 'si':
                           self._writeSic(openfile, specfile, compress=compress)
                        elif datatype == 'smi':
                           self._writeSmic(openfile, specfile, compress=compress)
                        elif datatype == 'sai':
                           self._writeSaic(openfile, specfile, compress=compress)

    def _writeCic(self, filelike, specfile, compress=True):
        #TODO: docstring
        aux.writeBinaryItemContainer(filelike, self.cic[specfile], compress=compress)

    def _writeSaic(self, filelike, specfile, compress=True):
        #TODO: docstring
        aux.writeBinaryItemContainer(filelike, self.saic[specfile], compress=compress)

    def _writeSmic(self, filelike, specfile, compress=True):
        #TODO: docstring
        aux.writeJsonZipfile(filelike, self.smic[specfile], compress=compress)

    def _writeSic(self, filelike, specfile, compress=True):
        #TODO: docstring
        aux.writeJsonZipfile(filelike, self.sic[specfile], compress=compress)

    def _writeRmc(self, filelike, specfile, compress=True):
        #TODO: docstring
        xmlString = ETREE.tostring(self.rmc[specfile], pretty_print=True)
        filelike.write(xmlString)

    def load(self, specfiles=None, rm=False, ci=False, smi=False, sai=False, si=False):
        """Import the specified datatypes from specfiles"""
        #TODO: docstring
        specfiles = [_ for _ in viewkeys(self.info)] if specfiles is None else specfiles
        datatypes = self._processDatatypes(rm, ci, smi, sai, si)
        if len(datatypes) == 0:
            datatypes = ['rm', 'ci', 'smi', 'sai', 'si']

        for specfile in aux.toList(specfiles):
            if specfile not in self.info:
                print(''.join(('Error while loading "', specfile, '", not present in msrunCountainer.info')))
                continue
            else:
                msrunInfo = self.info[specfile]
                specfilePath = msrunInfo['path']

            if 'rm' in datatypes:
                rmPath = aux.joinpath(specfilePath, specfile+'.mrc_rm')
                with open(rmPath, 'rb') as openfile:
                    xmlString = openfile.read()
                self.rmc[specfile] = ETREE.fromstring(xmlString)
                msrunInfo['status']['rm'] = True

            if 'ci' in datatypes:
                ciPath = aux.joinpath(specfilePath, specfile+'.mrc_ci')
                self.cic[specfile] = aux.loadBinaryItemContainer(ciPath, Ci.jsonHook)
                msrunInfo['status']['ci'] = True

            if 'smi' in datatypes:
                smiPath = aux.joinpath(specfilePath, specfile+'.mrc_smi')
                with zipfile.ZipFile(smiPath, 'r') as containerZip:
                    #Convert the zipfile data into a str object, necessary since containerZip.read() returns a bytes object.
                    jsonString = io.TextIOWrapper(containerZip.open('data'), encoding='utf-8').read()
                self.smic[specfile] = json.loads(jsonString, object_hook=Smi.jsonHook)
                msrunInfo['status']['smi'] = True

            if 'sai' in datatypes:
                saiPath = aux.joinpath(specfilePath, specfile+'.mrc_sai')
                self.saic[specfile] = aux.loadBinaryItemContainer(saiPath, Sai.jsonHook)
                msrunInfo['status']['sai'] = True

            if 'si' in datatypes:
                siPath = aux.joinpath(specfilePath, specfile+'.mrc_si')
                with zipfile.ZipFile(siPath, 'r') as containerZip:
                    #Convert the zipfile data into a str object, necessary since containerZip.read() returns a bytes object.
                    jsonString = io.TextIOWrapper(containerZip.open('data'), encoding='utf-8').read()
                self.sic[specfile] = json.loads(jsonString, object_hook=Si.jsonHook)
                msrunInfo['status']['si'] = True


class Ci(object):
    """Chromatogram item (Ci), representation of a mzML chromatogram.

    :ivar id: The unique id of this chromatogram. Typically descriptive
        for the chromatogram, eg. "TIC".
    :ivar dataProcessingRef: This attribute can optionally reference the 'id' of
        the appropriate dataProcessing, from mzML.
    :ivar precursor: The method of precursor ion selection and activation, from
        mzML.
    :ivar product: The method of product ion selection and activation in a
        precursor ion scan, from mzML.
    :ivar params: A list of parameter tuple, TODO: as described elsewhere
    :ivar arrays: dictionary of nummpy arrays containing the chromatogram data
        points. Keys are derived from the specified cvParam, see
        :func:`maspy.xml.findBinaryDataType`.
    :ivar arrayInfo: dictionary describing each dataType present in ``.arrays``. ::

            {dataType: {'dataProcessingRef': str,
                        'params': [paramTuple, paramTuple, ...]
                        }
             }

        code example::

            {u'i': {u'dataProcessingRef': None,
                    u'params': [('MS:1000521', '', None),
                                ('MS:1000574', '', None),
                                ('MS:1000515', '', 'MS:1000131')
                                ]
                    },
             u'rt': {u'dataProcessingRef': None,
                     u'params': [('MS:1000523', '', None),
                                 ('MS:1000574', '', None),
                                 ('MS:1000595', '', 'UO:0000031')
                                 ]
                     }
             }
    """
    __slots__ = ['id', 'dataProcessingRef', 'precursor', 'product', 'params', 'attrib', 'arrays', 'arrayInfo']

    def __init__(self):
        self.id = str()
        self.dataProcessingRef = None
        self.precursor = None
        self.product = None
        self.params = list()
        self.attrib = dict()
        self.arrays = dict()
        self.arrayInfo = dict()

    def _reprJSON(self):
        return {'__Ci__': (self.id, self.dataProcessingRef, self.precursor, self.product,
                           self.params, self.attrib, self.arrayInfo
                           )}

    @classmethod
    def _fromJSON(cls, jsonobject):
        newInstance = cls()
        attribDict = {}
        attribDict['id'] = jsonobject[0]
        attribDict['dataProcessingRef'] = jsonobject[1]
        attribDict['precursor'] = jsonobject[2]
        attribDict['product'] = jsonobject[3]
        attribDict['params'] = [tuple(param) for param in jsonobject[4]]
        attribDict['attrib'] = jsonobject[5]
        attribDict['arrayInfo'] = dict()
        for arrayType in jsonobject[6]:
            attribDict['arrayInfo'][arrayType] = {'dataProcessingRef': jsonobject[6][arrayType]['dataProcessingRef'],
                                                  'params': [tuple(param) for param in jsonobject[6][arrayType]['params']]
                                                  }
        for key, value in viewitems(attribDict):
            setattr(newInstance, key, value)
        return newInstance

    @staticmethod
    def jsonHook(encoded):
        if '__Ci__' in encoded:
            return Ci._fromJSON(encoded['__Ci__'])
        elif '__MzmlProduct__' in encoded:
            return MzmlProduct._fromJSON(encoded['__MzmlProduct__'])
        elif '__MzmlPrecursor__' in encoded:
            return MzmlPrecursor._fromJSON(encoded['__MzmlPrecursor__'])
        else:
            return encoded


class Sai(object):
    """Spectrum array item (Sai)
    Includes all spectrum information provided by a mzML file, excluding the actual data arrays.

    :ivar id: The unique id of this spectrum, typically the scan number. Is used
        together with "specfile" as a key to access the spectrum in its container
        :class:`maspy.core.SaiContainer`. Should be derived from the spectrums
        nativeID format (MS:1000767)

    specfile: An id representing a group of spectra, typically of the same mzML file.
        Is used together with "identifier" as a key to access the spectrum in its container :class:`maspy.core.SiContainer`
    """
    __slots__ = ['id', 'specfile', 'arrays', 'arrayInfo']

    def __init__(self, identifier, specfile):
        self.id = identifier
        self.specfile = specfile

        self.arrays = dict()
        self.arrayInfo = dict()

    def _reprJSON(self):
        return {'__Sai__': (self.id, self.specfile, self.arrayInfo)}

    @classmethod
    def _fromJSON(cls, jsonobject):
        newInstance = cls(jsonobject[0], jsonobject[1])
        for arrayType in jsonobject[2]:
            newInstance.arrayInfo[arrayType] = {'dataProcessingRef': jsonobject[2][arrayType]['dataProcessingRef'],
                                                'params': [tuple(param) for param in jsonobject[2][arrayType]['params']]
                                                 }
        return newInstance

    @staticmethod
    def jsonHook(encoded):
        if '__Sai__' in encoded:
            return Sai._fromJSON(encoded['__Sai__'])
        else:
            return encoded


class Smi(object):
    """Spectrum metadata item (Smi)
    Includes all spectrum information provided by a mzML file, excluding the actual data arrays.

    identifier: The unique id of this spectrum, typically the scan number.
        Is used together with "specfile" as a key to access the spectrum in its container :class:`maspy.core.SiContainer`
        Should be derived from the spectrums nativeID format (MS:1000767)

    specfile: An id representing a group of spectra, typically of the same mzML file.
        Is used together with "identifier" as a key to access the spectrum in its container :class:`maspy.core.SiContainer`
    """
    __slots__ = ['id', 'specfile', 'attributes', 'params', 'scanListParams',
                 'scanList', 'precursorList', 'productList'
                 ]

    def __init__(self, identifier, specfile):
        #super(SpectrumItem, self).__init__(identifier, specfile)
        self.id = identifier
        self.specfile = specfile
        self.attributes = dict()
        self.params = dict()
        self.scanListParams = list()
        self.scanList = list()
        self.precursorList = list()
        self.productList = list()

    def _reprJSON(self):
        return {'__Smi__': (self.id, self.specfile, self.attributes, self.params, self.scanListParams,
                            self.scanList, self.precursorList, self.productList
                            )}

    @classmethod
    def _fromJSON(cls, jsonobject):
        newInstance = cls(None, None)
        attribDict = {}
        attribDict['id'] = jsonobject[0]
        attribDict['specfile'] = jsonobject[1]
        attribDict['attributes'] = jsonobject[2]
        attribDict['params'] = [tuple(param) for param in jsonobject[3]]
        attribDict['scanListParams'] = [tuple(param) for param in jsonobject[4]]
        attribDict['scanList'] = jsonobject[5]
        attribDict['precursorList'] = jsonobject[6]
        attribDict['productList'] = jsonobject[7]
        for key, value in viewitems(attribDict):
            setattr(newInstance, key, value)
        return newInstance

    @staticmethod
    def jsonHook(encoded):
        if '__Smi__' in encoded:
            return Smi._fromJSON(encoded['__Smi__'])
        elif '__MzmlScan__' in encoded:
            return MzmlScan._fromJSON(encoded['__MzmlScan__'])
        elif '__MzmlProduct__' in encoded:
            return MzmlProduct._fromJSON(encoded['__MzmlProduct__'])
        elif '__MzmlPrecursor__' in encoded:
            return MzmlPrecursor._fromJSON(encoded['__MzmlPrecursor__'])
        else:
            return encoded


class Si(object):
    """Spectrum item (Si) - this is the spectrum representation intended to be used in maspy.
    A simplified representation of spectrum metadata. Contains only specifically imported attributes,
    which are necessary for data analysis. Does not follow any PSI datastructure or name space rules.

    Attributes can be transferred from the corresponding :class:`Smi` entry.

    identifier: The unique id of this spectrum, typically the scan number.
        Is used together with "specfile" as a key to access the spectrum in its container :class:`maspy.core.SiContainer`
        Should be derived from the spectrums nativeID format (MS:1000767)

    specfile: An id representing a group of spectra, typically of the same mzML file.
        Is used together with "identifier" as a key to access the spectrum in its container :class:`maspy.core.SiContainer`
    """
    def __init__(self, identifier, specfile):
        #super(SpectrumItem, self).__init__(identifier, specfile)
        self.id = identifier
        self.specfile = specfile
        self.isValid = None
        self.msLevel = None

    def _reprJSON(self):
        return {'__Si__': self.__dict__}

    @classmethod
    def _fromJSON(cls, jsonobject):
        newInstance = cls(None, None)
        newInstance.__dict__.update(jsonobject)
        return newInstance

    @staticmethod
    def jsonHook(encoded):
        if '__Si__' in encoded:
            return Si._fromJSON(encoded['__Si__'])
        else:
            return encoded


class MzmlScan(object):
    #TODO: docstring
    """
    :ivar scanWindowList: ... stored as a tuple because this variable is describing the measurement and should not be changed
    Note: the attributes "sourceFileRef" and "externalSpectrumID" are not supported
    TODO: add attributes "instrumentConfigurationRef" and "spectrumRef"
    """
    ## kwargs to only take the arguments needed and ignore additionally specified ones like 'arrayLength' of binaryDataArray
    __slots__ = ['scanWindowList', 'params']
    def __init__(self, scanWindowList=(), params=None, **kwargs):
        self.scanWindowList = tuple(tuple(_) for _ in scanWindowList)
        self.params = params

    def __getitem__(self, key):
        return getattr(self, key)

    #def __setitem__(self, key, value):
    #    setattr(self, key, value)

    def _reprJSON(self):
        return {'__MzmlScan__': (self.scanWindowList, self.params)}

    @classmethod
    def _fromJSON(cls, jsonobject):
        scanWindowList = _mzmlListAttribToTuple(jsonobject[0])
        params = [tuple(param) for param in jsonobject[1]]
        return cls(scanWindowList, params)


class MzmlProduct(object):
    #TODO: docstring
    __slots__ = ['isolationWindow']
    def __init__(self, isolationWindow=None, **kwargs):
        self.isolationWindow = tuple(isolationWindow)

    def __getitem__(self, key):
        return getattr(self, key)

    #def __setitem__(self, key, value):
    #    setattr(self, key, value)

    def _reprJSON(self):
        return {'__MzmlProduct__': self.isolationWindow}

    @classmethod
    def _fromJSON(cls, jsonobject):
        isolationWindow =[tuple(param) for param in jsonobject]
        return cls(isolationWindow)


class MzmlPrecursor(object):
    #TODO: docstring
    """
    Note: the attributes "sourceFileRef" and "externalSpectrumID" are not supported
    """
    __slots__ = ['spectrumRef', 'activation', 'isolationWindow', 'selectedIonList']
    def __init__(self, spectrumRef=None, activation=None, isolationWindow=None, selectedIonList=[], **kwargs):
        self.spectrumRef = spectrumRef
        self.isolationWindow = tuple(isolationWindow)
        self.selectedIonList = [tuple(_) for _ in selectedIonList]
        self.activation = tuple(activation)

    def __getitem__(self, key):
        return getattr(self, key)

    #def __setitem__(self, key, value):
    #    setattr(self, key, value)

    def _reprJSON(self):
        return {'__MzmlPrecursor__': (self.spectrumRef, self.activation, self.isolationWindow, self.selectedIonList)}

    @classmethod
    def _fromJSON(cls, jsonobject):
        spectrumRef = jsonobject[0]
        activation = [tuple(param) for param in jsonobject[1]]
        isolationWindow =[tuple(param) for param in jsonobject[2]]
        selectedIonList = _mzmlListAttribToTuple(jsonobject[3])
        return cls(spectrumRef, activation, isolationWindow, selectedIonList)


def _mzmlListAttribToTuple(oldList):
    """ Turns the list element, elements into tuples,
    Note: only intended for a list of elements that contain params, eg. mzml Node selectedIonList."""
    newList = list()
    for oldEntry in oldList:
        newEntry = [tuple(param) for param in oldEntry]
        newList.append(newEntry)
    return newList


def addMsrunContainers(mainContainer, subContainer):
    #TODO: docstrings
    #Note: does not generate new items, all items of the merged container still point to the inital memory location
    for specfile in subContainer.info:
        if specfile in mainContainer.info:
            continue

        mainContainer.addSpecfile(specfile, subContainer.info[specfile]['path'])
        for datatype, status in listitems(subContainer.info[specfile]['status']):
            if not status:
                continue
            datatypeContainer = datatype+'c'
            dataTypeContainer = getattr(mainContainer, datatypeContainer)
            subContainerData = getattr(subContainer, datatypeContainer)[specfile]
            dataTypeContainer[specfile] = subContainerData
            mainContainer.info[specfile]['status'][datatype] = True


##########################################################################
### SpectrumIdentificationItem related classes and functions #############
##########################################################################
class Sii(object):
    """Sii (SpectrumIdentificationItem)
    Representation of an ion fragment spectrum annotation, also referred to as
    peptide spectrum match (PSM).

    :ivar identifier: The unique id of this spectrum, typically the scan number.
        Is used together with "specfile" as a key to access this element in a
        :class:`maspy.core.SiiContainer` or the corresponding spectrum in a
        :class:`maspy.core.SiContainer`.

    :ivar specfile: An id representing a group of spectra, typically of the same
        mzML file. Is used together with "identifier" as a unique key.

    :ivar rank: The rank of this Sii compared to others for the same MSn
        spectrum. The rank is based on a score defined in the SiiContainer. If
        multiple Sii have the same top score, they should all be assigned
        rank=1.

    :ivar isValid: bool or None if not specified
        this attribute can be used to flag if a Sii has passed a given quality
        threshold or been validated as correct. Is used to filter valid elements
        of :class:`maspy.core.SiContainer`.
    """
    def __init__(self, identifier, specfile):
        self.id = identifier
        self.specfile = specfile
        self.rank = None
        self.isValid = None

    def _reprJSON(self):
        return {'__Sii__': self.__dict__}

    @classmethod
    def _fromJSON(cls, jsonobject):
        newInstance = cls(None, None)
        newInstance.__dict__.update(jsonobject)
        return newInstance

    @staticmethod
    def jsonHook(encoded):
        if '__Sii__' in encoded:
            return Sii._fromJSON(encoded['__Sii__'])
        else:
            return encoded


class SiiContainer(object):
    """ItemContainer for MSn spectrum identifications (Peptide Spectrum Matches),
    SiiContainer = Spectrum Identification Item Container.

    :ivar container: Access :class:`ContainerItem` storage list via a specfile keyword: {specfile:[ContainerItem(), ContainerItem(), ...]}
    :ivar info: a dictionary containing information about the imported specfiles;
        key = specfilename, value = {"scoreAttr": str, "largerBetter": bool, "path": str}
        "scoreAttr" describes which :class:`Sii` attribute should be used for ranking.
        "largerBetter" specifies wheter a larger score signifies a better match.
        "path" contains a directory path used for saving and loading

    #Note: In the future this container may be integrated in an evidence or mzIdentML like container.
    """
    def __init__(self):
        self.container = dict()
        self.info = dict()

    def getArrays(self, attr=None, specfiles=None, sort=False, reverse=False,
                  selector=lambda sii: sii.isValid, defaultValue=None):
        """Return a condensed array of data selected from :class:`Sii` objects of :instance:`self.container`
        for fast and convenient data processing.

        :param attr: list of :class:`Sii` item attributes that should be added to the returned array.
            The attributes "id" and "specfile" are always included, in combination they serve as a unique id.
        :param defaultValue: if an item is missing an attribute, the "defaultValue" is
            added to the array instead.
        :param specfiles: filenames of msrun files - if specified return only items from those files
        :type specfiles: str or [str, str, ...]
        :param sort: if "sort" is specified the returned list of items is sorted according to the :class:`Si`
            attribute specified by "sort", if the attribute is not present the item is skipped.
        :param reverse: boolean to reverse sort order
        :param selector: a function which is called with each :class:`Si` item and returns
            True (include item) or False (discard item). If not specified all items are returned.
            By default only items with "isValid" == True are returned.

        return {'attribute1': numpy.array(), 'attribute1': numpy.array(), ...}
        """
        attr = attr if attr is not None else []
        attr = set(['id', 'specfile'] + aux.toList(attr))
        items = self.getItems(specfiles, sort, reverse, selector)
        return _getArrays(items, attr, defaultValue)

    def getItems(self, specfiles=None, sort=False, reverse=False, selector=lambda sii: sii.isValid):
        """Generator that yields filtered and/or sorted :class:`Sii` objects from :instance:`self.container`

        :param specfiles: filenames of msrun files - if specified return only items from those files
        :type specfiles: str or [str, str, ...]
        :param sort: if "sort" is specified the returned list of items is sorted according to the :class:`Si`
        attribute specified by "sort", if the attribute is not present the item is skipped.
        :param reverse: boolean to reverse sort order
        :param selector: a function which is called with each :class:`Si` item and returns
        True (include item) or False (discard item). If not specified all items are returned
        """
        specfiles = [_ for _ in viewkeys(self.info)] if specfiles is None else aux.toList(specfiles)
        return _getListItems(self.container, specfiles, sort, reverse, selector)

    def getValidItem(self, specfile, identifier):
        """Returns a valid item or None,
        assumes that self.container[specfile][identifier] is a sorted list"""
        for item in self.container[specfile][identifier]:
            if item.isValid:
                return item
        else:
            return None

    def addSpecfile(self, specfiles, path):
        """Adds specfile entries to self.info, but doesn't import any data yet.
        To actually import the SiiContainer files, use :func:`siiContainer.load()`.
        """
        for specfile in aux.toList(specfiles):
            if specfile not in self.info:
                self._addSpecfile(specfile, path)
            else:
                warnings.warn('Specfile is already present in the SiiContainer.\nname: %s \npath: %s' %(specfile, self.info[specfile]['path']))

    def _addSpecfile(self, specfile, path):
        """Adds a new specfile entry to SiiContainer.info. """
        self.info[specfile] = {'scoreAttr': None, 'largerBetter': None, 'path': path}
        self.container[specfile] = dict()

    def setPath(self, folderpath, specfiles=None):
        """Change the folderpath of the specified specfiles. If save is called and no "siic"
        (SpectrumIdentificationItemContainer) file is present in the specified path,
        a new file is generated.
        """
        specfiles = [_ for _ in viewkeys(self.info)] if specfiles is None else specfiles
        _containerSetPath(self, folderpath, specfiles)

    def removeSpecfile(self, specfiles):
        """Completely removes the specified specfiles from the SiiContainer."""
        for specfile in aux.toList(specfiles):
            del self.container[specfile]
            del self.info[specfile]

    def save(self, specfiles=None, compress=True, path=None):
        #TODO: docstring
        specfiles = [_ for _ in viewkeys(self.info)] if specfiles is None else specfiles
        for specfile in aux.toList(specfiles):
            if specfile not in self.info:
                print('Error while saving', specfile, ', not found in msrunCountainer!')
                continue
            else:
                specfilePath = self.info[specfile]['path'] if path is None else path

            with aux.PartiallySafeReplace() as msr:
                filename = specfile + '.siic'
                filepath = aux.joinpath(specfilePath, filename)
                with msr.open(filepath, mode='w+b') as openfile:
                    self._writeContainer(openfile, specfile, compress=compress)

    def _writeContainer(self, filelike, specfile, compress=True):
        #TODO: docstring
        aux.writeJsonZipfile(filelike, self.container[specfile], compress=compress)
        zipcomp = zipfile.ZIP_DEFLATED if compress else zipfile.ZIP_STORED
        with zipfile.ZipFile(filelike, 'a', allowZip64=True) as containerFile:
            infodata = {key: value for key, value in viewitems(self.info[specfile]) if key != 'path'}
            containerFile.writestr('info', json.dumps(infodata, zipcomp))

    def load(self, specfiles=None):
        """Import specfiles"""
        #TODO: docstring
        specfiles = [_ for _ in viewkeys(self.info)] if specfiles is None else specfiles
        for specfile in aux.toList(specfiles):
            if specfile not in self.info:
                print(''.join(('Error while loading "', specfile, '", not present in msrunCountainer.info')))
                continue
            else:
                siiPath = aux.joinpath(self.info[specfile]['path'], specfile+'.siic')
                with zipfile.ZipFile(siiPath, 'r') as containerZip:
                    #Convert the zipfile data into a str object, necessary since containerZip.read() returns a bytes object.
                    jsonString = io.TextIOWrapper(containerZip.open('data'), encoding='utf-8').read()
                    infoString = io.TextIOWrapper(containerZip.open('info'), encoding='utf-8').read()
                self.container[specfile] = json.loads(jsonString, object_hook=Sii.jsonHook)
                self.info[specfile].update(json.loads(infoString))

    def addSiInfo(self, msrunContainer, specfiles=None, attributes=['obsMz', 'rt', 'charge']):
        #TODO: use obsMz (observed)? or mz? in contrary to exMz (exact mz) or calcMz (calculated mz)
        """ Copy attributes into Sii from the corresponding Si in msrunContainer,
        if an attribute is not presend in the SpectrumItem the attribute value is set to None
        Attribute examples: 'obsMz', 'rt', 'charge', 'tic', 'iit', 'ms1Id'
        """
        specfiles = [_ for _ in viewkeys(self.info)] if specfiles == None else aux.toList(specfiles)

        for specfile in specfiles:
            if specfile not in self.info:
                print(specfile, 'not present in siiContainer.')
            elif specfile not in msrunContainer.info:
                print(specfile, 'not present in msrunContainer.')
            else:
                for identifier in self.container[specfile]:
                    si = msrunContainer.sic[specfile][identifier]
                    for sii in self.container[specfile][identifier]:
                        for attribute in attributes:
                            setattr(sii, attribute, getattr(si, attribute, None))

    def calcMz(self, specfiles=None, guessCharge=True, obsMzKey='mz'):
        #TODO: docstring
        # Guess charge uses the calculated mass and the observed m/z value to calculate the charge
        specfiles = [_ for _ in viewkeys(self.info)] if specfiles is None else aux.toList(specfiles)
        tempPeptideMasses = dict()
        for specfile in specfiles:
            if specfile not in self.info:
                print(specfile, 'not present in siiContainer.')
            else:
                for sii in self.getItems(specfiles=specfile):
                    charge = sii.charge
                    peptide = sii.peptide
                    if peptide not in tempPeptideMasses:
                        if hasattr(sii, 'diPeptide'):
                            tempPeptideMasses[peptide] = (maspy.peptidemethods.calcPeptideMass(sii.peptide1) +
                                                          maspy.peptidemethods.calcPeptideMass(sii.peptide2))
                        else:
                            tempPeptideMasses[peptide] = maspy.peptidemethods.calcPeptideMass(peptide)
                    peptideMass = tempPeptideMasses[peptide]
                    if charge is not None:
                        sii.calcMz = maspy.peptidemethods.calcMzFromMass(peptideMass, charge)
                    elif guessCharge:
                        guessedCharge = round(peptideMass / (getattr(sii, obsMzKey) - maspy.constants.atomicMassProton), 0)
                        sii.calcMz = maspy.peptidemethods.calcMzFromMass(peptideMass, guessedCharge)
                        sii.charge = guessedCharge
        del(tempPeptideMasses)


###########################################################
### FeatureItem related classes and functions #############
###########################################################
class Fi(object):
    """FeatureItem (Fi), representation of a peptide LC-MS feature.

    :ivar containerId: used to look up item in :attr:`ItemContainer.index`
    :ivar id: identifier in original file
    :ivar specfile: Keyword (filename) to represent the originating file
    :ivar isValid: this attribute can be used to filter data.
    Should be set to True or False, None if unspecified
    :ivar isMatched: None if unspecified, should be set to False on import, True if any Si or Sii elements could be matched
    :ivar isAnnotated: None if unspecified, should be set to False on import, True if any Sii elements could be matched
    :ivar siIds: tuple(specfile, id) of matched Si
    :ivar siiIds: tuple(specfile, id) of matched Sii
    :ivar peptide: peptide sequence of the best scoring Sii match
    :ivar sequence: plain amino acid sequence of best scoring Sii match, used to retrieve protein information
    :ivar score: score of best scoring Sii match
    """
    def __init__(self, identifier, specfile):
        self.id = identifier
        self.specfile = specfile
        self.isValid = None

        #Annotation information
        self.isMatched = None
        self.isAnnotated = None
        self.siIds = list()
        self.siiIds = list()
        self.peptide = None
        self.sequence = None
        self.score = None

    def _reprJSON(self):
        return {'__Fi__': self.__dict__}

    @classmethod
    def _fromJSON(cls, jsonobject):
        newInstance = cls(None, None)
        newInstance.__dict__.update(jsonobject)
        newInstance.siIds = [tuple(_) for _ in newInstance.siIds]
        newInstance.siiIds = [tuple(_) for _ in newInstance.siiIds]
        return newInstance

    @staticmethod
    def jsonHook(encoded):
        if '__Fi__' in encoded:
            return Fi._fromJSON(encoded['__Fi__'])
        else:
            return encoded


class FiContainer(object):
    """ItemContainer for peptide elution features :class`Fi` (FeatureItem).

    #TODO: docstring
    """
    def __init__(self):
        self.container = dict()
        self.info = dict()

    def getArrays(self, attr=None, specfiles=None, sort=False, reverse=False,
                  selector=lambda fi: fi.isValid, defaultValue=None):
        """Return a condensed array of data selected from :class:`Fi` objects of :instance:`self.container`
        for fast and convenient data processing.

        :param attr: list of :class:`Fi` item attributes that should be added to the returned array.
            The attributes "id" and "specfile" are always included, in combination they serve as a unique id.
        :param defaultValue: if an item is missing an attribute, the "defaultValue" is
            added to the array instead.
        :param specfiles: filenames of msrun files - if specified return only items from those files
        :type specfiles: str or [str, str, ...]
        :param sort: if "sort" is specified the returned list of items is sorted according to the :class:`Si`
            attribute specified by "sort", if the attribute is not present the item is skipped.
        :param reverse: boolean to reverse sort order
        :param selector: a function which is called with each :class:`Si` item and returns
            True (include item) or False (discard item). If not specified all items are returned.
            By default only items with "isValid" == True are returned.

        return {'attribute1': numpy.array(), 'attribute1': numpy.array(), ...}
        """
        attr = attr if attr is not None else []
        attr = set(['id', 'specfile'] + aux.toList(attr))
        items = self.getItems(specfiles, sort, reverse, selector)
        return _getArrays(items, attr, defaultValue)

    def getItems(self, specfiles=None, sort=False, reverse=False, selector=lambda fi: fi.isValid):
        """Generator that yields filtered and/or sorted :class:`Fi` objects from :instance:`self.container`

        :param specfiles: filenames of msrun files - if specified return only items from those files
        :type specfiles: str or [str, str, ...]
        :param sort: if "sort" is specified the returned list of items is sorted according to the :class:`Fi`
            attribute specified by "sort", if the attribute is not present the item is skipped.
        :param reverse: boolean to reverse sort order
        :param selector: a function which is called with each :class:`Si` item and returns
            True (include item) or False (discard item). If not specified all items are returned
        """
        specfiles = [_ for _ in viewkeys(self.info)] if specfiles is None else aux.toList(specfiles)
        return _getItems(self.container, specfiles, sort, reverse, selector)

    def addSpecfile(self, specfiles, path):
        """Adds specfile entries to self.info, but doesn't import any data yet.
        To actually import the FiContainer files, use :func:`msrunContainer.load()`.
        """
        for specfile in aux.toList(specfiles):
            if specfile not in self.info:
                self._addSpecfile(specfile, path)
            else:
                warnings.warn('Specfile is already present in the FiContainer.\nname: %s \npath: %s' %(specfile, self.info[specfile]['path']))

    def _addSpecfile(self, specfile, path):
        """Adds a new specfile entry to FiContainer.info. """
        self.info[specfile] = {'path': path}
        self.container[specfile] = dict()

    def setPath(self, folderpath, specfiles=None):
        """Change the folderpath of the specified specfiles. If save is called and no "fic"
        (FeatureItemContainer) file is present in the specified path, a new file is generated.
        """
        specfiles = [_ for _ in viewkeys(self.info)] if specfiles is None else specfiles
        _containerSetPath(self, folderpath, specfiles)

    def removeSpecfile(self, specfiles):
        """Completely removes the specified specfiles from the SiiContainer."""
        for specfile in aux.toList(specfiles):
            del self.container[specfile]
            del self.info[specfile]

    def save(self, specfiles=None, compress=True, path=None):
        #TODO: docstring
        specfiles = [_ for _ in viewkeys(self.info)] if specfiles is None else specfiles
        for specfile in aux.toList(specfiles):
            if specfile not in self.info:
                print('Error while saving', specfile, ', not found in msrunCountainer!')
                continue
            else:
                specfilePath = self.info[specfile]['path'] if path is None else path

            with aux.PartiallySafeReplace() as msr:
                filename = specfile + '.fic'
                filepath = aux.joinpath(specfilePath, filename)
                with msr.open(filepath) as openfile:
                    self._writeContainer(openfile, specfile, compress=compress)

    def _writeContainer(self, filelike, specfile, compress=True):
        aux.writeJsonZipfile(filelike, self.container[specfile], compress=compress)
        #zipcomp = zipfile.ZIP_DEFLATED if compress else zipfile.ZIP_STORED
        #with zipfile.ZipFile(filelike, 'a', allowZip64=True) as containerFile:
        #    infodata = {key: value for key, value in viewitems(self.info[specfile]) if key != 'path'}
        #    containerFile.writestr('info', json.dumps(infodata, zipcomp))

    def load(self, specfiles=None):
        """Import specfiles"""
        specfiles = [_ for _ in viewkeys(self.info)] if specfiles is None else specfiles
        for specfile in aux.toList(specfiles):
            if specfile not in self.info:
                print(''.join(('Error while loading "', specfile, '", not present in msrunCountainer.info')))
                continue
            else:
                fiPath = aux.joinpath(self.info[specfile]['path'], specfile+'.fic')
                with zipfile.ZipFile(fiPath, 'r') as containerZip:
                    #Convert the zipfile data into a str object, necessary since containerZip.read() returns a bytes object.
                    jsonString = io.TextIOWrapper(containerZip.open('data'), encoding='utf-8').read()
                    #infoString = io.TextIOWrapper(containerZip.open('info'), encoding='utf-8').read()
                self.container[specfile] = json.loads(jsonString, object_hook=Fi.jsonHook)
                #self.info[specfile].update(json.loads(infoString))

    def removeAnnotation(self, specfiles=None):
        """Remove all annotation information from :class:`FeatureItem` in :class:`FeatureContainer`."""
        specfiles = [_ for _ in viewkeys(self.info)] if specfiles is None else specfiles
        for specfile in aux.toList(specfiles):
            for item in viewvalues(self.container[specfile]):
                item.isMatched = False
                item.isAnnotated = False
                item.siIds = list()
                item.siiIds = list()
                item.peptide = None
                item.sequence = None
                item.score = None
