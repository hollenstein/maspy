"""
The core module contains python classes to representat spectra, peptide spectrum
matches and peptide LC-MS features, and containers which manage storage,
data access, saving and loading of these data types. 
"""
######################### Python 2 and 3 compatibility #########################
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
from collections import defaultdict as ddict
import io
import json
from operator import itemgetter as ITEMGETTER
import os
import warnings

from lxml import etree as ETREE
import numpy
import zipfile

import maspy.auxiliary as aux
import maspy.constants
import maspy.peptidemethods


############################################
### Common container functions #############
############################################
def _getArrays(items, attr, defaultValue):
    """Return arrays with equal size of item attributes from a list of sorted
    "items" for fast and convenient data processing.

    :param attr: list of item attributes that should be added to the returned
        array.
    :param defaultValue: if an item is missing an attribute, the "defaultValue"
        is added to the array instead.

    :returns: {'attribute1': numpy.array([attributeValue1, ...]), ...}
    """
    arrays = dict([(key, []) for key in attr])
    for item in items:
        for key in attr:
            arrays[key].append(getattr(item, key, defaultValue))
    for key in [_ for _ in viewkeys(arrays)]:
        arrays[key] = numpy.array(arrays[key])

    return arrays


def _getItems(container, containerKeys=None, sort=False, reverse=False,
              selector=lambda item: True):
    """Generator that yields filtered and/or sorted items from the specified
    "container".

    :param container: The container has to be a dictionary of dictionaries that
        contain some kind of items. Depending on the specified parameters all or
        a subset of these items are yielded.
        ``{containerKey1: {key1: item1, key2: item2, ...}, ...}``
    :param containerKeys: valid keys of the "container", if None all keys are
        considered.
    :type containerKeys: a single dictionary key or a list of keys
    :param sort: if "sort" is specified the returned list of items is sorted
        according to the item attribute specified by "sort", if the attribute is
        not present the item is skipped.
    :param reverse: bool, ``True`` reverses the sort order
    :param selector: a function which is called with each item and returns
        True (include item) or False (discard item). If not specified all items
        are returned

    :returns: items from container that passed the selector function
    """
    if containerKeys is None:
        containerKeys = [_ for _ in viewkeys(container)]
    else:
        containerKeys = aux.toList(containerKeys)

    if sort:
        sortIdentifier = list()
        for containerKey in containerKeys:
            for identifier in [_ for _ in viewkeys(container[containerKey])]:
                item = container[containerKey][identifier]
                if selector(item):
                    try:
                        sortIdentifier.append((getattr(item, sort),
                                               containerKey, identifier
                                               )
                                              )
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


def _getListItems(container, containerKeys=None, sort=False, reverse=False,
                  selector=None):
    """Generator that yields filtered and/or sorted items from the specified
    container.

    .. note::
        Use this function instead of :class:`_getItems()` if the item entries in
        the container is not the item itself but a list of items, see param
        container below.

    :param container: The container has to be a dictionary of dictionaries that
        contain a list of some kind of items. Depending on the specified
        parameters all or a subset of these items are yielded.
        ``{containerKey1: {key1: [item1, ...], key2: [item1, ...], ...}, ...}``
    :param containerKeys: valid keys of the "container", if None all keys are
        considered.
    :type containerKeys: a single dictionary key or a list of keys
    :param sort: if "sort" is specified the returned list of items is sorted
        according to the item attribute specified by "sort", if the attribute is
        not present the item is skipped.
    :param reverse: bool, ``True`` reverses the sort order
    :param selector: a function which is called with each item and returns
        True (include item) or False (discard item). If not specified all items
        are returned

    :returns: items from container that passed the selector function
    """
    selector = (lambda item: True) if selector is None else selector
    if containerKeys is None:
        containerKeys = [_ for _ in viewkeys(container)]
    else:
        containerKeys = aux.toList(containerKeys)

    if sort:
        sortIdentifier = list()
        for containerKey in containerKeys:
            for identifier in [_ for _ in viewkeys(container[containerKey])]:
                for itemPos, item in enumerate(container[containerKey][identifier]):
                    if selector(item):
                        try:
                            sortIdentifier.append((getattr(item, sort),
                                                   containerKey, identifier,
                                                   itemPos
                                                   )
                                                  )

                        except AttributeError:
                            pass
        sortIdentifier.sort(key=ITEMGETTER(0), reverse=reverse)
        for _, containerKey, identifier, itemPos in sortIdentifier:
            yield container[containerKey][identifier][itemPos]
    else:
        for containerKey in containerKeys:
            for identifier in [_ for _ in viewkeys(container[containerKey])]:
                for item in container[containerKey][identifier]:
                    if selector(item):
                        yield item


def _containerSetPath(container, folderpath, specfiles):
    """Helper function for :class:`MsrunContainer`, :class:`SiiContainer` and
    :class:`FiContainer`. Changes the folderpath of the specified specfiles in
    container.info: ``container.info[specfile]['path'] = folderpath``.

    :param container: a container like class that has an attribute ``.info``
    :param folderpath: a filedirectory
    :param specfiles: a list of ms-run names
    """
    if not os.path.exists(folderpath):
        warntext = 'Error while calling "_containerSetPath()": The specified '\
                   'directory "%s" does not exist!' %(folderpath, )
        warnings.warn(warntext)
    for specfile in specfiles:
        if specfile in container.info:
            container.info[specfile]['path'] = folderpath
        else:
            warntext = 'Error while calling "_containerSetPath()": The '\
                       'specfile "%s" is not present in the container!'\
                       %(specfile, )
            warnings.warn(warntext)


##############################################################
### MsrunContainer related classes and functions #############
##############################################################
class MsrunContainer(object):
    """Container for mass spectrometry data (eg MS1 and MS2 spectra), provides
    full support for mzML files, see `mzML schema
    documentation <http://www.peptideatlas.org/tmp/mzML1.1.0.html>`_.

    :ivar rmc: "run metadata container", contains mzML metadata elements from
        the mzML file as a ``lxml.etree.Element`` object. This comprises all
        ``mzML`` subelements, except for the `run`` element subelements
        ``spectrumList`` and ``chromatogramList``.
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
        datatypes are currently imported.

        code example::

            {u'JD_06232014_sample1_A': {u'path': u'C:/filedirectory',
                                        u'status': {u'ci': True, u'rm': True,
                                                    u'sai': True, u'si': True,
                                                    u'smi': True
                                                    }
                                        }
             }

    .. note::
        The structure of the containers ``rmc``, ``cic``, ``smic``, ``saic``
        and ``sic`` is always: ``{"specfilename": {"itemId": item, ...}, ...}``

    """
    def __init__(self):
        self.rmc = {}
        self.cic = {}
        self.smic = {}
        self.saic = {}
        self.sic = {}
        self.info = {}

    def getArrays(self, attr=None, specfiles=None, sort=False, reverse=False,
                  selector=None, defaultValue=None):
        """Return a condensed array of data selected from :class:`Si` instances
        from ``self.sic`` for fast and convenient data processing.

        :param attr: list of :class:`Si` item attributes that should be added to
            the returned array. The attributes "id" and "specfile" are always
            included, in combination they serve as a unique id.
        :param defaultValue: if an item is missing an attribute, the
            "defaultValue" is added to the array instead.
        :param specfiles: filenames of ms-run files, if specified return only
            items from those files
        :type specfiles: str or [str, str, ...]
        :param sort: if "sort" is specified the returned list of items is sorted
            according to the :class:`Si` attribute specified by "sort", if the
            attribute is not present the item is skipped.
        :param reverse: bool, set True to reverse sort order
        :param selector: a function which is called with each :class:`Si` item
            and has to return True (include item) or False (discard item).
            Default function is: ``lambda si: True``

        :returns: {'attribute1': numpy.array(),
                   'attribute2': numpy.array(),
                   ...
                   }
        """
        selector = (lambda si: True) if selector is None else selector
        attr = attr if attr is not None else []
        attr = set(['id', 'specfile'] + aux.toList(attr))
        items = self.getItems(specfiles, sort, reverse, selector)
        return _getArrays(items, attr, defaultValue)

    def getItems(self, specfiles=None, sort=False, reverse=False,
                 selector=None):
        """Generator that yields filtered and/or sorted :class:`Si` instances
        from ``self.sic``.

        :param specfiles: filenames of ms-run files - if specified return only
            items from those files
        :type specfiles: str or [str, str, ...]
        :param sort: if "sort" is specified the returned list of items is sorted
            according to the :class:`Si` attribute specified by "sort", if the
            attribute is not present the item is skipped.
        :param reverse: bool, ``True`` reverses the sort order
        :param selector: a function which is called with each ``Si`` item
            and returns True (include item) or False (discard item). Default
            function is: ``lambda si: True``

        :returns: items from container that passed the selector function
        """
        selector = (lambda si: True) if selector is None else selector
        if specfiles is None:
            specfiles = [_ for _ in viewkeys(self.info)]
        else:
            specfiles = aux.toList(specfiles)
        return _getItems(self.sic, specfiles, sort, reverse, selector)

    def getItem(self, specfile, identifier):
        """Returns a :class:`Si` instance from ``self.sic``.

        :param specfile: a ms-run file name
        :param identifier: item identifier ``Si.id``

        :returns: ``self.sic[specfile][identifier]``
        """
        return self.sic[specfile][identifier]

    def addSpecfile(self, specfiles, path):
        """Prepares the container for loading ``mrc`` files by adding specfile
        entries to ``self.info``. Use :func:`MsrunContainer.load()` afterwards
        to actually import the files

        :param specfiles: the name of an ms-run file or a list of names
        :type specfiles: str or [str, str, ...]
        :param path: filedirectory used for loading and saving ``mrc`` files
        """
        for specfile in aux.toList(specfiles):
            if specfile not in self.info:
                self._addSpecfile(specfile, path)
            else:
                warntext = 'Error while calling "MsrunContainer.addSpecfile()"'\
                           ': "%s" is already present "MsrunContainer.info"'\
                            % (specfile, )
                warnings.warn(warntext)

    def _addSpecfile(self, specfile, path):
        """Adds a new specfile entry to MsrunContainer.info. See also
        :class:`MsrunContainer.addSpecfile()`.

        :param specfile: the name of an ms-run file
        :param path: filedirectory used for loading and saving ``mrc`` files
        """
        datatypeStatus = {'rm': False, 'ci': False, 'smi': False, 'sai': False,
                          'si': False
                          }
        self.info[specfile] = {'path': path, 'status': datatypeStatus}

    def setPath(self, folderpath, specfiles=None):
        """Changes the folderpath of the specified specfiles. The folderpath is
        used for saving and loading of ``mrc`` files.

        :param specfiles: the name of an ms-run file or a list of names. If None
            all specfiles are selected.
        :type specfiles: None, str, [str, str]
        :param folderpath: a filedirectory
        """
        if specfiles is None:
            specfiles = [_ for _ in viewkeys(self.info)]
        else:
            specfiles = aux.toList(specfiles)

        _containerSetPath(self, folderpath, specfiles)

    def removeData(self, specfiles=None, rm=False, ci=False, smi=False,
                   sai=False, si=False):
        """Removes the specified datatypes of the specfiles from the
        msrunContainer. To completely remove a specfile use
        :func:`MsrunContainer.removeSpecfile`, which also removes the complete
        entry from ``self.info``.

        :param specfiles: the name of an ms-run file or a list of names. If None
            all specfiles are selected.
        :type specfiles: None, str, [str, str]
        :param rm: bool, True to select ``self.rmc``
        :param ci: bool, True to select ``self.cic``
        :param smi: bool, True to select ``self.smic``
        :param sai: bool, True to select ``self.saic``
        :param si: bool, True to select ``self.sic``
        """
        if specfiles is None:
            specfiles = [_ for _ in viewkeys(self.info)]
        else:
            specfiles = aux.toList(specfiles)
        #TODO: add a check if specfiles are present in the container
        typeToContainer = {'rm': 'rmc', 'ci': 'cic', 'smi': 'smic',
                           'sai': 'saic', 'si': 'sic'
                           }
        datatypes = self._processDatatypes(rm, ci, smi, sai, si)
        for specfile in specfiles:
            for datatype in datatypes:
                datatypeContainer = typeToContainer[datatype]
                dataContainer = getattr(self, datatypeContainer)
                try:
                    del dataContainer[specfile]
                except KeyError:
                    pass
                finally:
                    self.info[specfile]['status'][datatype] = False

    def removeSpecfile(self, specfiles):
        """Completely removes the specified specfiles from the
        ``msrunContainer``.

        :param specfiles: the name of an ms-run file or a list of names. If None
            all specfiles are selected.
        :type specfiles: str, [str, str]
        """
        for specfile in aux.toList(specfiles):
            for datatypeContainer in ['rmc', 'cic', 'smic', 'saic', 'sic']:
                dataContainer = getattr(self, datatypeContainer)
                try:
                    del dataContainer[specfile]
                except KeyError:
                    pass
            del self.info[specfile]

    def _processDatatypes(self, rm, ci, smi, sai, si):
        """Helper function that returns a list of datatype strings, depending
        on the parameters boolean value.

        :param rm: bool, True to add ``rm``
        :param ci: bool, True to add ``ci``
        :param smi: bool, True to add ``smi``
        :param sai: bool, True to add ``sai``
        :param si: bool, True to add ``si``

        :returns: [datatype1, ...]
        """
        datatypes = list()
        for datatype, value in [('rm', rm), ('ci', ci), ('smi', smi),
                                ('sai', sai), ('si', si)]:
            if value:
                datatypes.append(datatype)
        return datatypes

    def save(self, specfiles=None, rm=False, ci=False, smi=False, sai=False,
             si=False, compress=True, path=None):
        """Writes the specified datatypes to ``mrc`` files on the hard disk.

        .. note::
            If ``.save()`` is called and no ``mrc`` files are present in the
            specified path new files are generated, otherwise old files are
            replaced.

        :param specfiles: the name of an ms-run file or a list of names. If None
            all specfiles are selected.
        :type specfiles: None, str, [str, str]
        :param rm: bool, True to select ``self.rmc`` (run metadata)
        :param ci: bool, True to select ``self.cic`` (chromatogram items)
        :param smi: bool, True to select ``self.smic`` (spectrum metadata items)
        :param sai: bool, True to select ``self.saic`` (spectrum array items)
        :param si: bool, True to select ``self.sic`` (spectrum items)
        :param compress: bool, True to use zip file compression
        :param path: filedirectory to which the ``mrc`` files are written. By
            default the parameter is set to ``None`` and the filedirectory is
            read from ``self.info[specfile]['path']``
        """
        if specfiles is None:
            specfiles = [_ for _ in viewkeys(self.info)]
        else:
            specfiles = aux.toList(specfiles)
        datatypes = self._processDatatypes(rm, ci, smi, sai, si)
        if len(datatypes) == 0:
            datatypes = ['rm', 'ci', 'smi', 'sai', 'si']

        for specfile in specfiles:
            if specfile not in self.info:
                warntext = 'Error while calling "MsrunContainer.save()": "%s" '\
                           'is not present in "MsrunContainer.info"!'\
                            % (specfile, )
                warnings.warn(warntext)
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
                           self._writeRmc(openfile, specfile)
                        elif datatype == 'ci':
                           self._writeCic(openfile, specfile, compress)
                        elif datatype == 'si':
                           self._writeSic(openfile, specfile, compress)
                        elif datatype == 'smi':
                           self._writeSmic(openfile, specfile, compress)
                        elif datatype == 'sai':
                           self._writeSaic(openfile, specfile, compress)

    def _writeCic(self, filelike, specfile, compress):
        """Writes the ``.cic`` container entry of the specified specfile to the
        ``mrc_cic`` format. For details see
         :func:`maspy.auxiliary.writeBinaryItemContainer()`

        :param filelike: path to a file (str) or a file-like object
        :param specfile: name of an ms-run file present in ``self.info``
        :param compress: bool, True to use zip file compression
        """
        aux.writeBinaryItemContainer(filelike, self.cic[specfile], compress)

    def _writeSaic(self, filelike, specfile, compress):
        """Writes the ``.ssic`` container entry of the specified specfile to the
        ``mrc_saic`` format. For details see
         :func:`maspy.auxiliary.writeBinaryItemContainer()`

        :param filelike: path to a file (str) or a file-like object
        :param specfile: name of an ms-run file present in ``self.info``
        :param compress: bool, True to use zip file compression
        """
        aux.writeBinaryItemContainer(filelike, self.saic[specfile], compress)

    def _writeSmic(self, filelike, specfile, compress):
        """Writes the ``.smic`` container entry of the specified specfile to the
        ``mrc_smic`` format. For details see
         :func:`maspy.auxiliary.writeJsonZipfile()`

        :param filelike: path to a file (str) or a file-like object
        :param specfile: name of an ms-run file present in ``self.info``
        :param compress: bool, True to use zip file compression
        """
        aux.writeJsonZipfile(filelike, self.smic[specfile], compress)

    def _writeSic(self, filelike, specfile, compress):
        """Writes the ``.sic`` container entry of the specified specfile to the
        ``mrc_sic`` format. For details see
         :func:`maspy.auxiliary.writeJsonZipfile()`

        :param filelike: path to a file (str) or a file-like object
        :param specfile: name of an ms-run file present in ``self.info``
        :param compress: bool, True to use zip file compression
        """
        aux.writeJsonZipfile(filelike, self.sic[specfile], compress)

    def _writeRmc(self, filelike, specfile):
        """Writes the ``.rmc`` container entry of the specified specfile as an
        human readable and pretty formatted xml string.

        :param filelike: path to a file (str) or a file-like object
        :param specfile: name of an ms-run file present in ``self.info``
        """
        xmlString = ETREE.tostring(self.rmc[specfile], pretty_print=True)
        filelike.write(xmlString)

    def load(self, specfiles=None, rm=False, ci=False, smi=False, sai=False,
             si=False):
        """Import the specified datatypes from ``mrc`` files on the hard disk.

        :param specfiles: the name of an ms-run file or a list of names. If None
            all specfiles are selected.
        :type specfiles: None, str, [str, str]
        :param rm: bool, True to import ``mrc_rm`` (run metadata)
        :param ci: bool, True to import ``mrc_ci`` (chromatogram items)
        :param smi: bool, True to import ``mrc_smi`` (spectrum metadata items)
        :param sai: bool, True to import ``mrc_sai`` (spectrum array items)
        :param si: bool, True to import ``mrc_si`` (spectrum items)
        """
        if specfiles is None:
            specfiles = [_ for _ in viewkeys(self.info)]
        else:
            specfiles = aux.toList(specfiles)

        #Select only specfiles which are present in the ``self.info``.
        selectedSpecfiles = list()
        for specfile in specfiles:
            if specfile not in self.info:
                warntext = 'Error while calling "MsrunContainer.load()": "%s" '\
                           'not present in MsrunContainer.info' % specfile
                warnings.warn(warntext)
            else:
                selectedSpecfiles.append(specfile)

        datatypes = self._processDatatypes(rm, ci, smi, sai, si)
        if len(datatypes) == 0:
            datatypes = ['rm', 'ci', 'smi', 'sai', 'si']

        for specfile in selectedSpecfiles:
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
                self.cic[specfile] = aux.loadBinaryItemContainer(ciPath,
                                                                 Ci.jsonHook)
                msrunInfo['status']['ci'] = True

            if 'smi' in datatypes:
                smiPath = aux.joinpath(specfilePath, specfile+'.mrc_smi')
                with zipfile.ZipFile(smiPath, 'r') as containerZip:
                    #Convert the zipfile data into a str object,necessary since
                    #containerZip.read() returns a bytes object.
                    jsonString = io.TextIOWrapper(containerZip.open('data'),
                                                  encoding='utf-8'
                                                  ).read()
                self.smic[specfile] = json.loads(jsonString,
                                                 object_hook=Smi.jsonHook
                                                 )
                msrunInfo['status']['smi'] = True

            if 'sai' in datatypes:
                saiPath = aux.joinpath(specfilePath, specfile+'.mrc_sai')
                self.saic[specfile] = aux.loadBinaryItemContainer(saiPath,
                                                                  Sai.jsonHook
                                                                  )
                msrunInfo['status']['sai'] = True

            if 'si' in datatypes:
                siPath = aux.joinpath(specfilePath, specfile+'.mrc_si')
                with zipfile.ZipFile(siPath, 'r') as containerZip:
                    #Convert the zipfile data into a str object, necessary since
                    #containerZip.read() returns a bytes object.
                    jsonString = io.TextIOWrapper(containerZip.open('data'),
                                                  encoding='utf-8'
                                                  ).read()
                self.sic[specfile] = json.loads(jsonString,
                                                object_hook=Si.jsonHook
                                                )
                msrunInfo['status']['si'] = True


class Ci(object):
    """Chromatogram item (Ci), representation of a mzML ``chromatogram``.

    :ivar id: The unique id of this chromatogram. Typically descriptive for the
        chromatogram, eg "TIC" (total ion current). Is used together with
        ``self.specfile`` as a key to access the spectrum in its container
        :class:`MsrunContainer.cic <maspy.core.MsrunContainer>`.
    :ivar specfile: An id representing a group of spectra, typically of the same
        mzML file / ms-run.
    :ivar id:
    :ivar dataProcessingRef: This attribute can optionally reference the 'id' of
        the appropriate dataProcessing, from mzML.
    :ivar precursor: The method of precursor ion selection and activation, from
        mzML.
    :ivar product: The method of product ion selection and activation in a
        precursor ion scan, from mzML.
    :ivar params: A list of parameter tuple, #TODO: as described elsewhere
    :ivar arrays: a dictionary containing the binary data of a chromatogram as
        ``numpy.array``. Keys are derived from the specified mzML cvParam, see
        :func:`maspy.xml.findBinaryDataType()`. Typically contains at least a
        time parameter ``rt`` (retention time) ``Ci.arrays = {'rt':
        numpy.array(), ...}``
    :ivar arrayInfo: dictionary describing each data type present in
        ``.arrays``. ::

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
    __slots__ = ['id', 'specfile', 'dataProcessingRef', 'precursor', 'product',
                 'params', 'attrib', 'arrays', 'arrayInfo']

    def __init__(self, identifier, specfile):
        self.id = identifier
        self.specfile = specfile
        self.dataProcessingRef = None
        self.precursor = None
        self.product = None
        self.params = list()
        self.attrib = dict()
        self.arrays = dict()
        self.arrayInfo = dict()

    def __repr__(self):
        return 'maspy.core.Ci(id=%r, specfile=%r)' % (self.id, self.specfile)

    def _reprJSON(self):
        """Returns a JSON serializable represenation of a ``Ci`` class instance.
        Use :func:`maspy.core.Ci._fromJSON()` to generate a new ``Ci`` instance
        from the return value.

        :returns: a JSON serializable python object
        """
        return {'__Ci__': (self.id, self.specfile, self.dataProcessingRef,
                           self.precursor, self.product, self.params,
                           self.attrib, self.arrayInfo
                           )
                }

    @classmethod
    def _fromJSON(cls, jsonobject):
        """Generates a new instance of :class:`maspy.core.Ci` from a decoded
        JSON object (as generated by :func:`maspy.core.Ci._reprJSON()`).

        :param jsonobject: decoded JSON object

        :returns: a new instance of :class:`Ci`
        """
        newInstance = cls(jsonobject[0], jsonobject[1])
        attribDict = {}
        attribDict['dataProcessingRef'] = jsonobject[2]
        attribDict['precursor'] = jsonobject[3]
        attribDict['product'] = jsonobject[4]
        attribDict['params'] = [tuple(param) for param in jsonobject[5]]
        attribDict['attrib'] = jsonobject[6]
        attribDict['arrayInfo'] = dict()
        for arrayType, jsonEntry in viewitems(jsonobject[7]):
            arrayEntry = {'dataProcessingRef': jsonEntry['dataProcessingRef'],
                          'params': [tuple(_) for _ in jsonEntry['params']]
                          }
            attribDict['arrayInfo'][arrayType] = arrayEntry
        for key, value in viewitems(attribDict):
            setattr(newInstance, key, value)
        return newInstance

    @staticmethod
    def jsonHook(encoded):
        """Custom JSON decoder that allows construction of a new ``Ci`` instance
        from a decoded JSON object.

        :param encoded: a JSON decoded object literal (a dict)

        :returns: "encoded" or one of the these objects: :class:`Ci`,
            :class:`MzmlProduct`, :class:`MzmlPrecursor`
        """
        if '__Ci__' in encoded:
            return Ci._fromJSON(encoded['__Ci__'])
        elif '__MzmlProduct__' in encoded:
            return MzmlProduct._fromJSON(encoded['__MzmlProduct__'])
        elif '__MzmlPrecursor__' in encoded:
            return MzmlPrecursor._fromJSON(encoded['__MzmlPrecursor__'])
        else:
            return encoded


class Sai(object):
    """Spectrum array item (Sai), representation of the binary data arrays of an
    mzML ``spectrum``.

    :ivar id: The unique id of this spectrum, typically the scan number. Is used
        together with ``self.specfile`` as a key to access the spectrum in its
        container :class:`MsrunContainer.saic <maspy.core.MsrunContainer>`.
        Should be derived from the spectrums nativeID format (MS:1000767).
    :ivar specfile: An id representing a group of spectra, typically of the same
        mzML file / ms-run.
    :ivar arrays: a dictionary containing the binary data of the recorded ion
        spectrum as ``numpy.array``. Keys are derived from the specified mzML
        cvParam, see :func:`maspy.xml.findBinaryDataType()`. Must at least
        contain the keys ``mz`` (mass to charge ratio) and ``i`` (intensity).
        ``Sai.arrays = {'mz': numpy.array(), 'i': numpy.array(), ...}``
    :ivar arrayInfo: dictionary describing each data type present in
        ``.arrays``. ::

            {dataType: {'dataProcessingRef': str,
                        'params': [paramTuple, paramTuple, ...]
                        }
             }

        code example::

            {u'i': {u'dataProcessingRef': None,
                    u'params': [('MS:1000521', '', None),
                                ('MS:1000574', '', None),
                                ('MS:1000515', '', 'MS:1000131')]},
             u'mz': {u'dataProcessingRef': None,
                     u'params': [('MS:1000523', '', None),
                                 ('MS:1000574', '', None),
                                 ('MS:1000514', '', 'MS:1000040')]}}
    """
    __slots__ = ['id', 'specfile', 'arrays', 'arrayInfo']

    def __init__(self, identifier, specfile):
        self.id = identifier
        self.specfile = specfile

        self.arrays = dict()
        self.arrayInfo = dict()

    def __repr__(self):
        return 'maspy.core.Sai(id=%r, specfile=%r)' % (self.id, self.specfile)

    def _reprJSON(self):
        """Returns a JSON serializable represenation of a ``Sai`` class
        instance. Use :func:`maspy.core.Sai._fromJSON()` to generate a new
        ``Sai`` instance from the return value.

        :returns: a JSON serializable python object
        """
        return {'__Sai__': (self.id, self.specfile, self.arrayInfo)}

    @classmethod
    def _fromJSON(cls, jsonobject):
        """Generates a new instance of :class:`maspy.core.Sai` from a decoded
        JSON object (as generated by :func:`maspy.core.Sai._reprJSON()`).

        :param jsonobject: decoded JSON object

        :returns: a new instance of :class:`Sai`
        """
        newInstance = cls(jsonobject[0], jsonobject[1])
        for arrayType, jsonEntry in viewitems(jsonobject[2]):
            arrayEntry = {'dataProcessingRef': jsonEntry['dataProcessingRef'],
                          'params': [tuple(_) for _ in jsonEntry['params']]
                          }
            newInstance.arrayInfo[arrayType] = arrayEntry

        return newInstance

    @staticmethod
    def jsonHook(encoded):
        """Custom JSON decoder that allows construction of a new ``Sai``
        instance from a decoded JSON object.

        :param encoded: a JSON decoded object literal (a dict)

        :returns: "encoded" or :class:`Sai`
        """
        if '__Sai__' in encoded:
            return Sai._fromJSON(encoded['__Sai__'])
        else:
            return encoded


class Smi(object):
    """Spectrum metadata item (Smi), representation of all the metadata data of
    an mzML ``spectrum``, excluding the actual binary data.

    For details on the mzML ``spectrum`` element refer to the `documentation,
    <http://www.peptideatlas.org/tmp/mzML1.1.0.html#spectrum>`_.

    :ivar id: The unique id of this spectrum, typically the scan number. Is used
        together with ``self.specfile`` as a key to access the spectrum in its
        container :class:`MsrunContainer.smic <maspy.core.MsrunContainer>`.
        Should be derived from the spectrums nativeID format (MS:1000767).
    :ivar specfile: An id representing a group of spectra, typically of the same
        mzML file / ms-run.
    :ivar attributes: dict, attributes of an mzML ``spectrum`` element
    :ivar params: a list of parameter tuple (cvParam tuple, userParam tuple or
        referencableParamGroup tuple) of an mzML ``spectrum`` element.
    :ivar scanListParams: a list of parameter tuple (cvParam tuple, userParam
        tuple or referencableParamGroup tuple) of an mzML ``scanList`` element.
    :ivar scanList: a list of :class:`MzmlScan` elements, derived from elements
        of an an mzML ``scanList`` element.
    :ivar precursorList: a list of :class:`MzmlPrecursor` elements, derived from
        elements of an an mzML ``precursorList`` element.
    :ivar productList: a list of :class:`MzmlProduct` elements, derived from
        elements of an an mzML ``productList`` element.

    .. warning::
        The ``Smi`` is used to generate ``spectrum`` xml elements by using the
        function :func:`maspy.writer.xmlSpectrumFromSmi()`. In order to generate
        a valid mzML element all attributes of ``Smi`` have to be in the correct
        format. Therefore it is highly recommended to only use properly
        implemented and tested methods for making changes to any ``Smi``
        attribute.
    """
    __slots__ = ['id', 'specfile', 'attributes', 'params', 'scanListParams',
                 'scanList', 'precursorList', 'productList'
                 ]

    def __init__(self, identifier, specfile):
        self.id = identifier
        self.specfile = specfile
        self.attributes = dict()
        self.params = list()
        self.scanListParams = list()
        self.scanList = list()
        self.precursorList = list()
        self.productList = list()

    def __repr__(self):
        return 'maspy.core.Smi(id=%r, specfile=%r)' % (self.id, self.specfile)

    def _reprJSON(self):
        """Returns a JSON serializable represenation of a ``Smi`` class
        instance. Use :func:`maspy.core.Sai._fromJSON()` to generate a new
        ``Smi`` instance from the return value.

        :returns: a JSON serializable python object
        """
        return {'__Smi__': (self.id, self.specfile, self.attributes,
                            self.params, self.scanListParams, self.scanList,
                            self.precursorList, self.productList
                            )
                }

    @classmethod
    def _fromJSON(cls, jsonobject):
        """Generates a new instance of :class:`maspy.core.Smi` from a decoded
        JSON object (as generated by :func:`maspy.core.Smi._reprJSON()`).

        :param jsonobject: decoded JSON object

        :returns: a new instance of :class:`Smi`
        """
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
        """Custom JSON decoder that allows construction of a new ``Smi``
        instance from a decoded JSON object.

        :param encoded: a JSON decoded object literal (a dict)

        :returns: "encoded" or one of the these objects: :class:`Smi`,
            :class:`MzmlScan`, :class:`MzmlProduct`, :class:`MzmlPrecursor`
        """
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
    """Spectrum item (Si) - this is the spectrum representation intended to be
    used in maspy. A simplified representation of spectrum metadata. Contains
    only specifically imported attributes, which are necessary for data
    analysis. Does not follow any PSI data structure or name space rules.

    Additional attributes can be transferred from the corresponding :class:`Smi`
    entry. This is done by default when importing an mzML file by using the
    function :func:`maspy.reader.defaultFetchSiAttrFromSmi()`.

    :ivar id: The unique id of this spectrum, typically the scan number. Is used
        together with ``self.specfile`` as a key to access the spectrum in its
        container :class:`MsrunContainer.sic <maspy.core.MsrunContainer>`.
        Should be derived from the spectrums nativeID format (MS:1000767).
    :ivar specfile: An id representing a group of spectra, typically of the same
        mzML file / ms-run.
    :ivar isValid: bool, can be used for filtering.
    :ivar msLevel: stage of ms level in a multi stage mass spectrometry
        experiment.
    """
    def __init__(self, identifier, specfile):
        self.id = identifier
        self.specfile = specfile
        self.isValid = None
        self.msLevel = None

    def __repr__(self):
        return 'maspy.core.Si(id=%r, specfile=%r)' % (self.id, self.specfile)

    def __str__(self):
        output = 'maspy.core.Si(id=%s, specfile=%s)' % (self.id, self.specfile)
        output += '\n  ' + 'keys=%r' % (sorted(viewkeys(self.__dict__)))
        return output

    def _reprJSON(self):
        """Returns a JSON serializable represenation of a ``Si`` class instance.
        Use :func:`maspy.core.Si._fromJSON()` to generate a new ``Si`` instance
        from the return value.

        :returns: a JSON serializable python object
        """
        return {'__Si__': self.__dict__}

    @classmethod
    def _fromJSON(cls, jsonobject):
        """Generates a new instance of :class:`maspy.core.Si` from a decoded
        JSON object (as generated by :func:`maspy.core.Si._reprJSON()`).

        :param jsonobject: decoded JSON object

        :returns: a new instance of :class:`Si`
        """
        newInstance = cls(None, None)
        newInstance.__dict__.update(jsonobject)
        return newInstance

    @staticmethod
    def jsonHook(encoded):
        """Custom JSON decoder that allows construction of a new ``Si``
        instance from a decoded JSON object.

        :param encoded: a JSON decoded object literal (a dict)

        :returns: "encoded" or :class:`Si`
        """
        if '__Si__' in encoded:
            return Si._fromJSON(encoded['__Si__'])
        else:
            return encoded


class MzmlScan(object):
    """MasPy representation of an mzML ``Scan`` element, see `mzML schema
    documentation <http://www.peptideatlas.org/tmp/mzML1.1.0.html#scan>`_.

    :ivar scanWindowList: a list of mzML ``scanWindow`` elements, which are
        represented as a tuple of parm tuples. The mzML ``scanWindowList`` is
        describing the measurement and should not be changed.
    :ivar params: a list of parameter tuple (cvParam tuple, userParam tuple or
        referencableParamGroup tuple) of an mzML ``Scan`` element.

    .. note:: The attributes "sourceFileRef" and "externalSpectrumID" are not
        supported by MasPy on purpose, since they are only used to refere to
        scans which are external to the mzML file. The attribute "spectrumRef"
        could be included but seems kind of useless.

        The attribute "instrumentConfigurationRef" should be included though:
        #TODO.
    """
    #kwargs are used to only take the arguments needed and ignore additionally
    #specified ones like 'arrayLength' of binaryDataArray
    __slots__ = ['scanWindowList', 'params']
    def __init__(self, scanWindowList=None, params=None, **kwargs):
        if scanWindowList is not None:
            self.scanWindowList = tuple(tuple(_) for _ in scanWindowList)
        else:
            self.scanWindowList = None
        self.params = params

    def __getitem__(self, key):
        return getattr(self, key)

    #def __setitem__(self, key, value):
    #    setattr(self, key, value)

    def _reprJSON(self):
        """Returns a JSON serializable represenation of a ``MzmlScan`` class
        instance. Use :func:`maspy.core.MzmlScan._fromJSON()` to generate a new
        ``MzmlScan`` instance from the return value.

        :returns: a JSON serializable python object
        """
        return {'__MzmlScan__': (self.scanWindowList, self.params)}

    @classmethod
    def _fromJSON(cls, jsonobject):
        """Generates a new instance of :class:`maspy.core.MzmlScan` from a
        decoded JSON object (as generated by
        :func:`maspy.core.MzmlScan._reprJSON()`).

        :param jsonobject: decoded JSON object

        :returns: a new instance of :class:`MzmlScan`
        """
        scanWindowList = _mzmlListAttribToTuple(jsonobject[0])
        params = [tuple(param) for param in jsonobject[1]]
        return cls(scanWindowList, params)


class MzmlProduct(object):
    """MasPy representation of an mzML ``Product`` element, the `mzML schema
    documentation <http://www.peptideatlas.org/tmp/mzML1.1.0.html#product>`_
    does however not provide a lot of information how this element is intended
    to be used and which information can be present.

    :ivar isolationWindow: the mzML ``isolationWindow`` is represented as a
        tuple of parm tuples. It is describing the measurement and should not be
        changed.
    """
    __slots__ = ['isolationWindow']
    def __init__(self, isolationWindow=None, **kwargs):
        if isolationWindow is not None:
            self.isolationWindow = tuple(isolationWindow)
        else:
            self.isolationWindow = None

    def __getitem__(self, key):
        return getattr(self, key)

    #def __setitem__(self, key, value):
    #    setattr(self, key, value)

    def _reprJSON(self):
        """Returns a JSON serializable represenation of a ``MzmlProduct`` class
        instance. Use :func:`maspy.core.MzmlProduct._fromJSON()` to generate a
        new ``MzmlProduct`` instance from the return value.

        :returns: a JSON serializable python object
        """
        return {'__MzmlProduct__': self.isolationWindow}

    @classmethod
    def _fromJSON(cls, jsonobject):
        """Generates a new instance of :class:`maspy.core.MzmlProduct` from a
        decoded JSON object (as generated by
        :func:`maspy.core.MzmlProduct._reprJSON()`).

        :param jsonobject: decoded JSON object

        :returns: a new instance of :class:`MzmlProduct`
        """
        isolationWindow =[tuple(param) for param in jsonobject]
        return cls(isolationWindow)


class MzmlPrecursor(object):
    """MasPy representation of an mzML ``Scan`` element, see `mzML schema
    documentation <http://www.peptideatlas.org/tmp/mzML1.1.0.html#precursor>`_.

    :ivar spectrumRef: native id of the spectrum corresponding to the precursor
        spectrum
    :ivar activation: the mzML ``activation`` is represented as a tuple of param
        tuples. It is describing the type and energy level used for activation
        and should not be changed.
    :ivar isolationWindow: the mzML ``isolationWindow`` is represented as a
        tuple of parm tuples. It is describing the measurement and should not be
        changed.
    :ivar selectedIonList: a list of mzML ``selectedIon`` elements, which are
        represented as a tuple of param tuples.

    .. note:: The attributes "sourceFileRef" and "externalSpectrumID" are not
        supported by MasPy on purpose, since they are only used to refere to
        scans which are external to the mzML file.
    """
    __slots__ = ['spectrumRef', 'activation', 'isolationWindow',
                 'selectedIonList'
                 ]
    def __init__(self, spectrumRef=None, activation=None, isolationWindow=None,
                 selectedIonList=None, **kwargs):
        self.spectrumRef = spectrumRef

        if activation is not None:
            self.activation = tuple(activation)
        else:
            self.activation = None

        if isolationWindow is not None:
            self.isolationWindow = tuple(isolationWindow)
        else:
            self.isolationWindow = None

        if selectedIonList is not None:
            self.selectedIonList = [tuple(_) for _ in selectedIonList]
        else:
            self.selectedIonList = None

    def __getitem__(self, key):
        return getattr(self, key)

    #def __setitem__(self, key, value):
    #    setattr(self, key, value)

    def _reprJSON(self):
        """Returns a JSON serializable represenation of a ``MzmlPrecursor``
        class instance. Use :func:`maspy.core.MzmlPrecursor._fromJSON()` to
        generate a new ``MzmlPrecursor`` instance from the return value.

        :returns: a JSON serializable python object
        """
        return {'__MzmlPrecursor__': (self.spectrumRef, self.activation,
                                      self.isolationWindow, self.selectedIonList
                                      )
                }

    @classmethod
    def _fromJSON(cls, jsonobject):
        """Generates a new instance of :class:`maspy.core.MzmlPrecursor` from a
        decoded JSON object (as generated by
        :func:`maspy.core.MzmlPrecursor._reprJSON()`).

        :param jsonobject: decoded JSON object

        :returns: a new instance of :class:`MzmlPrecursor`
        """
        spectrumRef = jsonobject[0]
        activation = [tuple(param) for param in jsonobject[1]]
        isolationWindow =[tuple(param) for param in jsonobject[2]]
        selectedIonList = _mzmlListAttribToTuple(jsonobject[3])
        return cls(spectrumRef, activation, isolationWindow, selectedIonList)


def _mzmlListAttribToTuple(oldList):
    """Turns the param entries of elements in a list elements into tuples, used
    in :func:`MzmlScan._fromJSON()` and :func:`MzmlPrecursor._fromJSON()`.

    .. note:: only intended for a list of elements that contain params. For
        example the mzML element ``selectedIonList`` or ``scanWindowList``.

    :param oldList: [[paramList, paramList, ...], ...]

    :returns: [[paramTuple, paramTuple, ...], ...]
    """
    newList = list()
    for oldParamList in oldList:
        newParamLIst = [tuple(param) for param in oldParamList]
        newList.append(newParamLIst)
    return newList


def addMsrunContainers(mainContainer, subContainer):
    """Adds the complete content of all specfile entries from the subContainer
    to the mainContainer. However if a specfile of ``subContainer.info`` is
    already present in ``mainContainer.info`` its contents are not added to the
    mainContainer.

    :param mainContainer: :class:`MsrunContainer`
    :param subContainer: :class:`MsrunContainer`

    .. warning:: does not generate new items, all items added to the
        ``mainContainer`` are still present in the ``subContainer`` and changes
        made to elements of one container also affects the elements of the other
        one (ie elements share same memory location).
    """
    typeToContainer = {'rm': 'rmc', 'ci': 'cic', 'smi': 'smic',
                       'sai': 'saic', 'si': 'sic'
                       }
    for specfile in subContainer.info:
        if specfile in mainContainer.info:
            continue

        mainContainer.addSpecfile(specfile, subContainer.info[specfile]['path'])
        for datatype, status in listitems(subContainer.info[specfile]['status']):
            if not status:
                continue
            datatypeContainer = typeToContainer[datatype]
            dataTypeContainer = getattr(mainContainer, datatypeContainer)
            subContainerData = getattr(subContainer,
                                       datatypeContainer
                                       )[specfile]
            dataTypeContainer[specfile] = subContainerData
            mainContainer.info[specfile]['status'][datatype] = True


##########################################################################
### SpectrumIdentificationItem related classes and functions #############
##########################################################################
class Sii(object):
    """Spectrum identification item (Sii) - representation of an MSn fragment
    spectrum annotation, also referred to as peptide spectrum match (PSM).

    :ivar id: The unique id of this spectrum, typically the scan number. Is used
        together with ``self.specfile`` as a key to access the spectrum in its
        container :class:`SiiContainer` or the corresponding spectrum in a
        :class:`MsrunContainer`.
    :ivar specfile: An id representing an mzML file / ms-run filename.
    :ivar rank: The rank of this ``Sii`` compared to others for the same MSn
        spectrum. The rank is based on a score defined in the ``SiiContainer``.
        If multiple Sii have the same top score, they should all be assigned
        ``self.rank = 1``.
    :ivar isValid: bool or None if not specified
        this attribute can be used to flag if a Sii has passed a given quality
        threshold or has been validated as correct. Is used to filter valid
        elements in eg :func:`SiiContainer.getArrays()`.
    """
    def __init__(self, identifier, specfile):
        self.id = identifier
        self.specfile = specfile
        self.rank = None
        self.isValid = None

    def _reprJSON(self):
        """Returns a JSON serializable represenation of a ``Sii`` class
        instance. Use :func:`maspy.core.Sii._fromJSON()` to generate a new
        ``Sii`` instance from the return value.

        :returns: a JSON serializable python object
        """
        return {'__Sii__': self.__dict__}

    @classmethod
    def _fromJSON(cls, jsonobject):
        """Generates a new instance of :class:`maspy.core.Sii` from a decoded
        JSON object (as generated by :func:`maspy.core.Sii._reprJSON()`).

        :param jsonobject: decoded JSON object

        :returns: a new instance of :class:`Sii`
        """
        newInstance = cls(None, None)
        newInstance.__dict__.update(jsonobject)
        return newInstance

    @staticmethod
    def jsonHook(encoded):
        """Custom JSON decoder that allows construction of a new ``Sii``
        instance from a decoded JSON object.

        :param encoded: a JSON decoded object literal (a dict)

        :returns: "encoded" or :class:`Sii`
        """
        if '__Sii__' in encoded:
            return Sii._fromJSON(encoded['__Sii__'])
        else:
            return encoded


class SiiContainer(object):
    """Conainer for :class:`Sii` elements.

    :ivar container: contains the stored :class:`Sii <maspy.core.Sii>`
        elements. ::

            {specfilename: {'Sii identifier': [Sii, ...], ...}

    :ivar info: a dictionary containing information about imported specfiles. ::

            {specfilename: {'path': str, 'qcAttr': str, 'qcLargerBetter': bool,
                            'qcCutoff': float, 'rankAttr': str,
                            'rankLargerBetter': bool
                            },
             ...
             }

        **path**: folder location used by the ``SiiContainer`` to save and load
        data to the hard disk.

        **qcAttr**: name of the parameter to define a quality cutoff. Typically
        this is some sort of a global false positive estimator (eg FDR)

        **qcLargerBetter**: bool, True if a large value for the ``.qcAttr``
        means a higher confidence.

        **qcCutoff**: float, the quality threshold for the specifed ``.qcAttr``

        **rankAttr**: name of the parameter used for ranking ``Sii`` according
        to how well they match to a fragment ion spectrum, in the case when
        their are multiple ``Sii`` present for the same spectrum.

        **rankLargerBetter**: bool, True if a large value for the ``.rankAttr``
        means a better match to the fragment ion spectrum

    .. note:: In the future this container may be integrated in an evidence or
        an mzIdentML like container.
    """
    def __init__(self):
        self.container = dict()
        self.info = dict()

    def getArrays(self, attr=None, specfiles=None, sort=False, reverse=False,
                  selector=None, defaultValue=None):
        """Return a condensed array of data selected from :class:`Sii` instances
        from ``self.container`` for fast and convenient data processing.

        :param attr: list of :class:`Sii` item attributes that should be added
            to the returned array. The attributes "id" and "specfile" are always
            included, in combination they serve as a unique id.
        :param defaultValue: if an item is missing an attribute, the
            "defaultValue" is added to the array instead.
        :param specfiles: filenames of ms-run files - if specified return only
            items from those files
        :type specfiles: str or [str, str, ...]
        :param sort: if "sort" is specified the returned list of items is sorted
            according to the :class:`Sii` attribute specified by "sort", if the
            attribute is not present the item is skipped.
        :param reverse: bool, set True to reverse sort order
        :param selector: a function which is called with each `Sii` item and has
            to return True (include item) or False (discard item).
            Default function is: ``lambda si: True``. By default only items with
            ``Sii.isValid == True`` are returned.

        :returns: {'attribute1': numpy.array(),
                   'attribute2': numpy.array(),
                   ...
                   }
        """
        selector = (lambda sii: sii.isValid) if selector is None else selector
        attr = attr if attr is not None else []
        attr = set(['id', 'specfile'] + aux.toList(attr))
        items = self.getItems(specfiles, sort, reverse, selector)
        return _getArrays(items, attr, defaultValue)

    def getItems(self, specfiles=None, sort=False, reverse=False,
                 selector=None):
        """Generator that yields filtered and/or sorted :class:`Sii` instances
        from ``self.container``.

        :param specfiles: filenames of ms-run files - if specified return only
            items from those files
        :type specfiles: str or [str, str, ...]
        :param sort: if "sort" is specified the returned list of items is sorted
            according to the :class:`Sii` attribute specified by "sort", if the
            attribute is not present the item is skipped.
        :param reverse: bool, ``True`` reverses the sort order
        :param selector: a function which is called with each ``Sii`` item and
            has to return True (include item) or False (discard item). By
            default only items with ``Sii.isValid == True`` are returned.

        :returns: items from container that passed the selector function
        """
        selector = (lambda sii: sii.isValid) if selector is None else selector
        if specfiles is None:
            specfiles = [_ for _ in viewkeys(self.info)]
        else:
            specfiles = aux.toList(specfiles)
        return _getListItems(self.container, specfiles, sort, reverse, selector)

    def getValidItem(self, specfile, identifier):
        """Returns a ``Sii`` instance from ``self.container`` if it is valid,
        if all elements of ``self.container[specfile][identifier] are
        ``Sii.isValid == False`` then ``None`` is returned.

        :param specfile: a ms-run file name
        :param identifier: item identifier ``Sii.id``

        :returns: ``Sii`` or ``None``
        """
        for item in self.container[specfile][identifier]:
            if item.isValid:
                return item
        else:
            return None

    def addSpecfile(self, specfiles, path):
        """Prepares the container for loading ``siic`` files by adding specfile
        entries to ``self.info``. Use :func:`SiiContainer.load()` afterwards to
        actually import the files.

        :param specfiles: the name of an ms-run file or a list of names
        :type specfiles: str or [str, str, ...]
        :param path: filedirectory used for loading and saving ``siic`` files
        """
        for specfile in aux.toList(specfiles):
            if specfile not in self.info:
                self._addSpecfile(specfile, path)
            else:
                warntext = 'Error while calling "SiiContainer.addSpecfile(): "'\
                           '"%s" is already present "SiiContainer.info"'\
                            % (specfile, )
                warnings.warn(warntext)

    def _addSpecfile(self, specfile, path):
        """Adds a new specfile entry to SiiContainer.info. See also
        :class:`SiiContainer.addSpecfile()`.

        :param specfile: the name of an ms-run file
        :param path: filedirectory for loading and saving the ``siic`` files
        """
        self.info[specfile] = {'path': path, 'qcAttr': None, 'qcCutoff': None,
                               'qcLargerBetter': None, 'rankAttr': None,
                               'rankLargerBetter': None
                               }
        self.container[specfile] = dict()

    def setPath(self, folderpath, specfiles=None):
        """Changes the folderpath of the specified specfiles. The folderpath is
        used for saving and loading of ``siic`` files.

        :param folderpath: a filedirectory
        :param specfiles: the name of an ms-run file or a list of names. If None
            all specfiles are selected.
        :type specfiles: None, str, [str, str]
        """
        if specfiles is None:
            specfiles = [_ for _ in viewkeys(self.info)]
        else:
            specfiles = aux.toList(specfiles)

        _containerSetPath(self, folderpath, specfiles)


    def removeSpecfile(self, specfiles):
        """Completely removes the specified specfiles from the ``SiiContainer``.

        :param specfiles: the name of an ms-run file or a list of names.
        """
        for specfile in aux.toList(specfiles):
            del self.container[specfile]
            del self.info[specfile]

    def save(self, specfiles=None, compress=True, path=None):
        """Writes the specified specfiles to ``siic`` files on the hard disk.

        .. note::
            If ``.save()`` is called and no ``siic`` files are present in the
            specified path new files are generated, otherwise old files are
            replaced.

        :param specfiles: the name of an ms-run file or a list of names. If None
            all specfiles are selected.
        :param compress: bool, True to use zip file compression
        :param path: filedirectory to which the ``siic`` files are written. By
            default the parameter is set to ``None`` and the filedirectory is
            read from ``self.info[specfile]['path']``
        """
        if specfiles is None:
            specfiles = [_ for _ in viewkeys(self.info)]
        else:
            specfiles = aux.toList(specfiles)

        for specfile in specfiles:
            if specfile not in self.info:
                warntext = 'Error while calling "SiiContainer.save()": "%s" is'\
                           ' not present in "SiiContainer.info"!'\
                            % (specfile, )
                warnings.warn(warntext)
                continue
            else:
                path = self.info[specfile]['path'] if path is None else path

            with aux.PartiallySafeReplace() as msr:
                filename = specfile + '.siic'
                filepath = aux.joinpath(path, filename)
                with msr.open(filepath, mode='w+b') as openfile:
                    self._writeContainer(openfile, specfile, compress)

    def _writeContainer(self, filelike, specfile, compress):
        """Writes the ``self.container`` entry of the specified specfile to the
        ``siic`` format. In addition it also dumps the ``self.info`` entry to
        the zipfile with the filename ``info``. For details see
        :func:`maspy.auxiliary.writeJsonZipfile()`

        :param filelike: path to a file (str) or a file-like object
        :param specfile: name of an ms-run file present in ``self.info``
        :param compress: bool, True to use zip file compression
        """
        aux.writeJsonZipfile(filelike, self.container[specfile],
                             compress=compress
                             )
        zipcomp = zipfile.ZIP_DEFLATED if compress else zipfile.ZIP_STORED
        with zipfile.ZipFile(filelike, 'a', allowZip64=True) as containerFile:
            infodata = {key: value for key, value in
                        viewitems(self.info[specfile]) if key != 'path'
                        }
            containerFile.writestr('info', json.dumps(infodata, zipcomp))

    def load(self, specfiles=None):
        """Imports ``siic`` files from the hard disk.

        :param specfiles: the name of an ms-run file or a list of names. If None
            all specfiles are selected.
        :type specfiles: None, str, [str, str]
        """
        if specfiles is None:
            specfiles = [_ for _ in viewkeys(self.info)]
        else:
            specfiles = aux.toList(specfiles)

        for specfile in specfiles:
            if specfile not in self.info:
                warntext = 'Error while calling "SiiContainer.load()": "%s" is'\
                           ' not present in "SiiContainer.info"!'\
                            % (specfile, )
                warnings.warn(warntext)
                continue
            else:
                siiPath = aux.joinpath(self.info[specfile]['path'],
                                       specfile+'.siic'
                                       )

            with zipfile.ZipFile(siiPath, 'r') as containerZip:
                #Convert the zipfile data into a str object, necessary since
                #containerZip.read() returns a bytes object.
                jsonString = io.TextIOWrapper(containerZip.open('data'),
                                              encoding='utf-8'
                                              ).read()
                infoString = io.TextIOWrapper(containerZip.open('info'),
                                              encoding='utf-8'
                                              ).read()
            self.container[specfile] = json.loads(jsonString,
                                                  object_hook=Sii.jsonHook
                                                  )
            self.info[specfile].update(json.loads(infoString))

    def addSiInfo(self, msrunContainer, specfiles=None,
                  attributes=['obsMz', 'rt', 'charge']):
        """Transfer attributes to :class:`Sii` elements from the corresponding
        :class`Si` in :class:`MsrunContainer.sic <MsrunContainer>`. If an
        attribute is not present in the ``Si`` the attribute value in the
        ``Sii``is set to ``None``.

        Attribute examples: 'obsMz', 'rt', 'charge', 'tic', 'iit', 'ms1Id'

        :param msrunContainer: an instance of :class:`MsrunContainer` which has
            imported the corresponding specfiles
        :param specfiles: the name of an ms-run file or a list of names. If None
            all specfiles are selected.
        :param attributes: a list of ``Si`` attributes that should be
            transfered.
        """
        if specfiles is None:
            specfiles = [_ for _ in viewkeys(self.info)]
        else:
            specfiles = aux.toList(specfiles)

        for specfile in specfiles:
            if specfile not in self.info:
                warntext = 'Error while calling "SiiContainer.addSiInfo()": '\
                           '"%s" is not present in "SiiContainer.info"!'\
                            % (specfile, )
                warnings.warn(warntext)
            elif specfile not in msrunContainer.info:
                warntext = 'Error while calling "SiiContainer.addSiInfo()": '\
                           '"%s" is not present in "MsrunContainer.info"'\
                            % (specfile, )
                warnings.warn(warntext)
            else:
                for identifier in self.container[specfile]:
                    si = msrunContainer.sic[specfile][identifier]
                    for sii in self.container[specfile][identifier]:
                        for attribute in attributes:
                            setattr(sii, attribute,
                                    getattr(si, attribute, None)
                                    )

    def calcMz(self, specfiles=None, guessCharge=True, obsMzKey='obsMz'):
        """Calculate the exact mass for ``Sii`` elements from the
        ``Sii.peptide`` sequence.

        :param specfiles: the name of an ms-run file or a list of names. If None
            all specfiles are selected.
        :param guessCharge: bool, True if the charge should be guessed if the
            attribute ``charge`` is missing from ``Sii``. Uses the calculated
            peptide mass and the observed m/z value to calculate the charge.
        :param obsMzKey: attribute name of the observed m/z value in ``Sii``.
        """
        #TODO: important to test function, since changes were made
        _calcMass = maspy.peptidemethods.calcPeptideMass
        _calcMzFromMass = maspy.peptidemethods.calcMzFromMass
        _massProton = maspy.constants.atomicMassProton
        _guessCharge = lambda mass, mz: round(mass / (mz - _massProton), 0)

        if specfiles is None:
            specfiles = [_ for _ in viewkeys(self.info)]
        else:
            specfiles = aux.toList(specfiles)

        tempMasses = dict()
        for specfile in specfiles:
            if specfile not in self.info:
                warntext = 'Error while calling "SiiContainer.calcMz()": '\
                           '"%s" is not present in "SiiContainer.info"!'\
                            % (specfile, )
                warnings.warn(warntext)
            else:
                for sii in self.getItems(specfiles=specfile):
                    peptide = sii.peptide
                    if peptide not in tempMasses:
                        if hasattr(sii, 'diPeptide'):
                            tempMasses[peptide] = (_calcMass(sii.peptide1) +
                                                   _calcMass(sii.peptide2)
                                                   )
                        else:
                            tempMasses[peptide] = _calcMass(peptide)
                    peptideMass = tempMasses[peptide]
                    if sii.charge is not None:
                        sii.excMz = _calcMzFromMass(peptideMass, sii.charge)
                    elif guessCharge:
                        guessedCharge = _guessCharge(peptideMass,
                                                     getattr(sii, obsMzKey)
                                                     )
                        sii.excMz = _calcMzFromMass(peptideMass, guessedCharge)
                        sii.charge = guessedCharge
                    else:
                        sii.excMz = None
        del(tempMasses)


###########################################################
### FeatureItem related classes and functions #############
###########################################################
class Fi(object):
    """Feature item (Fi), representation of a peptide LC-MS feature.

    :ivar id: the unique identifier of a LC-MS feature, as generated by the
        software used for extracting features from MS1 spectra.
    :ivar specfile: An id representing an mzML file / ms-run filename.
    :ivar rt: a representative retention time value of the ``Fi`` (in seconds).
        For example the retention time of the feature apex.
    :ivar mz: a representative mass to charge value of the ``Fi`` (in Dalton /
        charge). For example the average m/z value of all data points.
    :ivar charge: the ``Fi`` charge state
    :ivar intensity: a meassure for the ``Fi`` abundance, used for relative
        quantification. Typically the area of the feature intensity over time.
    :ivar isValid:  bool or None if not specified
        this attribute can be used to flag if a ``Fi`` has passed a given
        quality threshold. Can be used to filter valid elements in eg
        :func:`FiContainer.getArrays()`.
    :ivar isMatched: bool or None if not specified
        True if any ``Si`` or ``Sii`` elements could be matched. Should be set
        to False on import.
    :ivar isAnnotated: bool or None if not specified
        True if any ``Sii`` elements could be matched. Should be set to False on
        import. Not sure yet how to handle annotation from other features.
    :ivar siIds: list of tuple(specfile, id) from matched :class:`Si`
    :ivar siiIds: list of tuple(specfile, id) from matched :class:`Sii`
    :ivar peptide: peptide sequence containing amino acid modifications. If
        multiple peptide sequences are possible due to multiple ``Sii`` matches
        the most likely must be chosen. A simple and accepted way to do this is
        by choosing the ``Sii`` identification with the best score.
    :ivar sequence: the plain amino acid sequence of ``self.peptide``
    :ivar bestScore: the score of the acceppted ``Sii`` for annotation
    """
    def __init__(self, identifier, specfile):
        self.id = identifier
        self.specfile = specfile
        self.isValid = None
        self.rt = float()
        self.mz = float()
        self.charge = int()

        #Annotation information
        self.isMatched = None
        self.isAnnotated = None
        self.siIds = list()
        self.siiIds = list()
        self.peptide = None
        self.sequence = None

    def _reprJSON(self):
        """Returns a JSON serializable represenation of a ``Fi`` class instance.
        Use :func:`maspy.core.Fi._fromJSON()` to generate a new ``Fi`` instance
        from the return value.

        :returns: a JSON serializable python object
        """
        return {'__Fi__': self.__dict__}

    @classmethod
    def _fromJSON(cls, jsonobject):
        """Generates a new instance of :class:`maspy.core.Fi` from a decoded
        JSON object (as generated by :func:`maspy.core.Fi._reprJSON()`).

        :param jsonobject: decoded JSON object

        :returns: a new instance of :class:`Fi`
        """
        newInstance = cls(None, None)
        newInstance.__dict__.update(jsonobject)
        return newInstance

    @staticmethod
    def jsonHook(encoded):
        """Custom JSON decoder that allows construction of a new ``Fi``
        instance from a decoded JSON object.

        :param encoded: a JSON decoded object literal (a dict)

        :returns: "encoded" or :class:`Fi`
        """
        if '__Fi__' in encoded:
            return Fi._fromJSON(encoded['__Fi__'])
        else:
            return encoded


class FiContainer(object):
    """Conainer for :class:`Fi` elements.

    :ivar container: contains the stored :class:`Fi <maspy.core.Fi>`
        elements. ::

            {specfilename: {'Fi identifier': [Fi, ...], ...}

    :ivar info: a dictionary containing information about imported specfiles. ::

            {specfilename: {'path': str},
             ...
             }

        **path**: folder location used by the ``FiContainer`` to save and load
        data to the hard disk.

    """
    def __init__(self):
        self.container = dict()
        self.info = dict()

    def getItem(self, specfile, identifier):
        """Returns a :class:`Fi` instance from ``self.container``.

        :param specfile: a ms-run file name
        :param identifier: item identifier ``Fi.id``

        :returns: ``self.container[specfile][identifier]``
        """
        return self.container[specfile][identifier]

    def getArrays(self, attr=None, specfiles=None, sort=False, reverse=False,
                  selector=None, defaultValue=None):
        """Return a condensed array of data selected from :class:`Fi` instances
        from ``self.container`` for fast and convenient data processing.

        :param attr: list of :class:`Fi` item attributes that should be added
            to the returned array. The attributes "id" and "specfile" are always
            included, in combination they serve as a unique id.
        :param defaultValue: if an item is missing an attribute, the
            "defaultValue" is added to the array instead.
        :param specfiles: filenames of ms-run files - if specified return only
            items from those files
        :type specfiles: str or [str, str, ...]
        :param sort: if "sort" is specified the returned list of items is sorted
            according to the :class:`Fi` attribute specified by "sort", if the
            attribute is not present the item is skipped.
        :param reverse: bool, set True to reverse sort order
        :param selector: a function which is called with each `Fi` item and has
            to return True (include item) or False (discard item).
            Default function is: ``lambda si: True``. By default only items with
            ``Fi.isValid == True`` are returned.

        :returns: {'attribute1': numpy.array(),
                   'attribute2': numpy.array(),
                   ...
                   }
        """
        selector = (lambda fi: fi.isValid) if selector is None else selector
        attr = attr if attr is not None else []
        attr = set(['id', 'specfile'] + aux.toList(attr))
        items = self.getItems(specfiles, sort, reverse, selector)
        return _getArrays(items, attr, defaultValue)

    def getItems(self, specfiles=None, sort=False, reverse=False,
                 selector=None):
        """Generator that yields filtered and/or sorted :class:`Fi` instances
        from ``self.container``.

        :param specfiles: filenames of ms-run files - if specified return only
            items from those files
        :type specfiles: str or [str, str, ...]
        :param sort: if "sort" is specified the returned list of items is sorted
            according to the :class:`Fi` attribute specified by "sort", if the
            attribute is not present the item is skipped.
        :param reverse: bool, ``True`` reverses the sort order
        :param selector: a function which is called with each ``Fi`` item and
            has to return True (include item) or False (discard item). By
            default only items with ``Fi.isValid == True`` are returned.

        :returns: items from container that passed the selector function
        """
        selector = (lambda fi: fi.isValid) if selector is None else selector
        if specfiles is None:
            specfiles = [_ for _ in viewkeys(self.info)]
        else:
            specfiles = aux.toList(specfiles)
        return _getItems(self.container, specfiles, sort, reverse, selector)

    def addSpecfile(self, specfiles, path):
        """Prepares the container for loading ``fic`` files by adding specfile
        entries to ``self.info``. Use :func:`FiContainer.load()` afterwards
        to actually import the files.

        :param specfiles: the name of an ms-run file or a list of names
        :type specfiles: str or [str, str, ...]
        :param path: filedirectory used for loading and saving ``fic`` files
        """
        for specfile in aux.toList(specfiles):
            if specfile not in self.info:
                self._addSpecfile(specfile, path)
            else:
                warntext = 'Error while calling "FiContainer.addSpecfile()": '\
                           '"%s" is already present "FiContainer.info"'\
                            % (specfile, )
                warnings.warn(warntext)

    def _addSpecfile(self, specfile, path):
        """Adds a new specfile entry to FiContainer.info. See also
        :class:`FiContainer.addSpecfile()`.

        :param specfile: the name of an ms-run file
        :param path: filedirectory used for loading and saving ``fic`` files
        """
        self.info[specfile] = {'path': path}
        self.container[specfile] = dict()

    def setPath(self, folderpath, specfiles=None):
        """Changes the folderpath of the specified specfiles. The folderpath is
        used for saving and loading of ``fic`` files.

        :param specfiles: the name of an ms-run file or a list of names. If None
            all specfiles are selected.
        :type specfiles: None, str, [str, str]
        :param folderpath: a filedirectory
        """
        if specfiles is None:
            specfiles = [_ for _ in viewkeys(self.info)]
        else:
            specfiles = aux.toList(specfiles)

        _containerSetPath(self, folderpath, specfiles)

    def removeSpecfile(self, specfiles):
        """Completely removes the specified specfiles from the ``FiContainer``.

        :param specfiles: the name of an ms-run file or a list of names.
        """
        for specfile in aux.toList(specfiles):
            del self.container[specfile]
            del self.info[specfile]

    def save(self, specfiles=None, compress=True, path=None):
        """Writes the specified specfiles to ``fic`` files on the hard disk.

        .. note::
            If ``.save()`` is called and no ``fic`` files are present in the
            specified path new files are generated, otherwise old files are
            replaced.

        :param specfiles: the name of an ms-run file or a list of names. If None
            all specfiles are selected.
        :param compress: bool, True to use zip file compression
        :param path: filedirectory to which the ``fic`` files are written. By
            default the parameter is set to ``None`` and the filedirectory is
            read from ``self.info[specfile]['path']``
        """
        if specfiles is None:
            specfiles = [_ for _ in viewkeys(self.info)]
        else:
            specfiles = aux.toList(specfiles)

        for specfile in specfiles:
            if specfile not in self.info:
                warntext = 'Error while calling "FiContainer.save()": "%s" is'\
                           ' not present in "FiContainer.info"!'\
                            % (specfile, )
                warnings.warn(warntext)
                continue
            else:
                path = self.info[specfile]['path'] if path is None else path

            with aux.PartiallySafeReplace() as msr:
                filename = specfile + '.fic'
                filepath = aux.joinpath(path, filename)
                with msr.open(filepath) as openfile:
                    self._writeContainer(openfile, specfile, compress)

    def _writeContainer(self, filelike, specfile, compress):
        """Writes the ``self.container`` entry of the specified specfile to the
        ``fic`` format.

        :param filelike: path to a file (str) or a file-like object
        :param specfile: name of an ms-run file present in ``self.info``
        :param compress: bool, True to use zip file compression

        .. note::
            In addition it could also dump the ``self.info`` entry to the
            zipfile with the filename ``info``, but this is not used at the
            moment. For details see :func:`maspy.auxiliary.writeJsonZipfile()`
        """
        aux.writeJsonZipfile(filelike, self.container[specfile],
                             compress=compress
                             )
        #zipcomp = zipfile.ZIP_DEFLATED if compress else zipfile.ZIP_STORED
        #with zipfile.ZipFile(filelike, 'a', allowZip64=True) as containerFile:
        #    infodata = {key: value for key, value in
        #                viewitems(self.info[specfile]) if key != 'path'
        #                }
        #    containerFile.writestr('info', json.dumps(infodata, zipcomp))

    def load(self, specfiles=None):
        """Imports the specified ``fic`` files from the hard disk.

        :param specfiles: the name of an ms-run file or a list of names. If None
            all specfiles are selected.
        :type specfiles: None, str, [str, str]
        """
        if specfiles is None:
            specfiles = [_ for _ in viewkeys(self.info)]
        else:
            specfiles = aux.toList(specfiles)

        for specfile in specfiles:
            if specfile not in self.info:
                warntext = 'Error while calling "FiContainer.load()": "%s" is'\
                           ' not present in "FiContainer.info"!'\
                            % (specfile, )
                warnings.warn(warntext)
                continue
            else:
                fiPath = aux.joinpath(self.info[specfile]['path'],
                                      specfile+'.fic'
                                      )
                with zipfile.ZipFile(fiPath, 'r') as containerZip:
                    #Convert the zipfile data into a str object, necessary since
                    #containerZip.read() returns a bytes object.
                    jsonString = io.TextIOWrapper(containerZip.open('data'),
                                                  encoding='utf-8'
                                                  ).read()
                    #infoString = io.TextIOWrapper(containerZip.open('info'),
                    #                              encoding='utf-8'
                    #                              ).read()
                self.container[specfile] = json.loads(jsonString,
                                                      object_hook=Fi.jsonHook
                                                      )
                #self.info[specfile].update(json.loads(infoString))

    def removeAnnotation(self, specfiles=None):
        """Remove all annotation information from :class:`Fi` elements.

        :param specfiles: the name of an ms-run file or a list of names. If None
            all specfiles are selected.
        :type specfiles: None, str, [str, str]
        """
        if specfiles is None:
            specfiles = [_ for _ in viewkeys(self.info)]
        else:
            specfiles = aux.toList(specfiles)

        for specfile in aux.toList(specfiles):
            for item in viewvalues(self.container[specfile]):
                item.isMatched = False
                item.isAnnotated = False
                item.siIds = list()
                item.siiIds = list()
                item.peptide = None
                item.sequence = None
                item.bestScore = None
