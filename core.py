from __future__ import print_function, division

from collections import defaultdict as ddict
import functools
import itertools
import operator
import os

import cPickle as pickle
import numpy

import maspy.auxiliary as aux
from maspy.auxiliary import lazyAttribute
import maspy.constants
import maspy.peptidemethods


# --- Container and item main classes --- #
class ContainerItem(object):
    """Mass spectrometry data elemtents, derived from specfiles.

    :ivar containerId: used to look up item in :attr:`ItemContainer.index`
    :ivar id: identifier in original file
    :ivar specfile: Keyword (filename) to represent the originating file
    :ivar isValid: this attribute can be used to filter data.
    Should be set to True or False, None if unspecified

    See also :class:`ItemContainer`
    """
    def __init__(self, identifier, specfile):
        self.containerId  = (specfile, identifier)
        self.id = identifier
        self.specfile = specfile
        self.isValid = None

    def __str__(self):
        maxStrLength = max([len(str(key)) for key in self.__dict__.keys()])

        primaryKeys = ['id', 'specfile', 'isValid', 'containerId']
        secondaryKeys = sorted(list(set(self.__dict__.keys()).difference(set(primaryKeys))))

        output = [str(self.__class__)]
        for key in primaryKeys:
            value = getattr(self, key)
            output.append(' '.join([str(key).ljust(maxStrLength), repr(value)]))
        output.append('')
        for key in secondaryKeys:
            value = getattr(self, key)
            output.append(' '.join([str(key).ljust(maxStrLength), repr(value)]))

        return '\n'.join(output)

    def __repr__(self):
        return self.__str__()

    def copy(self):
        """Returns a copy of itself.

        Copies all key, value pairs of self.__dict__, CAUTION: doesn't generate a new instance for values like dict, objects,...
        """
        newItem = self.__class__(self.id, self.specfile)
        newItem.__dict__ = self.__dict__.copy()
        return newItem


class ItemContainer(object):
    """Storage container for :class:`ContainerItem`.

    :ivar index: Use :attr:`ContainerItem.containerId` to look up :class:`ContainerItem` inside the ItemContainer
    :ivar container: Access :class:`ContainerItem` storage list via a specfile keyword: {specfile:[ContainerItem(), ContainerItem(), ...]}
    :ivar specfiles: list of keywords (filenames) representing files
    """
    def __init__(self):
        self.index = dict()
        self.container = dict()
        self.specfiles = list()

    def __getitem__(self, key):
        """Return an item from index using the containerId."""
        return self.index[key]

    def getItems(self, specfiles=None, sort=None, reverse=False, filterAttribute='isValid', filterTargetValue=True, selector=None):
        """Return a filter and/or sorted set of items. By default only valid items are returned.

        :param specfiles: filenames of spectrum files - return only items from those files.
        :type specfiles: str or [str, str, ...]
        :param sort: if sort is specified the returned list of items is sorted according to the
        item attribute specified by sort.
        :param reverse: boolean to reverse sort order

        :param selector: a function which is called with the value of
        the attribute specified by filterAttribute and has to return
        True (include value) or False (discard value). if no selector
        is specified, equality to filterTargetValue is used.

        :param filterAttribute: name of an attribute to use for
        filtering, if attribute does not exist or is set to None
        the value "True" is silently returned -> returns all items

        :param filterTargetValue: if the value of the attribute
        specified by filterAttribute is equal to filterTargetValue the
        element is returned. Is overwritten by specifying selector
        """
        if selector is None:
            if filterAttribute is None:
                selector = lambda x: True
                filterAttribute = str()
            else:
                selector = functools.partial(operator.eq, filterTargetValue)

        specfiles = self.specfiles if specfiles == None else aux.toList(specfiles)
        items = list()
        for specfile in specfiles:
            items.extend(self.container[specfile])

        if sort is not None:
            sortArray = list()
            for item in items:
                sortValue = getattr(item, sort, 0)
                sortArray.append((sortValue, item))
            items = [item for sortValue, item in sorted(sortArray, reverse=reverse)]

        for item in items:
            attributeValue = getattr(item, filterAttribute, True)
            if selector(attributeValue):
                yield item

    def getArrays(self, attributes, specfiles=None, sort=None, reverse=False, filterAttribute='isValid', filterTargetValue=True, selector=None):
        """Return a condensed array of data selected from ContainerItems for faster data processing.

        :param attributes: list of item attributes that should be written to the returned array.
        for the other parameters see :func:`ItemContainer.getItems`

        :returns: dict(key1 from keylist: numpy.array, key2 from keylist: numpy.array, ..., indexPos: numpy.array, id: numpy.array, specfile: numpy.array),
        i.e. returns the columns of the table specified by the list of keys, each numpy.array has the dimensions Nx1. If a value is not present, None, is substituted.
        """
        items = self.getItems(specfiles=specfiles, sort=sort, reverse=reverse,
                              filterAttribute=filterAttribute, filterTargetValue=filterTargetValue,
                              selector=selector
                              )
        arrays = dict()
        attributes = set(['containerId', 'id', 'specfile'] + aux.toList(attributes))
        for key in attributes:
            arrays[key] = list()

        for item in items:
            for key in attributes:
                arrays[key].append(getattr(item, key, None))

        for key in arrays.keys():
            if isinstance(arrays[key][0], tuple):
                emptyArray = numpy.empty(len(arrays[key]), dtype='object')
                emptyArray[:] = arrays[key]
                arrays[key] = emptyArray
            else:
                if type(arrays[key][0]) is float:
                    arrays[key] = numpy.array(arrays[key], dtype='float64')
                elif type(arrays[key][0]) is int:
                    arrays[key] = numpy.array(arrays[key], dtype='int64')
                else:
                    arrays[key] = numpy.array(arrays[key])
        return arrays

    def save(self, filefolder, filename):
        """Store a pickled version of self using  __class__.__name__.lower() as file-appendix.

        :ivar filefolder: folder where the container should be saved
        :ivar filename: filename to store the container, should not contain any appendix like .txt
        """
        #Note: delete the index for storing containers and rewrite upon import to prevent duplication of ContainerItem instances.
        tempIndex = self.index
        del(self.index)
        self._save(filefolder, filename)
        self.index = tempIndex

    def _save(self, filefolder, filename):
        filename = '.'.join((filename, self.__class__.__name__.lower()))
        filepath = aux.joinpath(filefolder, filename)
        with open(filepath, 'w') as openFile:
            pickle.dump(self, openFile)

    @classmethod
    def load(cls, filefolder, filename):
        """Load a pickled version of self using __class__.__name__.lower() as file-appendix.

        :ivar filefolder: folder where the container has been saved
        :ivar filename: filename of the stores container, without file appendix
        """
        classInstance = cls._load(filefolder, filename)
        #Note: delete the index for storing containers and rewrite upon import to prevent duplication of ContainerItem instances.
        classInstance.index = dict()
        for items in classInstance.container.values():
            for item in items:
                classInstance.index[item.containerId] = item
        return classInstance

    @classmethod
    def _load(cls, filefolder, filename):
        filename = '.'.join((filename, cls.__name__.lower()))
        filepath = aux.joinpath(filefolder, filename)
        with open(filepath, 'r') as openFile:
            return pickle.load(openFile)

    def __str__(self):
        numSpecfiles = len(self.specfiles)
        numItems = sum([len(container) for container in self.container.values()])

        output = [str(self.__class__)]
        output.append(' '.join(['Containing', str(numSpecfiles), 'specfiles and', str(numItems), 'items.']))
        for specfile in self.specfiles:
            output.append(''.join([specfile, ', ', str(len(self.container[specfile])), ' items.']))
        return '\n'.join(output)

    def __repr__(self):
        return self.__str__()

    def itemStats(self):
        """Prints the names and the number of occurences of all item attributes that are set
        in items stored in the container instance. """
        itemAttributes = ddict(int)
        for item in self.getItems(filterAttribute=None):
            for key in item.__dict__.keys():
                itemAttributes[key] += 1
        attributeNames = itemAttributes.keys() + ['Attribute name']
        maxStrLength = max([len(str(key)) for key in attributeNames])

        itemCounts = ['  '.join(['Attribute name'.ljust(maxStrLength), 'item counts']),
                      '  '.join(['--------------'.ljust(maxStrLength), '-----------'])
                      ]
        itemCounts.extend(['  '.join([key.ljust(maxStrLength), str(itemAttributes[key])]) for key in sorted(itemAttributes.keys())])
        print('\n'.join(itemCounts))


# --- Auxiliary functions for container class --- #
def addContainer(baseContainer, *newContainers):
    """Merge the content of multiple instances of :class:`ItemContainer` or its subclasses, containers must be of same type.

    :param baseContainer: append newContainers to the baseContainer, has to a class instance
    :param newContainer: one or multiple containers to be appended to the baseContainer

    CAUTION, order of :class:`SpectrumIdentificationItem` in :attr:`SiiContainer.index` can be changed by merging
    #TODO
    """
    for newContainer in newContainers:
        if not isinstance(newContainer, type(baseContainer)):
            print('Cannot combine different container classes, ',
                  repr(baseContainer), ' and ',
                  repr(newContainer)
                  )
            continue

        for specfile in newContainer.specfiles:
            if specfile not in baseContainer.specfiles:
                baseContainer.specfiles.append(specfile)
                baseContainer.container[specfile] = list()

                for item in newContainer.container[specfile]:
                    newItem = item.copy()
                    baseContainer.container[specfile].append(newItem)

                    if baseContainer.__class__.__name__ == 'SiiContainer':
                        if newItem.containerId not in baseContainer.index:
                            baseContainer.index[newItem.containerId] = list()
                        baseContainer.index[newItem.containerId].append(newItem)
                    else:
                        baseContainer.index[newItem.containerId] = newItem

                    if baseContainer.__class__.__name__ == 'SiContainer':
                        if len(newContainer.ionLists) > 0:
                            baseContainer.ionLists[newItem.containerId] = dict(newContainer.ionLists[newItem.containerId])
            else:
                print(specfile, 'already present in baseContainer.')
    return baseContainer


# --- Container and item subclasses --- #
class SpectrumItem(ContainerItem):
    """Representation of a spectrum."""
    def __init__(self, identifier, specfile):
        super(SpectrumItem, self).__init__(identifier, specfile)
        self.msLevel = None


class SiContainer(ItemContainer):
    """ItemContainer for mass spectrometry data (eg Ms1 and Ms2 spectra),
    SiContainer = Spectrum Item Container.

    for parameter and method description see :class:`ItemContainer`
    see also :class:`SiiContainer` (Spectrum Identification Item Container) which contains sequence data.

    :ivar ionLists: spectrum ion m/z and intensity information, not loaded by default
    dict(containerId=dict(mz=numpy.array([mass / charge, ...]), i=numpy.array([intensity, ...]))).
    """
    def __init__(self):
        super(SiContainer, self).__init__()
        self.ionLists = dict()

    def save(self, filefolder, filename, saveIonList=True):
        """Store a pickled version of the self, using '.SiContainer' as file-appendix.

        Stores the ionList in a separate file with appendix '.ionlist'.
        """
        try:
            ionLists = self.ionLists
            del self.ionList
            super(self.__class__, self).save(filefolder, filename)
        finally:
            self.ionLists = ionLists

        if saveIonList and len(ionLists) != 0:
            keyList = list()
            mzList = list()
            iList = list()
            for key, value in self.ionLists.items():
                keyList.append(key)
                mzList.append(value['mz'])
                iList.append(value['i'])
            ionArray = numpy.array([keyList, numpy.array(mzList), numpy.array(iList)])

            ionListFileName = '.'.join((filename, 'ionlist'))
            ionListFilePath = os.path.join(filefolder, ionListFileName).replace('\\', '/')
            with open(ionListFilePath, 'wb') as openFile:
                numpy.save(openFile, ionArray)

    @classmethod
    def load(cls, filefolder, filename, importIonList=True):
        """Load a pickled version of the self, using the __class__.__name__ as file-appendix."""
        siContainer = super(cls, cls).load(filefolder, filename)
        siContainer.ionLists = dict()

        ionListFileName = '.'.join((filename, 'ionlist'))
        ionListFilePath = os.path.join(filefolder, ionListFileName).replace('\\', '/')
        if importIonList and os.path.isfile(ionListFilePath):
            importedArray = numpy.load(ionListFilePath)
            for key, mzList, iList in itertools.izip(importedArray[0], importedArray[1], importedArray[2]):
                siContainer.ionLists[key] = {'mz': mzList, 'i':iList}
        return siContainer


class SpectrumIdentificationItem(ContainerItem):
    """Representation of a peptide fragment spectrum annotation (Peptide Spectrum Match)."""
    def __init__(self, identifier, specfile):
        super(SpectrumIdentificationItem, self).__init__(identifier, specfile)
        self.diPeptide = None


class SiiContainer(ItemContainer):
    """ItemContainer for msn spectrum identifications (Peptide Spectrum Matches),
    SiiContainer = Spectrum Identification Item Container.

    for parameter and method description see :class:`ItemContainer`
    see also :class:`SiContainer` (Spectrum Item Container) which contains spectrum data.
    """
    def __init__(self):
        super(SiiContainer, self).__init__()

    def addSiInfo(self, siContainer, specfiles=None, attributes=['obsMz', 'rt', 'charge']):
        """ Copy attributes into sii from the corresponding SpectrumItem in siContainer,
        if an attribute is not presend in the SpectrumItem the attribute value is set to None
        Attribute examples: 'obsMz', 'rt', 'charge', 'tic', 'iit', 'ms1Id'
        """
        specfiles = self.specfiles if specfiles == None else aux.toList(specfiles)

        for specfile in specfiles:
            if specfile not in self.specfiles:
                print(specfile, 'not present in siiContainer.')
            elif specfile not in siContainer.specfiles:
                print(specfile, 'not present in siContainer.')
            else:
                for sii in self.container[specfile]:
                    si = siContainer.index[sii.containerId]
                    for attribute in attributes:
                        setattr(sii, attribute, getattr(si, attribute, None))

    def getValidItem(self, key):
        """Returns one or a tuple of only valid items from index, usually the rank1 PSM.

        :param key: is the value of :attr:`SpectrumIdentificationItem.containerId` to retrieve
        """
        if not isinstance(key, tuple):
            key = tuple(key)
        items = [item for item in self.index[key] if item.isValid]
        items = items[0] if len(items) == 1 else tuple(items)
        return items

    def calcMz(self, specfiles=None, guessCharge=True):
        # Guess charge uses the calculated mass and the observed m/z value to calculate the charge
        specfiles = self.specfiles if specfiles is None else aux.toList(specfiles)
        tempPeptideMasses = dict()
        for specfile in specfiles:
            if specfile not in self.specfiles:
                print(specfile, 'not present in siiContainer.')
            else:
                for sii in self.getItems(specfiles=specfile, filterAttribute=None):
                    charge = sii.charge
                    peptide = sii.peptide
                    if peptide not in tempPeptideMasses:
                        if sii.diPeptide:
                            tempPeptideMasses[peptide] = (maspy.peptidemethods.calcPeptideMass(sii.peptide1) +
                                                          maspy.peptidemethods.calcPeptideMass(sii.peptide2))
                        else:
                            tempPeptideMasses[peptide] = maspy.peptidemethods.calcPeptideMass(peptide)
                    peptideMass = tempPeptideMasses[peptide]
                    if charge is not None:
                        sii.calcMz = maspy.peptidemethods.calcMzFromMass(peptideMass, charge)
                    elif guessCharge:
                        guessedCharge = round(peptideMass / (sii.obsMz - maspy.constants.atomicMassProton), 0)
                        sii.calcMz = maspy.peptidemethods.calcMzFromMass(peptideMass, guessedCharge)
                        sii.charge = guessedCharge
        del(tempPeptideMasses)

    @classmethod
    def load(cls, filefolder, filename):
        """Load a pickled version of self using __class__.__name__.lower() as file-appendix.

        :ivar filefolder: folder where the container has been saved
        :ivar filename: filename of the stores container, without file appendix
        """
        classInstance = cls._load(filefolder, filename)
        classInstance.index = dict()
        for items in classInstance.container.values():
            for item in items:
                if item.containerId not in classInstance.index:
                    classInstance.index[item.containerId] = list()
                classInstance.index[item.containerId].append(item)
        return classInstance


class FeatureItem(ContainerItem):
    """Representation of a peptide LC-MS feature.

    :ivar isMatched: None if unspecified, should be set to False on import, True if any Si or Sii elements could be matched
    :ivar isAnnotated: None if unspecified, should be set to False on import, True if any Sii elements could be matched
    :ivar siIds: containerId values of matched Si entries
    :ivar siiIds: containerId values of matched Sii entries
    :ivar peptide: peptide sequence of best scoring Sii match
    :ivar sequence: plain amino acid sequence of best scoring Sii match, used to retrieve protein information
    :ivar score: score of best scoring Sii match
    """
    def __init__(self, identifier, specfile):
        super(FeatureItem, self).__init__(identifier, specfile)
        self.isMatched = None
        self.isAnnotated = None
        self.siIds = list()
        self.siiIds = list()
        self.peptide = None
        self.sequence = None
        self.score = None


class FeatureContainer(ItemContainer):
    """ItemContainer for peptide elution features :class`FeatureItem`.

    for parameter and method description see :class:`ItemContainer`
    see also :class:`SiContainer` (Spectrum Item Container) which contains spectrum data.
    see also :class:`SiiContainer` (Spectrum Identification Item Container) which contains sequence data.
    """
    def __init__(self):
        super(FeatureContainer, self).__init__()

    def removeAnnotation(self):
        """Remove all annotation information from :class:`FeatureItem` in :class:`FeatureContainer`."""
        for items in self.container.values():
            for item in items:
                item.isMatched = False
                item.isAnnotated = False
                item.siIds = list()
                item.siiIds = list()
                item.peptide = None
                item.sequence = None
                item.score = None
