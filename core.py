from __future__ import print_function

import bisect
import csv
import cPickle as pickle
from collections import defaultdict as ddict
import functools
import itertools
import operator
import os
import re

import numpy
from lxml import etree as etree
from matplotlib import pyplot as plt

import pymzml
import pyteomics.mass

import pyms.auxiliary as aux


def lazyAttribute(fn):
    attributeName = '_lazy_' + fn.__name__
    @property
    def _lazyAttribute(self):
        if not hasattr(self, attributeName):
            setattr(self, attributeName, fn(self))
        return getattr(self, attributeName)
    return _lazyAttribute


class ContainerItem(object):
    """Mass spectrometry data elemtents, derived from specfiles.

    :ivar containerId: used to look up item in :attr:`ItemContainer.index`
    :ivar id: identifier in original file
    :ivar specfile: Keyword (filename) to represent the originating file
    :ivar isValid: should evaluate to True or False, None if unspecified - used to filter data

    See also :class:`ItemContainer`
    """
    def __init__(self, identifier, specfile):
        self.containerId  = (specfile, identifier)
        self.id = identifier
        self.specfile = specfile
        self.isValid = None

    def copy(self):
        """Returns a copy of itself.

        Copies all key, value pairs of self.__dict__, doesn't generate a new instance for values like dict, objects,...
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
        :param sort: if sort is None items are returned in the import order, otherwise the items are sorted according to the
        item attribute specified by sort
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
                arrays[key] = numpy.array(arrays[key])

        return arrays

    def save(self, filefolder, filename):
        """Store a pickled version of self using  __class__.__name__.lower() as file-appendix.

        :ivar filefolder: folder where the container should be saved
        :ivar filename: filename to store the container, should not contain any appendix like .txt
        """
        del(self.index)
        self._save(filefolder, filename)

    def _save(self, filefolder, filename):
        filename = '.'.join((filename, self.__class__.__name__.lower()))
        filepath = os.path.join(filefolder, filename).replace('\\', '/')
        with open(filepath, 'w') as openFile:
            pickle.dump(self, openFile)

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
                classInstance.index[item.containerId] = item
        return classInstance

    @classmethod
    def _load(cls, filefolder, filename):
        filename = '.'.join((filename, cls.__name__.lower()))
        filepath = os.path.join(filefolder, filename).replace('\\', '/')
        with open(filepath, 'r') as openFile:
            return pickle.load(openFile)


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

    :ivar ionLists: spectrum ion information, not loaded by default
    dict(containerId=dict(mz=nump.array([mass / charge, ...]), i=numpy.array([intensity, ...]))).
    """
    def __init__(self):
        super(SiContainer, self).__init__()
        self.ionLists = dict()

    def save(self, filefolder, filename, saveIonList=True):
        """Store a pickled version of the self, using '.SiContainer' as file-appendix.

        Stores the ionList in a separate file with appendix '.ionlist'.
        """
        ionLists = self.ionLists
        del(self.ionLists)
        try:
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
                siContainer.ionLists[key] = dict()
                siContainer.ionLists[key]['mz'] = mzList
                siContainer.ionLists[key]['i'] = iList
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
        items = list()
        for item in self.index[key]:
            if item.isValid:
                items.append(item)
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
                for sii in self.getItems(specfiles=specfile):
                    charge = sii.charge
                    peptide = sii.peptide
                    if peptide not in tempPeptideMasses:
                        if sii.diPeptide:
                            tempPeptideMasses[peptide] = calcPeptideMass(sii.peptide1) + calcPeptideMass(sii.peptide2)
                        else:
                            tempPeptideMasses[peptide] = calcPeptideMass(peptide)
                    peptideMass = tempPeptideMasses[peptide]
                    if charge is not None:
                        sii.calcMz = aux.calcMzFromMass(peptideMass, charge)
                    elif guessCharge:
                        guessedCharge = round(peptideMass / (sii.obsMz - aux.atomicMassProton), 0)
                        sii.calcMz = aux.calcMzFromMass(peptideMass, guessedCharge)
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
    """Representation of a peptide elution feature.

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


class FeatureGroupItem(ContainerItem):
    """Representation of a group of :class:`FeatureItem`.

    :ivar isMatched: False by default, True if any :class:`FeatureItem` in the group are matched.
    :ivar isAnnotated: False by default, True if any :class:`FeatureItem` in the group are annotated.
    :ivar siIds: containerId values of matched Si entries
    :ivar siiIds: containerId values of matched Sii entries
    :ivar featureIds: containerId values of :class:`FeatureItem` in the feature group
    :ivar peptide: peptide sequence of best scoring Sii match
    :ivar sequence: plain amino acid sequence of best scoring Sii match, used to retrieve protein information
    :ivar score: score of best scoring Sii match
    :ivar matchMatrix: structured representation of :attr:`FeatureItem.containerId` in the feature group.
    :ivar intensityMatrix: similar to :attr:`matchMatrix` but contains :attr:`FeatureItem.intensity` values.
    {chargeState: 2d numpy.array with specfiles as 1st dimension and labelState as 2nd dimension}

    See also :class:`ItemContainer`
    """
    def __init__(self):
        super(FeatureGroupItem, self).__init__(None, None)
        del(self.id)
        del(self.specfile)
        del(self.containerId)
        self.isMatched = None
        self.isAnnotated = None
        self.siIds = list()
        self.siiIds = list()
        self.featureIds = list()
        self.peptide = None
        self.sequence = None
        self.score = None
        self.matchMatrix = dict()
        self.intensityMatrix = dict()


class FeatureGroupContainer(ItemContainer):
    """ItemContainer for peptide feature groups :class`FeatureGroupItem`.

    :ivar container: Storage list of :class:`FeatureGroupItem`
    :ivar index: Use :attr:`FeatureItem.containerId` to which :class:`FeatureGroupItem` the feature was grouped
    :ivar labelDescriptor: :class:`LabelDescriptor` describes the label setup of an experiment
    :ivar specfiles: List of keywords (filenames) representing files
    :ivar specfilePositions: {specfile:arrayPosition, ...}
    arrayPosition respresents the array position of a specfile in :attr:`FeatureGroupItem.matchMatrix`

    for further parameter and method description see :class:`ItemContainer`
    see also :class:`SiContainer` (Spectrum Item Container) which contains spectrum data.
    see also :class:`SiiContainer` (Spectrum Identification Item Container) which contains sequence data.
    see also :class:`FeatureContainer` (Feature Container) which contains peptide features.
    """
    def __init__(self, specfiles, labelDescriptor=None):
        super(FeatureGroupContainer, self).__init__()
        self.specfiles = specfiles
        self.container = list()
        self.labelDescriptor = LabelDescriptor() if labelDescriptor is None else labelDescriptor

        self.specfilePositions = dict()
        for position, specfile in enumerate(self.specfiles):
            self.specfilePositions[specfile] = position

    def getItems(self, sort=None, reverse=False, filterAttribute='isValid', filterTargetValue=True, selector=None):
        """Return a filter and/or sorted set of items. By default only valid items are returned.

        :param sort: if sort is None items are returned in the container order, otherwise the items are sorted according to the
        item attribute specified by sort
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

        items = self.container
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

    def getArrays(self, attributes=['mh', 'rt', 'peptide', 'sequence'], sort=None, reverse=False,
                  filterAttribute='isValid', filterTargetValue=True, selector=None, report='lfq'
                  ):
        """Return a condensed array of data selected from ContainerItems for faster data processing.

        :param attributes: ('score', 'peptide', 'sequence', 'siIds', 'siiIds', 'featureIds', 'rt', 'mh', 'isAnnotated', 'isMatched', 'isValid')
        :param attributes: list of item attributes that should be written to the returned array.
        for the other parameters see :func:`FeatureGroupContainer.getItems`

        :returns: dict(key1 from keylist: numpy.array, key2 from keylist: numpy.array, ..., indexPos: numpy.array, id: numpy.array, specfile: numpy.array),
        i.e. returns the columns of the table specified by the list of keys, each numpy.array has the dimensions Nx1. If a value is not present, None, is substituted.
        """
        items = self.getItems(sort=sort, reverse=reverse, filterAttribute=filterAttribute,
                              filterTargetValue=filterTargetValue, selector=selector
                              )
        if report.lower() in ['lfq', 'labelfree']:
            report = 'lfq'
        elif report.lower() in ['label', 'labeled', 'sil']:
            report = 'sil'

        arrays = dict()
        for key in attributes:
            arrays[key] = list()

        reportAttributes = list()
        if report == 'lfq':
            arrays['charge'] = list()
            arrays['labelState'] = list()
            for specfile in self.specfiles:
                arrays[specfile] = list()
                reportAttributes.append(specfile)
        elif report == 'sil':
            arrays['charge'] = list()
            arrays['specfile'] = list()
            for labelState in self.labelDescriptor.labels.keys() + [-1]:
                labelAttributeName = ' '.join(('label:', str(labelState)))
                arrays[labelAttributeName] = list()
                reportAttributes.append(labelAttributeName)

        for item in items:
            if report == 'sil':
                for charge in item.intensityMatrix.keys():
                    for specfile, specfilePosition in self.specfilePositions.items():
                        for key in attributes:
                            arrays[key].append(getattr(item, key, None))
                        arrays['charge'].append(charge)
                        arrays['specfile'].append(specfile)
                        for labelState in self.labelDescriptor.labels.keys() + [-1]:
                            labelAttributeName = ' '.join(('label:', str(labelState)))
                            arrays[labelAttributeName].append(item.intensityMatrix[charge][specfilePosition, labelState])
            if report == 'lfq':
                for charge in item.intensityMatrix.keys():
                    for labelState in self.labelDescriptor.labels.keys() + [-1]:
                        for key in attributes:
                            arrays[key].append(getattr(item, key, None))
                        arrays['charge'].append(charge)
                        arrays['labelState'].append(labelState)
                        for specfile, specfilePosition in self.specfilePositions.items():
                            arrays[specfile].append(item.intensityMatrix[charge][specfilePosition, labelState])

        for key in arrays.keys():
            if key in reportAttributes:
                arrays[key] = numpy.array(arrays[key], dtype='float64')
            else:
                arrays[key] = numpy.array(arrays[key])
        return arrays

    @classmethod
    def load(cls, filefolder, filename):
        """Load a pickled version of self using __class__.__name__.lower() as file-appendix.

        :ivar filefolder: folder where the container has been saved
        :ivar filename: filename of the stores container, without file appendix
        """
        classInstance = super(cls, cls)._load(filefolder, filename)
        classInstance.index = dict()
        for item in classInstance.container:
            for featureId in item.featureIds:
                if featureId not in classInstance.index:
                    classInstance.index[featureId] = list()
                classInstance.index[featureId].append(item)
        return classInstance


class PeptideSequence(object):
    """Describes a peptide as derived by digestion of one or multiple proteins, can't contain any modified amino acids.
    see also :class:`PeptideEvidence`

    :ivar sequence: amino acid sequence of the peptide
    :ivar missedCleavage: number of missed cleavages, dependens on enzyme specificity
    :ivar proteinList: protein ids that generate this peptide under certain digest condition
    :ivar proteinPositions: {proteinId:(startPosition, endPositions) ...} startposition and endposition of peptide in protein
    eg. 'AADITSLYK' IN 'TAKAADITSLYKEETR':(4,12)
    :ivar mass: peptide mass in Daltons
    :ivar length: number of amino acids
    """
    def __init__(self, sequence, mc=None):
        self.sequence = sequence

        self.missedCleavage = mc
        self.isUnique = None
        self.proteins = set()
        self.proteinPositions = dict()

    @lazyAttribute
    def length(self):
        return len(self.sequence)
    @lazyAttribute
    def mass(self):
        return pyteomics.mass.calculate_mass(self.sequence, charge=0)


class PeptideEvidence(PeptideSequence):
    """Summarizes all the evidence (:class:`SpectrumIdentificationItem`) for a certain peptide.
    for other parameters see :class:`PeptideSequence`

    :ivar peptide: amino acid sequence of the peptide including modifications (written in brackets 'AADIT[modification]SLYK')
    :ivar sequence: amino acid sequence of the peptide, corresponds to peptideRef of mzidentml files
    :ivar bestId: containerId of best scoring Sii item
    :ivar siiIds: containerIds of all added Sii items
    :ivar score: best score of all added Sii items
    :ivar scores: scores of all added Sii items
    """
    def __init__(self, peptide, sequence=None):
        sequence = removeModifications(peptide) if sequence is None else sequence
        super(PeptideEvidence, self).__init__(sequence)

        self.peptide = peptide
        self.sequence = sequence
        self.bestId = tuple()
        self.score = float()
        self.siiIds = list()
        self.scores = list()

    @lazyAttribute
    def mass(self):
        return calcPeptideMass(peptide)

    @classmethod
    def fromPeptideSequence(cls, peptide, peptideSequence):
        newInstance = cls(peptide, peptideSequence.sequence)
        newInstance.isUnique = peptideSequence.isUnique
        newInstance.missedCleavage = peptideSequence.missedCleavage
        newInstance.proteins = peptideSequence.proteins
        newInstance.proteinPositions = peptideSequence.proteinPositions
        return newInstance


class ProteinSequence(object):
    """Describes a protein.

    :ivar id: identifier of the protein eg. UniprotId
    :ivar name: name of the protein
    :ivar sequence: amino acid sequence of the protein
    :ivar isUnique: boolean, True if at least one unique peptide can be assigned to the protein
    :ivar uniquePeptides: a set of peptides which can be unambiguously assigned to this protein
    :ivar sharedPeptides: a set of peptides which are shared between different proteins
    :ivar mass: protein mass in Daltons
    :ivar length: number of amino acids
    :ivar coverageUnique: the number of amino acids in the protein sequence that are coverd by unique peptides
    :ivar coverageShared: the number of amino acids in the protein sequence that are coverd by unique or shared peptides
    """
    def __init__(self, identifier, sequence, name=str()):
        self.id = identifier
        self.name = name
        self.sequence = sequence

        self.isUnique = None
        self.uniquePeptides = set()
        self.sharedPeptides = set()

    @lazyAttribute
    def mass(self):
        return pyteomics.mass.calculate_mass(self.sequence, charge=0)

    @lazyAttribute
    def length(self):
        return len(self.sequence)


class ProteinEvidence(ProteinSequence):
    """ Summarizes all the PeptideEvidence information for a certain protein
    for other paremeters see :class:`ProteinSequence`

    :ivar uniquePsmCount: the sum of PSMs of all unique peptides
    :ivar sharedPsmCount: the sum of PSMs of all shared peptides
    :ivar isValid: should evaluate to True or False, None if unspecified - used to filter data.
    """
    def __init__(self, identifier, sequence, name=None):
        super(ProteinEvidence, self).__init__(identifier, sequence, name)
        self.uniquePsmCount = int()
        self.sharedPsmCount = int()
        self.isValid = None


class ProteinDatabase(object):
    """Describes proteins and peptides derived by an in silico digest.

    :ivar peptides: {sequence:PeptideSequence(), ...} contains elements of :class:`PeptideSequence` derived by an in silico digest of the proteins
    :ivar proteins: {proteinId:Protein(), proteinId:Protein()}, used to access :class:`ProteinSequence` elements by their id
    :ivar proteinNames: {proteinName:Protein(), proteinName:Protein()}, alternative way to access :class:`ProteinSequence` elements by their names
    """
    # A container for Protein or ProteinEvidence objects
    def __init__(self):
        self.peptides = dict()
        self.proteins = dict()
        self.proteinNames = dict()

    def __getitem__(self, key):
        """Uses key to return protein entries :class:`Protein`.

        :ivar key: either a protein id or a protein name
        """
        if key in self.proteins:
            return self.proteins[key]
        elif key in self.proteinNames:
            return self.proteinNames[key]
        else:
            raise KeyError(key)

    def calculateCoverage(self):
        """Calcualte the sequence coverage masks for all ProteinEvidence() elements.

        For detailed description see :func:`_calculateCoverageMasks`
        """
        _calculateCoverageMasks(self.proteins, self.peptides)


class EvidenceContainer(object):
    """Container to collect peptide evidence from :class:`SiiContainer` and summarize to protein evidence.

    :ivar db: :class:`ProteinDatabase`, representation of an in silico digest of a fasta file.
    :ivar proteinEvidence: {proteinId: :class:`ProteinEvidence`, ...}
    :ivar peptideEvidence: {peptide: :class:`PeptideEvidence`, ...}

    :ivar uniqueProteins: list of protein ids which have at least one unique peptideEvidence entry.
    :ivar scoreKey: score attribtue name of :class:`SiiItem`
    :ivar largerBetter: boolean, True if a larger score is better
    """
    def __init__(self, proteinDatabase):
        self.db = proteinDatabase
        self.proteinEvidence = dict()
        self.peptideEvidence = dict()

        self.uniqueProteins = list()
        self.scoreKey = None
        self.largerBetter = None

    def calculateCoverage(self):
        """Calcualte the sequence coverage masks for all ProteinEvidence() elements.

        For detailed description see :func:`_calculateCoverageMasks`
        """
        _calculateCoverageMasks(self.proteinEvidence, self.peptideEvidence)


#####################################################
### Auxiliary functions for ItemContainer classes ###
#####################################################
def addContainer(baseContainer, *newContainers):
    """Merge the content of multiple instances of :class:`ItemContainer` or its subclasses, containers must be of same type.

    :param baseContainer: append newContainers to the baseContainer, can be an class or an instance
    :param newContainer: one or multiple containers to be appended to the baseContainer

    Caution, order of :class:`SpectrumIdentificationItem` in :attr:`SiiContainer.index` can be changed by merging
    """
    if isinstance(baseContainer, type):
        baseContainer = baseContainer()
    containerClass = baseContainer.__class__

    for newContainer in newContainers:
        if not isinstance(newContainer, baseContainer.__class__):
            print('Cannot combine different container classes, ',
                  baseContainer.__class__.__name__, ' and ',
                  newContainer.__class__.__name__
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
                            baseContainer.ionLists[newItem.containerId] = dict()
                            baseContainer.ionLists[newItem.containerId]['i'] = newContainer.ionLists[newItem.containerId]['i']
                            baseContainer.ionLists[newItem.containerId]['mz'] = newContainer.ionLists[newItem.containerId]['mz']
            else:
                print(specfile, 'already present in baseContainer.')
    return baseContainer


def importSpecfiles(specfiles, fileDirectory, importIonList=False, siContainer=None):
    """Auxiliary function to conveniently import a group of specfiles.

    :ivar specfiles: Filenames which should be imported
    :type specfiles: str() or [str(), str(), ...]
    :ivar fileDirectory: Filenames are searched in this folder and its folder
    :ivar bool importIonList: True if ion arrays (mz and intensity) should be imported
    :ivar siContainer: Add imported specfiles to siContainer, if None a new instance of :class:`SiContainer` is returned
    """
    if siContainer is None:
        siContainer = SiContainer()
    specfiles = aux.toList(specfiles)
    for specfile in specfiles:
        specfilePath = aux.searchFileLocation(specfile, 'sicontainer', fileDirectory)
        if specfilePath:
            fileFolder = os.path.dirname(specfilePath)
            addContainer(siContainer, SiContainer.load(fileFolder, specfile, importIonList=importIonList))
        else:
            specfilePath = aux.searchFileLocation(specfile, 'mzML', fileDirectory)
            if specfilePath is None:
                print('File not found: ', specfile)
            else:
                importSpectrumItems(siContainer, specfilePath, specfile, importIonList=importIonList)
    return siContainer


def generateSiContainerFiles(fileDirectory, report=True):
    """Generate .sicontainer and .ionlist files for all mzML files in the fileDirectory and its subfolders.

    see also :func:`removeSiContainerFiles` and :meth:`SiContainer.save`
    """
    for filePath in aux.matchingFilePaths('', fileDirectory, targetFileExtension='mzML', selector=lambda x: True):
        dotPosition = [x for x in aux.findAllSubstrings(filePath, '.')][-1]

        fileFolder = os.path.dirname(filePath)
        fileName = os.path.basename(filePath[:dotPosition])
        targetFilePath = '.'.join((filePath[:dotPosition], 'sicontainer'))
        if not os.path.isfile(targetFilePath):
            siContainer = SiContainer()
            importSpectrumItems(siContainer, filePath, fileName, importIonList=True)
            if report:
                print('Saving SiContainer() / ionList :', fileName)
            siContainer.save(fileFolder, fileName)


def removeSiContainerFiles(fileDirectory, report=True):
    """Remove all .sicontainer and .ionlist files in the fileDirectory and its subfolders.

    :ivar fileDirectory: target directory
    :ivar report: boolean, if True print path of removed files

    see also :func:`generateSiContainerFiles` and :meth:`SiContainer.save`
    """
    for filePath in aux.matchingFilePaths('', fileDirectory, targetFileExtension='sicontainer', selector=lambda x: True):
        os.remove(filePath)
        if report:
            print('removed:\n', filePath)
    for filePath in aux.matchingFilePaths('', fileDirectory, targetFileExtension='ionlist', selector=lambda x: True):
        os.remove(filePath)
        if report:
            print('removed:\n', filePath)


###############################################
### Import functions for various file types ###
###############################################
def pymzmlReadMzml(mzmlPath):
    """Auxiliary function to specify extra accesions when reading an mzml file with :func:`pymzml.run.Reader`.

    :ivar mzmlPath: File path
    """
    extraAccessions = [('MS:1000827', ['value']),
                       ('MS:1000828', ['value']),
                       ('MS:1000829', ['value']),
                       ('MS:1000744', ['value']),
                       ('MS:1000016', ['value']),
                       ('MS:1000927', ['value']),
                       ('MS:1000512', ['value'])
                       ]
    return pymzml.run.Reader(mzmlPath, extraAccessions=extraAccessions)


def importSpectrumItems(siContainer, specfilePath, specfile, msLevel=[1, 2], importIonList=True, mgfType=None):
    """ Import spectra from mzml or mgf files.

    :param siContainer: Spectra are added to to this instance of :class:`SiContainer`
    :param specfilePath: Actual file path (file type has to be .mzml or .mgf)
    :param specfile: Keyword (filename) to represent file in the :class:`SiContainer`. Each filename
    can only occure once, therefore importing the same filename again is prevented
    :param importIonList: bool if ion arrays (mz and intensity) should be imported
    :param mgfType: if the file is of type '.mgf', and the mgf was generated by pParse set to mgfType="pParse" because of header information ambiguity.
    """
    if specfilePath is None:
        raise Exception('SpecfilePath not specified, pyms.core.importSpectrumItems()') #TODO
    if not os.path.isfile(specfilePath):
        print('File does not exits:', specfilePath)
    elif not specfilePath.lower().endswith('.mzml') and not specfilePath.lower().endswith('.mgf'):
        print('File is not an "mzml" or "mgf" file:', specfilePath)
    else:
        if specfile not in siContainer.specfiles:
            siContainer.specfiles.append(specfile)
            siContainer.container[specfile] = list()
            if specfilePath.lower().endswith('.mzml'):
                msrun = pymzmlReadMzml(specfilePath)
                importMzmlSpectrumItems(siContainer, msrun, specfile, importIonList=importIonList)
            elif specfilePath.lower().endswith('.mgf'):
                importMgfSpectrumItems(siContainer, specfilePath, specfile, importIonList=importIonList, mgfType=mgfType)
        else:
            print(specfile, 'is already present in the SiContainer, import interrupted.')


def importMzmlSpectrumItems(siContainer, msrun, specfile, importIonList=False):
    """Import spectrum information from mzML files. Mostly used as a private function by :func:`importSpectrumItems`.

    :param siContainer: Spectra are added to this instance of :class:`SiContainer`
    :param msrun: :object:`pymzml.run.Reader` iterator for mzml spectra of the module :mod:`pymzml`
    :param specfile: Keyword (filename) to represent file in the :class:`SiContainer`
    :param importIonList: bool if ion arrays (mz and intensity) should be imported
    """
    def _generateSpectrumItem(spectrum):
        si = SpectrumItem(str(spectrum['id']), specfile)
        si.msLevel = spectrum['ms level']
        si.isValid = True
        for attributeName, accession, msLevel in [('iit', 'MS:1000927', None),
                                                  ('tic', 'MS:1000285', None),
                                                  ('rt', 'MS:1000016', None),
                                                  ('targetWindowMz', 'MS:1000827', 2),
                                                  ('lowWindowOffset', 'MS:1000828', 2),
                                                  ('highWindowOffset', 'MS:1000829', 2),
                                                  ('obsMz', 'MS:1000744', 2),
                                                  ('charge', 'MS:1000041', 2),
                                                  ('filterString', 'MS:1000512', None)
                                                  ]:
            if msLevel is None or si.msLevel == msLevel:
                if accession in spectrum:
                    setattr(si, attributeName, spectrum[accession])
                    if attributeName == 'rt':
                        #change from minutes to seconds
                        setattr(si, attributeName, spectrum[accession] * 60)
                else:
                    setattr(si, attributeName, None)
        return si

    currMsnContainerIdList = list()
    for spectrum in msrun:
        if spectrum['ms level'] >= 1:
            si = _generateSpectrumItem(spectrum)
            siContainer.container[specfile].append(si)
            siContainer.index[si.containerId] = si

            if importIonList:
                siContainer.ionLists[si.containerId] = dict()
                siContainer.ionLists[si.containerId]['mz'] = numpy.array(spectrum.mz, dtype='float64')
                siContainer.ionLists[si.containerId]['i'] = numpy.array(spectrum.i, dtype='float64')

            if si.msLevel == 1:
                # Add currMsnIndexList to last ms1 entry and reset currMsnContainerIdList
                setattr(si, 'msnIdList', currMsnContainerIdList)
                currMsnContainerIdList = list()
                # Change lastMs1Item to current SpectrumItem
                lastMs1Item = si

            else:
                currMsnContainerIdList.append(si.containerId)
                setattr(si, 'ms1Id', lastMs1Item.containerId)
    setattr(lastMs1Item, 'msnIdList', currMsnContainerIdList)


def importMgfSpectrumItems(siContainer, specfilePath, specfile, importIonList=False, mgfType=None):
    """Import spectrum information from mgf files.

    Mostly used as a private function by :func:`importSpectrumItems`.

    :param siContainer: Spectra are added to this instance of :func:`SiContainer`
    :param specfilePath: Actual path to file
    :param specfile: Keyword (filename) to represent file in the :class:`SiContainer`
    :param importIonList: bool if ion arrays (mz and intensity) should be imported
    :param mgfType: If mgf was generated by pParse set to "Parse" because of header information ambiguity
    """
    def _generateSpectrumItem(mgfEntry):
        mgfEntry = mgfEntry.replace('END IONS','').strip()
        ionList = list()
        ionPos = False
        attributes = dict()
        for line in mgfEntry.split('\n'):
            if ionPos:
                ionList.append(line)
            elif line[0] not in ['0','1','2','3','4','5','6','7','8','9']:
                key = line.split('=')[0].lower()
                value = line.split('=')[1]
                attributes[key] = value
            else:
                if importIonList:
                    ionPos = True
                    ionList.append(line)
                else:
                    break

        # Read scanNr from various mgf formats #
        if mgfType == 'pParse':
            nativeScanNr = attributes['title'].split('.')[1]
            scanNrExtension = attributes['title'].split('.')[4]
            scanNr = nativeScanNr+'.'+scanNrExtension
        elif 'SCANS' in attributes:
            scanNr = attributes['scans']
        else:
            scanNr = attributes['title'].split('.')[1]

        si = SpectrumItem(str(scanNr), specfile)
        si.msLevel = 2
        si.isValid = True
        si.charge = int(attributes['charge'].strip('+').strip('-')) if 'charge' in attributes else None
        si.obsMz = float(attributes['pepmass'].split(' ')[0]) if 'pepmass' in attributes else None
        si.rt = float(attributes['rtinseconds']) if 'rtinseconds' in attributes else None
        return si, ionList

    def _splitMgfIonList(ionList):
        mzList = list()
        iList = list()
        for ionEntry in ionList:
            ionEntry = ionEntry.split(' ')
            if len(ionEntry) == 2:
                mzList.append(ionEntry[0])
                iList.append(ionEntry[1])
        return mzList, iList

    with open(specfilePath,'rb') as mgfFile:
        mgfRead  = mgfFile.read()
        mgfSplit = mgfRead.split('BEGIN IONS\n')
        for mgfEntry in mgfSplit:
            if mgfEntry != '' and mgfEntry.find('END IONS') != -1:
                si, ionList = _generateSpectrumItem(mgfEntry)
                siContainer.container[specfile].append(si)
                siContainer.index[si.containerId] = si

                if importIonList:
                    mzList, iList = _splitMgfIonList(ionList)
                    siContainer.ionLists[si.containerId] = dict()
                    siContainer.ionLists[si.containerId]['mz'] = numpy.array(mzList, dtype='float64')
                    siContainer.ionLists[si.containerId]['i'] = numpy.array(iList, dtype='float64')


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
    if specfile not in siiContainer.container:
        siiContainer.container[specfile] = list()
    if specfile not in siiContainer.specfiles:
        siiContainer.specfiles.append(specfile)

    if psmType == 'percolator':
        _psmArrays = _importPercolatorResults(fileLocation, psmEngine=psmEngine)
        _importFromPercolatorArrays(siiContainer, _psmArrays, specfile, qValueCutOff=qValue)


def _importFromPercolatorArrays(siiContainer, psmArrays, specfile, qValueCutOff=None):
    """Writes :class:`SpectruIdentificationItem` into :class:`SiContainer`.

    :ivar siiContainer: Add PSMs to this instance
    :ivar psmArrays: contains PSM information, generated by :func:`_importFromPercolatorArray`
    :ivar specfile: Keyword (filename) to represent file in the :class:`SiContainer`
    :ivar qValueCutOff: define a qValue cut off for valid items

    See also :func:`importPsmResults`
    """
    sortMask = psmArrays['score'].argsort()[::-1]
    for key in psmArrays:
        psmArrays[key] = psmArrays[key][sortMask]

    for currPosition in range(0, len(psmArrays['scanNr'])):
        peptide = psmArrays['peptide'][currPosition]
        if peptide.find('.') != -1:
            peptide = peptide.split('.')[1]

        sequence = removeModifications(peptide)
        scanNr = psmArrays['scanNr'][currPosition]
        qValue = psmArrays['q-value'][currPosition]
        score = psmArrays['score'][currPosition]
        psmId = psmArrays['PSMId'][currPosition]
        pep = psmArrays['posterior_error_prob'][currPosition]

        sii = SpectrumIdentificationItem(scanNr, specfile)
        sii.peptide = peptide
        sii.sequence = sequence
        sii.qValue = qValue
        sii.score = score
        sii.pep = pep
        sii.isValid = False

        if sii.containerId in siiContainer.index:
            siiList = siiContainer.index[sii.containerId]
            sii.rank = len(siiList) + 1
        else:
            sii.rank = 1
            siiContainer.index[sii.containerId] = list()

        if sii.rank == 1:
            if qValueCutOff is not None:
                if sii.qValue <= qValueCutOff:
                    sii.isValid = True
            else:
                sii.isValid = True

        siiContainer.index[sii.containerId].append(sii)
        siiContainer.container[specfile].append(sii)


def _importPercolatorResults(fileLocation, psmEngine=None):
    """Reads percolator PSM results from txt file.

    :ivar fileLocation: File path
    :ivar psmEngine: Specifies the used peptide spectrum matching engine ('comet', 'msgf', 'xtandem')

    :return: {attribute:numpy.array(), attribute:numpy.array()}

    See also :func:`importPsmResults` and :func:`_importFromPercolatorArray`
    """
    #HEADERLINE: xtandem seperates proteins with ';', msgf separates proteins by a tab
    with open(fileLocation,'rb') as openFile:
        tsvreader = csv.reader(openFile, delimiter='\t')
        headerLine = tsvreader.next()
        headerDict = dict([[y,x] for (x,y) in enumerate(headerLine)])
        scanEntryList = list()
        for line in tsvreader:
            entryDict = dict()
            for headerName,headerPos in headerDict.items():
                entryDict[headerName] = line[headerPos]
            if psmEngine == 'msgf':
                entryDict['proteinIds'] = list(line[headerDict['proteinIds']:])
            elif psmEngine == 'xtandem':
                entryDict['proteinIds'] = entryDict['proteinIds'].split(';')
            scanEntryList.append(entryDict)

    scanArrDict = dict()
    for headerName in headerDict.keys():
        scanArrDict[headerName] = list()

    # Define list of headers #
    for scanEntryDict in scanEntryList:
        for headerName,entry in scanEntryDict.items():
            if headerName in ['score','q-value','posterior_error_prob']:
                scanArrDict[headerName].append( float(entry) )
            else:
                scanArrDict[headerName].append( entry )

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

    for headerName in scanArrDict.keys():
        scanArrDict[headerName] = numpy.array( scanArrDict[headerName] )
    return scanArrDict


def importPeptideFeatures(featureContainer, filelocation, specfile):
    """ Import peptide features from a featureXml file (eg. generated by OPENMS featureFinderCentroided).

    :param featureContainer: Spectra are added to to this instance of :class:`FeatureContainer`
    :param filelocation: Actual file path
    :param specfile: Keyword (filename) to represent file in the :class:`FeatureContainer`. Each filename
    can only occure once, therefore importing the same filename again is prevented.
    """
    if not os.path.isfile(filelocation):
        print('File does not exits:', filelocation)
    elif not filelocation.lower().endswith('.featurexml'):
        print('File is not a featurexml file:', filelocation)
    else:
        if specfile not in featureContainer.specfiles:
            featureContainer.specfiles.append(specfile)
            featureContainer.container[specfile] = list()
            featureDict = _importFeatureXml(filelocation)

            for featureId, featureEntryDict in featureDict.items():
                rtArea = set()
                for convexHullEntry in featureEntryDict['convexHullDict']['0']:
                    rtArea.update([convexHullEntry[0]])

                featureItem = FeatureItem(featureId, specfile)
                featureItem.rt = featureEntryDict['rt']
                featureItem.rtArea = max(rtArea) - min(rtArea)
                featureItem.rtLow = min(rtArea)
                featureItem.rtHigh = max(rtArea)
                featureItem.charge = featureEntryDict['charge']
                featureItem.mz = featureEntryDict['mz']
                featureItem.mh = aux.calcMhFromMz(featureEntryDict['mz'], featureEntryDict['charge'])
                featureItem.intensity = featureEntryDict['intensity']
                featureItem.quality = featureEntryDict['overallquality']
                featureItem.isMatched = False
                featureItem.isAnnotated = False
                featureItem.isValid = True

                featureContainer.index[featureItem.containerId] = featureItem
                featureContainer.container[specfile].append(featureItem)
        else:
            print(specfile, 'is already present in the SiContainer, import interrupted.')


def _importFeatureXml(fileLocation):
    """Reads a featureXml file.

    :return: {featureKey1: {attribute1:value1, attribute2:value2, ...}, ...}

    See also :func:`importPeptideFeatures`
    """
    with open(fileLocation, 'r') as openFile:
        readingFeature = False
        readingHull = False
        featureDict = dict()

        for i,line in enumerate(openFile):
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


def importProteinDatabase(filePath, proteindb=None, minLength=5, maxLength=40, missedCleavage=2,
                          removeNtermM=True, ignoreIsoleucine=False, fastaType='sgd'
                          ):
    """Generates a :class:`ProteinContainer` and :class:`PeptideContainer` by in silico digestion of proteins from a fasta file.

    :param filePath: File path
    :param ignoreIsoleucine: If True, treat I and L in peptide sequence as indistinguishable
    :param missedCleavages: number of allowed missed digestion sites
    :param removeNtermM: If True, consider peptides with the n-terminal methionine of the protein removed
    :param minLength: only yield peptides with length >= minLength
    :param maxLength: only yield peptides with length <= maxLength
    :param fastaType: see :func:`_importFasta`

    See also :func:`digestInSilico`
    """
    fastaRead = _importFasta(filePath, fastaType=fastaType)
    proteindb = ProteinDatabase() if proteindb is None else proteindb

    for fastaEntry in fastaRead:
        proteinName = fastaEntry['stdName'] if 'stdName' in fastaEntry else fastaEntry['sysName']
        if fastaEntry['sysName'] not in proteindb.proteins:
            protein = ProteinSequence(fastaEntry['sysName'], fastaEntry['sequence'], proteinName)
            proteindb.proteins[fastaEntry['sysName']] = protein
            proteindb.proteinNames[proteinName] = protein

        for unmodPeptide, info in digestInSilico(fastaEntry['sequence'], missedCleavage,
                                                  removeNtermM=True, minLength=minLength,
                                                  maxLength=maxLength
                                                  ):
            if ignoreIsoleucine:
                unmodPeptideNoIsoleucine = unmodPeptide.replace('I', 'L')
                if unmodPeptideNoIsoleucine in proteindb.peptides:
                    currPeptide = proteindb.peptides[unmodPeptideNoIsoleucine]
                else:
                    currPeptide = PeptideSequence(unmodPeptideNoIsoleucine, mc=info['missedCleavage'])
                    proteindb.peptides[unmodPeptideNoIsoleucine] = currPeptide

                if unmodPeptide not in proteindb.peptides:
                    proteindb.peptides[unmodPeptide] = currPeptide
            else:
                if unmodPeptide in proteindb.peptides:
                    currPeptide = proteindb.peptides[unmodPeptide]
                else:
                    currPeptide = PeptideSequence(unmodPeptide, mc=info['missedCleavage'])
                    proteindb.peptides[unmodPeptide] = currPeptide

            if fastaEntry['sysName'] not in currPeptide.proteins:
                currPeptide.proteins.add(fastaEntry['sysName'])
                currPeptide.proteinPositions[fastaEntry['sysName']] = (info['startPos'], info['endPos'])

    for peptide, peptideEntry in proteindb.peptides.items():
        numProteinMatches = len(peptideEntry.proteins)
        if numProteinMatches == 1:
            peptideEntry.isUnique = True
        elif numProteinMatches > 1:
            peptideEntry.isUnique = False
        else:
            print('No protein matches in proteindb for peptide sequence: ', peptide)

        for proteinId in peptideEntry.proteins:
            if peptideEntry.isUnique:
                proteindb.proteins[proteinId].uniquePeptides.add(peptide)
            else:
                proteindb.proteins[proteinId].sharedPeptides.add(peptide)

    for proteinEntry in proteindb.proteins.values():
        if len(proteinEntry.uniquePeptides) > 0:
            proteinEntry.isUnique = True
        else:
            proteinEntry.isUnique = False
    return proteindb


def _importFasta(fastaFileLocation, fastaType='sgd'):
    """Imports fasta files. (Could be merged with or substituted by :class:`pyteomics.fasta.read`)

    :param fastaType: Used to choose which regular expression pattern should be used to read the fasta header line
    :type fastaType: 'sgd' or 'contaminations' or 'uniprot' or 'kustnerPeptideLibrary'

    See also :func:`returnDigestedFasta` and :func:`digestInSilico`
    """
    if fastaType == 'sgd':
        geneAccessionPattern = ">(?P<sysName>[\S]+)\s(?P<stdName>[\S]+).+(?P<description>\".+\")\n(?P<sequence>[A-Z\n]+\*)"
        outputGroups = ['sysName', 'stdName', 'description', 'sequence']
    elif fastaType == 'contaminations':
        geneAccessionPattern = ">(?P<sysName>[\S]+)\s(?P<description>.+)\n(?P<sequence>[A-Z\n]+\*)"
        outputGroups = ['sysName', 'description', 'sequence']
    elif fastaType == 'kustnerPeptideLibrary':
        geneAccessionPattern = ">(?P<sysName>IPI:[^\|\s]+)(\n|(.+\n))(?P<sequence>[A-Z\n]+)"
        #sysName = "match many [non white space characters, not pipe]
        #description = "match many [non newline characters]" (?P<description>.+)\n
        #sequence = "match many[Letters or newLine]
        outputGroups = ['sysName', 'sequence']
    elif fastaType == 'ipi':
        geneAccessionPattern = ">(?P<sysName>[\S]+)(.+\n)(?P<sequence>[A-Z\n]+)"
        outputGroups = ['sysName', 'sequence']
    elif fastaType == 'uniprot':
        geneAccessionPattern = ">[\S]+\|(?P<sysName>[\S]+)\|(?P<stdName>[\S^\|]+)\s(?P<description>.+)\n(?P<sequence>[A-Z\n]+)"
        outputGroups = ['sysName', 'stdName', 'description', 'sequence']

    with open(fastaFileLocation, 'r') as openFastaFile:
        readFastaFile = openFastaFile.read()
        proteinList = list()

        regexpPattern = re.compile(geneAccessionPattern, re.VERBOSE)
        regexpResult = regexpPattern.finditer(str(readFastaFile))
        for entry in regexpResult:
            outputDict = dict()
            outputDict['description'] = ''
            for outputGroup in outputGroups:
                outputDict[outputGroup] = str(entry.group(outputGroup))
            outputDict['sequence'] = outputDict['sequence'].replace('\n', '').replace('*', '')
            proteinList.append(outputDict)
    return proteinList


def digestInSilico(proteinSequence, missedCleavages, removeNtermM=True, minLength=5, maxLength=40):
    """Yields peptides derived from an in silico digest of a protein.

    (Could be merged with or substituted by :class:`pyteomics.fasta.read`)

    :param proteinSequence: amino acid sequence of the protein to be digested
    :param missedCleavages: number of allowed missed cleavage sites
    :param removeNtermM: If True, consider peptides with the n-terminal methionine of the protein removed
    :param minLength: only yield peptides with length >= minLength
    :param maxLength: only yield peptides with length <= maxLength

    NOTE: at the moment it only works for trypsin (K/R) and c-terminal cleavage
    """
    # Return in silico digested peptides, peptide start position, peptide end position
    # Peptide position start at 1 and end at len(proteinSequence)
    passFilter = lambda startPos, endPos: (endPos - startPos >= minLength and endPos - startPos <= maxLength)

    cleavagePosList = list()
    # Position +1 if cut c terminal of amino acid
    cleavagePosList.extend([m.start()+1 for m in re.finditer('K', proteinSequence)])
    cleavagePosList.extend([m.start()+1 for m in re.finditer('R', proteinSequence)])
    cleavagePosList.sort()

    # Add end of protein as cleavage site if protein doesn't end with specififed cleavage positions
    if proteinSequence[-1] != 'K' and proteinSequence[-1] != 'R':
        cleavagePosList.append(len(proteinSequence))
    numCleavageSites = len(cleavagePosList)

    if missedCleavages >= numCleavageSites:
        missedCleavages = numCleavageSites -1

    #Yield protein n-terminal peptides
    if cleavagePosList[0] != 0:
        for cleavagePos in range(0,missedCleavages+1):
            startPos = 0
            endPos = cleavagePosList[cleavagePos]
            if passFilter(startPos, endPos):
                sequence = proteinSequence[startPos:endPos]
                info = dict()
                info['startPos'] = startPos+1
                info['endPos'] = endPos
                info['missedCleavage'] = cleavagePos
                yield sequence, info

    #Yield protein n-terminal peptides after methionine removal
    if removeNtermM and proteinSequence[0] == 'M':
        for cleavagePos in range(0,missedCleavages+1):
            startPos = 1
            endPos = cleavagePosList[cleavagePos]
            if passFilter(startPos, endPos):
                sequence = proteinSequence[startPos:endPos]
                info = dict()
                info['startPos'] = startPos+1
                info['endPos'] = endPos
                info['missedCleavage'] = cleavagePos
                yield sequence, info

    #Yield all remaining peptides, including the c-terminal peptides
    lastCleavagePos = 0
    while lastCleavagePos < numCleavageSites:
        for missedCleavage in range(0, missedCleavages+1):
            nextCleavagePos = lastCleavagePos + missedCleavage + 1
            if nextCleavagePos < numCleavageSites:
                startPos = cleavagePosList[lastCleavagePos]
                endPos = cleavagePosList[nextCleavagePos]
                if passFilter(startPos, endPos):
                    sequence = proteinSequence[startPos:endPos]
                    info = dict()
                    info['startPos'] = startPos+1
                    info['endPos'] = endPos
                    info['missedCleavage'] = missedCleavage
                    yield sequence, info
        lastCleavagePos += 1


def _calculateCoverageMasks(proteindb, peptidedb):
    """Calcualte the sequence coverage masks for all proteindb elements.
    Private method used by :class:`ProteinDatabase` and :class:`EvidenceContainer`

    A coverage mask is a numpy boolean array with the length of the protein sequence.
    Each protein position that has been covered in at least one peptide is set to True.
    Coverage masks are calculated for unique and for shared peptides. Peptides are
    matched to proteins according to positions derived by the digestion of the FASTA file.

    Alternatively peptides could also be matched to proteins just by sequence as it is
    done in :func:`pyteomics.parser.coverage`, but this is not the case here.

    :ivar :attr:`PeptideSequence.coverageMaskUnique`: coverage mask of unique peptides
    :ivar :attr:`PeptideSequence.coverageMaskShared`: coverage mask of shared peptides
    """
    for proteinId, proteinEntry in proteindb.items():
        coverageMaskUnique = numpy.zeros(proteinEntry.length, dtype='bool')
        for peptide in proteinEntry.uniquePeptides:
            startPos, endPos = peptidedb[peptide].proteinPositions[proteinId]
            coverageMaskUnique[startPos-1:endPos] = True
        coverageMaskShared = numpy.zeros(proteinEntry.length, dtype='bool')
        for peptide in proteinEntry.sharedPeptides:
            startPos, endPos = peptidedb[peptide].proteinPositions[proteinId]
            coverageMaskShared[startPos-1:endPos] = True
        setattr(proteinEntry, 'coverageMaskUnique', coverageMaskUnique)
        setattr(proteinEntry, 'coverageMaskShared', coverageMaskShared)


def generateEvidence(evContainer, siiContainer, scoreKey, largerBetter, peptideKey='peptide'):
    """Summarize all experimental evidence from :class:`SiiContainer` and generate an :class:`EvidenceContainer`.

    :ivar evContainer: :class:`EvidenceContainer` to store the evidence information
    :ivar siiContainer: :class`SiiContainer` instance that contains the experimental evidence
    :ivar peptideKey: 'sequence' -> ignore modification, uses sequence as unique key in
    peptide evidence; 'peptide' -> consider modified peptides as unique entries
    :ivar scoreKey: score attribtue name of :class:`SiiItem`
    :ivar largerBetter: boolean, True if a larger score is better

    NOTE: Peptides which sequence is not present in evContainer.db.peptides are ignored and generate an error message.
    """
    evContainer.scoreKey = scoreKey
    evContainer.largerBetter = largerBetter
    evContainer.peptideKey = peptideKey

    #Summarize psm informatino into unique peptides
    for sii in siiContainer.getItems(sort=scoreKey, reverse=largerBetter):
        peptide = getattr(sii, peptideKey)
        siiScore = getattr(sii, scoreKey)
        if peptide in evContainer.peptideEvidence:
            evContainer.peptideEvidence[peptide].siiIds.append(sii.containerId)
            evContainer.peptideEvidence[peptide].scores.append(siiScore)
        else:
            try:
                pepEv = PeptideEvidence.fromPeptideSequence(peptide, evContainer.db.peptides[sii.sequence])
                pepEv.bestId = sii.containerId
                pepEv.siiIds.append(sii.containerId)
                pepEv.score = siiScore
                pepEv.scores.append(siiScore)
                evContainer.peptideEvidence[peptide] = pepEv
            except KeyError:
                print('Sequence not found in evContainer.db.peptides: ', sii.sequence, sii.containerId)
                pass

    #Assemble peptide evidence into protein evidence
    for peptide, pepEv in evContainer.peptideEvidence.items():
        for protein in pepEv.proteins:
            if protein in evContainer.proteinEvidence:
                proteinEv = evContainer.proteinEvidence[protein]
            else:
                proteinEv = ProteinEvidence(protein, evContainer.db.proteins[protein].sequence,
                                            evContainer.db.proteins[protein].name
                                            )
                evContainer.proteinEvidence[protein] = proteinEv

            if pepEv.isUnique:
                proteinEv.uniquePeptides.add(peptide)
                proteinEv.uniquePsmCount += len(pepEv.siiIds)
            else:
                proteinEv.sharedPeptides.add(peptide)
                proteinEv.sharedPsmCount += len(pepEv.siiIds)

    #Define proteins which have unique evidence
    evContainer.uniqueProteins = list()
    for proteinEv in evContainer.proteinEvidence.values():
        if len(proteinEv.uniquePeptides) > 0:
            proteinEv.isUnique = True
            evContainer.uniqueProteins.append(proteinEv.id)
        else:
            proteinEv.isUnique = False


################################################
### Functions to work with peptide sequences ###
################################################
def calcPeptideMass(peptide):
    """Calculate the mass of a peptide. (Should be changed to allow for modifications not present in unimod.org)

    :ivar peptide: peptide sequence, modifications have to be written in the format "[modificationName]"
    and 'modificationName' has to be present in aux.unimodToMassDict
    """
    unimodMassDict = aux.unimodToMassDict

    additionalModMass = float()
    unmodPeptide = peptide
    for unimodNumber, unimodMass in unimodMassDict.items():
        try:
            int(unimodNumber)
        except ValueError:
            unimodSymbol = '[' + unimodNumber + ']'
        else:
            unimodSymbol = '[UNIMOD:' + unimodNumber + ']'
        numMod = peptide.count(unimodSymbol)
        unmodPeptide = unmodPeptide.replace(unimodSymbol, '')
        additionalModMass += unimodMass * numMod

    if unmodPeptide.find('[') != -1:
        print(unmodPeptide)
        raise Exception()

    unmodPeptideMass = pyteomics.mass.calculate_mass(unmodPeptide, charge=0)
    modPeptideMass = unmodPeptideMass + additionalModMass
    return modPeptideMass


def removeModifications(peptide):
    """Removes all modifications from a peptide string

    :ivar peptide: peptide sequence, modifications have to be written in the format "[modificationName]"
    :type peptide: str
    """
    while peptide.find('[') != -1:
        peptide = peptide.split('[', 1)[0] + peptide.split(']', 1)[1]
    return peptide


def returnModPositions(peptide, indexStart=1, removeModString='UNIMOD:'):
    """Determines the amino acid positions of all present modifications.

    :ivar peptide: peptide sequence, modifications have to be written in the format "[modificationName]"
    :ivar indexStart: returned amino acids positions of the peptide start with this number (1st amino acid position = indexStart)
    :ivar removeModString: string to remove from the returned modification name

    :return: {modificationName:[position1, position2, ...], ...}
    """
    unidmodPositionDict = dict()
    while peptide.find('[') != -1:
        currModification = peptide.split('[')[1].split(']')[0]
        currPosition = peptide.find('[') - 1
        if currPosition == -1: # move n-terminal modifications to first position
            currPosition = 0
        currPosition += indexStart

        peptide = peptide.replace('['+currModification+']', '', 1)

        if isinstance(removeModString, str):
            currModification = currModification.replace(removeModString, '')
        unidmodPositionDict.setdefault(currModification,list())
        unidmodPositionDict[currModification].append(currPosition)
    return unidmodPositionDict


############################################
## Functions dealing with isotopic labels ##
############################################
class LabelDescriptor(object):
    """Describes a MS1 stable isotope label setup for quantification.

    :ivar labels: Contains a dictionary with all possible label states, keys (=labelStates) are increasing integers starting from 0
    :ivar excludingModifictions: bool, set to True if any label has specified excludingModifications
    """
    def __init__(self):
        self.labels = dict()
        self.excludingModifictions = False
        self._labelCounter = 0

    def addLabel(self, aminoAcidLabels, excludingModifications=None):
        """Adds a new labelstate.

        :ivar aminoAcidsLabels: Describes which amino acids can bear which labels
        possible keys amino acids in one letter code and ('nTerm', 'cTerm')
        possible values are keys from :var:`pyms.auxiliary.unimodToMassDict` as strings or list of strings
        eg. {'nTerm':'188', 'K':['188', '188']} for one expected label at the nterminus and two expected labels at Lysine
        :ivar excludingModifications: Describes which modifications can prevent the addition of labels
        keys and values have to be keys from :var:`pyms.auxiliary.unimodToMassDict` written as a string.
        eg. {'1':'188'} For each modification '1' that is present at an amino acid or terminus of a peptide
        the number of expected labels at this position is reduced by one
        """
        if excludingModifications is not None:
            self.excludingModifictions = True

        self.labels[self._labelCounter] = dict()
        self.labels[self._labelCounter]['aminoAcidLabels'] = aminoAcidLabels
        self.labels[self._labelCounter]['excludingModifications'] = excludingModifications
        self._labelCounter += 1


def returnLabelStateMassDifferences(peptide, labelDescriptor, labelState=None, sequence=None):
    """Calculates the mass difference for alternative possible label states of a given peptide.

    :ivar peptide: Peptide to calculate alternative label states
    :ivar labelDescriptor: :class:`LabelDescriptor` describes the label setup of an experiment
    :ivar labelState: label state of the peptide, if None it is calculated by :func:`returnLabelState`
    :ivar sequence: unmodified amino acid sequence of :var:`peptide`, if None it is calculated with :func:`removeModifications`

    :return: {alternativeLabelSate: massDifference, ...} or {} if the peptide label state is -1
    (massDifference + peptide mass = expected mass of alternatively labeled peptide)

    See also :class:`LabelDescriptor`, :func:`returnLabelState`
    """
    labelState = returnLabelState(peptide, labelDescriptor) if labelState is None else labelState
    sequence = removeModifications(peptide) if sequence is None else sequence

    if labelState < 0:
        # special case for mixed label... #
        return dict()

    # define type and number of labels of the peptide
    labelModNumbers = dict()
    for labelStateModList in expectedLabelPosition(peptide, labelDescriptor.labels[labelState], sequence=sequence).values():
        for labelMod in labelStateModList:
            labelModNumbers.setdefault(labelMod, int())
            labelModNumbers[labelMod] += 1

    # calculate the combined labels mass of the peptide
    labelMass = int()
    for labelMod, modCounts in labelModNumbers.items():
        labelMass += aux.unimodToMassDict[labelMod] * modCounts

    # calculate mass differences to all other possible label states
    labelStateMassDifferences = dict()
    for possibleLabelState in labelDescriptor.labels.keys():
        if possibleLabelState == labelState:
            continue

        labelModNumbers = dict()
        for labelStateModList in expectedLabelPosition(peptide, labelDescriptor.labels[possibleLabelState], sequence=sequence).values():
            for labelMod in labelStateModList:
                labelModNumbers.setdefault(labelMod, int())
                labelModNumbers[labelMod] += 1

        possibleLabelMass = int()
        for labelMod, modCounts in labelModNumbers.items():
            possibleLabelMass += aux.unimodToMassDict[labelMod] * modCounts

        possibleLabelMassDifference = possibleLabelMass - labelMass
        labelStateMassDifferences[possibleLabelState] = possibleLabelMassDifference
    return labelStateMassDifferences


def returnLabelState(peptide, labelDescriptor, labelSymbols=None, labelAminoacids=None):
    """Calculates the label state of a given peptide for the label setup described in labelDescriptor

    :ivar peptide: peptide which label state should be calcualted
    :ivar labelDescriptor: :class:`LabelDescriptor` describes the label setup of an experiment
    :ivar labelSymbols: modifications that show a label, calculated by :func:`modSymbolsFromLabelInfo`
    :ivar labelAminoacids: amino acids that can bear a label, calculated by :func:`modAminoacidsFromLabelInfo`

    :return: integer that shows the label state
    >=0: predicted label state of the peptide
     -1: peptide sequence can't bear any labelState modifications
     -2: peptide modifications don't fit to any predicted labelState
     -3: peptide modifications fit to a predicted labelState, but not all predicted labelStates are distinguishable
    """
    labelSymbols = modSymbolsFromLabelInfo(labelDescriptor) if labelSymbols is None else labelSymbols
    labelAminoacids = modAminoacidsFromLabelInfo(labelDescriptor) if labelAminoacids is None else labelAminoacids

    sequence = removeModifications(peptide)
    modPositions = returnModPositions(peptide, indexStart=0)

    labelState = None

    # No amino acids in sequence which can bear a label modification (ignores presence of excluding modifications)
    if all([(True if sequence.find(labelAminoacid) == -1 else False) for labelAminoacid in labelAminoacids]):
        # No terminal label modifications specified by labelDescriptor
        if 'nTerm' not in labelAminoacids and 'cTerm' not in labelAminoacids:
            labelState = -1

    # Check if the peptide mofidifcations fit to any predicted label state
    if labelState is None:
        peptideLabelPositions = dict()
        for labelSymbol in labelSymbols:
            if labelSymbol in modPositions.keys():
                for sequencePosition in modPositions[labelSymbol]:
                    peptideLabelPositions.setdefault(sequencePosition, list())
                    peptideLabelPositions[sequencePosition].append(labelSymbol)
        for sequencePosition in peptideLabelPositions.keys():
            peptideLabelPositions[sequencePosition] = sorted(peptideLabelPositions[sequencePosition])

        predictedLabelStates = dict()
        for predictedLabelState, labelStateInfo in labelDescriptor.labels.items():
            expectedLabelMods = expectedLabelPosition(peptide, labelStateInfo, sequence=sequence, modPositions=modPositions)
            predictedLabelStates[predictedLabelState] = expectedLabelMods
            if peptideLabelPositions == expectedLabelMods:
                # If another expectedLabel state has already been matched, there is an ambiguity between label states...
                labelState = predictedLabelState

    if labelState is None:
        # Peptide mofidifcations don't fit to any predicted label state
        labelState = -2
    elif labelState != -1:
        # Check if all predicted label states are distinguishable
        for labelState1, labelState2 in set(itertools.combinations(range(len(predictedLabelStates)), 2)):
            if predictedLabelStates[labelState1] == predictedLabelStates[labelState2]:
                labelState = -3
                break

    return labelState


def modSymbolsFromLabelInfo(labelDescriptor):
    """Returns a set of all modiciation symbols which were used in the labelDescriptor

    :ivar labelDescriptor: :class:`LabelDescriptor` describes the label setup of an experiment
    """
    modSymbols = set()
    for labelStateEntry in labelDescriptor.labels.values():
        for labelPositionEntry in labelStateEntry['aminoAcidLabels'].values():
            for modSymbol in aux.toList(labelPositionEntry):
                if modSymbol != '':
                    modSymbols.add(modSymbol)
    return modSymbols


def modAminoacidsFromLabelInfo(labelDescriptor):
    """Returns a set of all amino acids and termini which can bear a label, as described in labelDescriptor

    :ivar labelDescriptor: :class:`LabelDescriptor` describes the label setup of an experiment
    """
    modAminoacids = set()
    for labelStateEntry in labelDescriptor.labels.values():
        for labelPositionEntry in labelStateEntry['aminoAcidLabels'].keys():
            for modAminoacid in aux.toList(labelPositionEntry):
                if modAminoacid != '':
                    modAminoacids.add(modAminoacid)
    return modAminoacids


def expectedLabelPosition(peptide, labelStateInfo, sequence=None, modPositions=None):
    """Returns a modification description of a certain label state of a peptide.

    :ivar peptide: Peptide sequence used to calculat the expected label state modifications
    :ivar labelStateInfo: An entry of :attr:`LabelDescriptor.labels` that describes a label state
    :ivar sequence: unmodified amino acid sequence of :var:`peptide`, if None it is calculated with :func:`removeModifications`
    :ivar modPositions: dictionary describing the modification state of :var:`peptide`,
    if None it is calculated with :func:`returnModPositions`

    :return: {sequence position: sorted list of expected label modifications on that position, ...}
    """
    modPositions = returnModPositions(peptide, indexStart=0) if modPositions is None else modPositions
    sequence = removeModifications(peptide) if sequence is None else sequence

    currLabelMods = dict()
    for labelPosition, labelSymbols in labelStateInfo['aminoAcidLabels'].items():
        labelSymbols = aux.toList(labelSymbols)
        if labelSymbols == ['']:
            pass
        elif labelPosition == 'nTerm':
            currLabelMods.setdefault(0, list())
            currLabelMods[0].extend(labelSymbols)
        else:
            for sequencePosition in aux.findAllSubstrings(sequence, labelPosition):
                currLabelMods.setdefault(sequencePosition, list())
                currLabelMods[sequencePosition].extend(labelSymbols)

    if labelStateInfo['excludingModifications'] is not None:
        for excludingModification, excludedLabelSymbol in labelStateInfo['excludingModifications'].items():
            if excludingModification in modPositions:
                for excludingModPosition in modPositions[excludingModification]:
                    if excludingModPosition in currLabelMods:
                        if excludedLabelSymbol in currLabelMods[excludingModPosition]:
                            if len(currLabelMods[excludingModPosition]) == 1:
                                del(currLabelMods[excludingModPosition])
                            else:
                                excludedModIndex = currLabelMods[excludingModPosition].index(excludedLabelSymbol)
                                currLabelMods[excludingModPosition].pop(excludedModIndex)

    for sequencePosition in currLabelMods.keys():
        currLabelMods[sequencePosition] = sorted(currLabelMods[sequencePosition])
    return currLabelMods


##########################################
### Functions to work with FeatureItem ###
##########################################
def matchToFeatures(featureContainer, specContainer, specfiles=None, fMassKey='mz', sMassKey='obsMz', isotopeErrorList=(0, 1),
                    precursorTolerance=5, toleranceUnit='ppm', rtExpansionUp=0.10, rtExpansionDown=0.05, matchCharge=True,
                    scoreKey='pep', largerBetter=False):
    """Annotate :class:`FeatureItem` in :class:`FeatureContainer` by matching :class:`SpectrumItem` (Si) or :class:`SpectrumIdentificationItem` (Sii).

    :ivar featureContainer: :class:`FeatureItem` of this instance of :class:`FeatureContainer` are annotated
    :type specContainer: :class:`SiContainer` or :class:`SiiContainer`
    :ivar specfiles: Annotate only items of :attr:`FeatureContainer.container[specfile]`,
    if None all specfiles present in featureContainer and specContainer are processed
    :type specfiles: str, list or None
    :ivar fMassKey: mass attribute key in :attr:`FeatureItem.__dict__`
    :ivar sMassKey: mass attribute key in :attr:`SpectrumItem.__dict__` or :attr:`SpectrumIdentificationItem.__dict__` (eg 'obsMz', 'calcMz')
    :ivar isotopeErrorList: allowed isotope errors relative to the spectrum mass
    eg. [0, 1] if no feature has been matched with isotope error 0, the spectrum mass is increased by 1*C13 and matched again
    the different isotope error values are tested in the specified order therefore 0 should normally be the 1st value of the tuple
    :type isotopeErrorList: list or tuple of int
    :ivar precursorTolerance: is used to calculate the mass window to match Si or Sii to :class:`FeatureItem`
    :ivar toleranceUnit: defines how the precursorTolerance is applied to the mass value, 'ppm' * (1 +/- tolerance*1E-6) or 'da': +/- value
    :ivar rtExpansionUp: relative upper expansion of :class:`FeatureItem` retention time areas
    :ivar rtExpansionDown: relative lower expansion of :class:`FeatureItem` retention time areas
    Ad rtExpansion; lower or upper :class:`FeatureItem` retention time position is expanded with rtExpansion * rtArea
    if Si or Sii is matched to multiple :class:`FeatureItem`, the rt expansion is removed and the matching repeated.
    :ivar matchCharge: boolean, True if :class:`FeatureItem` and Si/Sii must have the same charge state to be matched
    :ivar scoreKey: Sii attribute name of the PSM score
    :ivar largerBetter: boolean, True if a larger PSM score is better

    If specContainer is :class:`SiiContainer` then matched features are annotated with :attr:`SpectrumIdentificationItem.peptide`,
    if multiple Sii are matched to the same :class:`FeatureItem` the one with the best score is used
    """
    isotopeErrorList = aux.toList(isotopeErrorList)

    if specContainer.__class__.__name__ == 'SiiContainer':
        listKeySpecIds = 'siiIds'
    else:
        listKeySpecIds = 'siIds'

    if specfiles is not None:
        specfiles = aux.toList(specfiles)
    else:
        specfiles = list(set(featureContainer.specfiles).intersection(set(specContainer.specfiles)))

    for specfile in specfiles:
        multiMatchCounter = int()
        isotopeErrorMatchCounter = int()
        specArrays = specContainer.getArrays([sMassKey, 'rt', 'charge'], specfiles=specfile)
        featureArrays = featureContainer.getArrays(['rtHigh', 'rtLow', 'charge', fMassKey],
                                                   specfiles=specfile, sort=fMassKey
                                                   )
        featureArrays['rtHighExpanded'] = featureArrays['rtHigh'] + (featureArrays['rtHigh'] - featureArrays['rtLow']) * rtExpansionUp
        featureArrays['rtLowExpanded'] = featureArrays['rtLow'] - (featureArrays['rtHigh'] - featureArrays['rtLow']) * rtExpansionDown

        specFeatureDict = dict() ## key = scanNr, value = set(featureKeys)
        featureSpecDict = dict() ## key = featureKey, value = set(scanNrs)

        for specPos, specId in enumerate(specArrays['containerId']):
            specMass = specArrays[sMassKey][specPos]
            specRt = specArrays['rt'][specPos]
            specZ = specArrays['charge'][specPos]
            if specZ is None:
                continue

            matchComplete = False
            isotopeErrorPos = 0

            while not matchComplete:
                isotopeError = isotopeErrorList[isotopeErrorPos]

                # calculate mass limits for each isotope error
                if toleranceUnit.lower() == 'ppm':
                    specMassHigh = (specMass + isotopeError * 1.003355 / specZ) * (1 + precursorTolerance*1E-6)
                    specMassLow = (specMass + isotopeError * 1.003355 / specZ) * (1 - precursorTolerance*1E-6)
                elif toleranceUnit.lower() == 'da':
                    specMassHigh = (specMass + isotopeError * 1.003355 / specZ) + precursorTolerance
                    specMassLow  = (specMass + isotopeError * 1.003355 / specZ) - precursorTolerance

                posL = bisect.bisect_left(featureArrays[fMassKey], specMassLow)
                posR = bisect.bisect_right(featureArrays[fMassKey], specMassHigh)

                if matchCharge:
                    chargeMask = (featureArrays['charge'][posL:posR] == specZ)

                fRtHighKey = 'rtHighExpanded'
                fRtLowKey = 'rtLowExpanded'
                for fRtHighKey, fRtLowKey in [('rtHighExpanded', 'rtLowExpanded'), ('rtHigh', 'rtLow')]:
                    rtMask = ((featureArrays[fRtLowKey][posL:posR] <= specRt) &
                              (featureArrays[fRtHighKey][posL:posR] >= specRt)
                              )
                    if matchCharge:
                        matchedFeatureIds = featureArrays['containerId'][posL:posR][rtMask & chargeMask]
                    else:
                        matchedFeatureIds = featureArrays['containerId'][posL:posR][rtMask]

                    if len(matchedFeatureIds) <= 1:
                        break

                # if exactly one feature has been matched,
                if len(matchedFeatureIds) > 0:
                    if len(matchedFeatureIds) == 1:
                        matchComplete = True
                        if isotopeErrorList[isotopeErrorPos] != 0:
                            isotopeErrorMatchCounter += 1
                    else:
                        #Stop if Spectrum can be matched to multiple features
                        multiMatchCounter += 1
                        break

                isotopeErrorPos += 1
                if isotopeErrorPos >= len(isotopeErrorList):
                    #Stop if all allowed isotope errors have been tested
                    break

            if matchComplete:
                for featureId in matchedFeatureIds:
                    getattr(featureContainer[featureId], listKeySpecIds).append(specId)
                    featureContainer[featureId].isMatched = True
                    specFeatureDict[specId] = featureId
                    featureSpecDict[featureId] = specId

        stats = dict()
        stats['totalFeatures'] = len(featureArrays['containerId'])
        stats['matchedFeatures'] = len(featureSpecDict)
        stats['relMatchedFeatures'] = round(1.*stats['matchedFeatures']/stats['totalFeatures'], 3)
        stats['totalSpectra'] = len(specArrays['containerId'])
        stats['matchedSpectra'] = len(specFeatureDict)
        stats['relMatchedSpectra'] = round(1.*stats['matchedSpectra']/stats['totalSpectra'], 3)

        print('------', specfile, '------')
        print('Annotated features:\t\t\t', stats['matchedFeatures'], '/', stats['totalFeatures'], '=', stats['relMatchedFeatures'], '%')
        print('Spectra matched to features:\t\t', stats['matchedSpectra'], '/', stats['totalSpectra'], '=', stats['relMatchedSpectra'], '%')
        if multiMatchCounter != 0:
                print('Discarded because of multiple matches:\t', multiMatchCounter)
        if isotopeErrorMatchCounter != 0:
                print('Isotope error matched spectra:\t\t', isotopeErrorMatchCounter)

        if specContainer.__class__.__name__ == 'SiiContainer':
            for featureId in featureSpecDict.keys():
                matches = list()
                for specId in featureContainer[featureId].siiIds:
                    _sii = specContainer.getValidItem(specId)
                    score = getattr(_sii, scoreKey)
                    peptide = _sii.peptide
                    sequence = _sii.sequence
                    matches.append([score, peptide, sequence])
                matches.sort(reverse=largerBetter)

                featureContainer[featureId].isAnnotated = True
                featureContainer[featureId].score = matches[0][0]
                featureContainer[featureId].peptide = matches[0][1]
                featureContainer[featureId].sequence = matches[0][2]


def removeFeatureAnnotation(featureContainer):
    """Remove all annotation information from :class:`FeatureItem` in :class:`featureContainer`."""
    for items in featureContainer.container.values():
        for item in items:
            item.isMatched = False
            item.isAnnotated = False
            item.siIds = list()
            item.siiIds = list()
            item.peptide = None
            item.sequence = None
            item.score = None


def rtCalibration(featureContainer, allowedRtDev=60, allowedMzDev=2.5, showPlots=False, reference=None, specfiles=None):
    """Performs a retention time calibration between :class:`FeatureItem` of multiple specfiles.

    :ivar featureContainer: Perform alignment on :class:`FeatureItem` in :attr:`FeatureContainer.specfiles`
    :ivar allowedRtDev: maxium retention time difference of two features in two runs to be matched
    :ivar allowedMzDev: maxium relative m/z difference (in ppm) of two features in two runs to be matched
    :ivar showPlots: boolean, True if a plot should be generated which shows to results of the calibration
    :ivar reference: Can be used to specifically specify a reference specfile
    :ivar specfiles: Limit alignment to those specfiles in the featureContainer
    """
    specfiles = featureContainer.specfiles if specfiles is None else specfiles
    matchCharge = True

    refMzKey = 'mz'
    mzKey = 'mz'
    minIntensity = 1E5

    if reference is not None:
        if reference in specfiles:
            specfiles = [reference] + list(set(specfiles).difference(set([reference])))
        else:
            print('Specified reference specfile not present, using reference: ', specfiles[0])

    for featureItem in featureContainer.getItems(specfiles=specfiles, filterAttribute=None):
        if not hasattr(featureItem, 'obsRt'):
            setattr(featureItem, 'obsRt', featureItem.rt)

    referenceArrays = None
    for specfile in specfiles:
        featureArrays = featureContainer.getArrays(['rt', 'charge', 'mz', 'intensity'],
                                                   specfiles=specfile, sort='rt'
                                                   )
        if minIntensity is not None:
            intensityMask = (featureArrays['intensity'] > minIntensity)
            for key in featureArrays.keys():
                featureArrays[key] = featureArrays[key][intensityMask]

        if referenceArrays is None:
            referenceArrays = featureArrays
            if showPlots:
                print('Reference: '+specfile)
            continue

        rtPosList = list()
        rtDevList = list()
        mzDevRelList = list()
        mzDevAbsList = list()

        for featurePos in xrange(len(featureArrays[mzKey])):
            currRt = featureArrays['rt'][featurePos]
            currMz = featureArrays[mzKey][featurePos]
            currZ = featureArrays['charge'][featurePos]
            mzLimitUp = currMz*(1+allowedMzDev*1E-6)
            mzLimitLow = currMz*(1-allowedMzDev*1E-6)
            rtLimitUp = currRt+allowedRtDev
            rtLimitLow = currRt-allowedRtDev

            posL = bisect.bisect_left(referenceArrays['rt'], rtLimitLow)
            posU = bisect.bisect_right(referenceArrays['rt'], rtLimitUp)

            refMask = (referenceArrays[refMzKey][posL:posU] <= mzLimitUp) & (referenceArrays[refMzKey][posL:posU] >= mzLimitLow)
            if matchCharge:
                refMask = refMask & (referenceArrays['charge'][posL:posU] == currZ)

            currMzDev = abs(referenceArrays[refMzKey][posL:posU][refMask] - currMz)
            bestHitMask = currMzDev.argsort()
            for refRt, refMz in zip(referenceArrays['rt'][posL:posU][refMask][bestHitMask],
                                    referenceArrays[refMzKey][posL:posU][refMask][bestHitMask]):
                rtPosList.append(currRt)
                rtDevList.append(currRt - refRt)
                mzDevRelList.append((1 - currMz / refMz)*1E6)
                mzDevAbsList.append(currMz - refMz)
                break

        rtPosList = numpy.array(rtPosList)
        rtDevList = numpy.array(rtDevList)

        splineInitialKnots = int(max(rtPosList)-min(rtPosList))
        dataFit = aux.DataFit(rtDevList, rtPosList)
        dataFit.splineInitialKnots=splineInitialKnots
        dataFit.splineTerminalExpansion=0.2
        dataFit.processInput(dataAveraging='median', windowSize=10)
        dataFit.generateSplines()

        if showPlots:
            corrDevArr = rtDevList - dataFit.corrArray(rtPosList)
            timePoints = [min(rtPosList) + x for x in range(int(max(rtPosList)-min(rtPosList)))]
            corrValues  = dataFit.corrArray(timePoints)
            plt.supfig, ( ax ) = plt.subplots(3, 2, sharex=False, sharey=False, figsize=(20, 18))
            plt.suptitle(specfile)
            ax[0][0].hist(rtDevList, bins = 100, color='grey', alpha=0.5, label='observed')
            ax[0][0].hist(corrDevArr, bins = 100, color='red', alpha=0.5, label='corrected')
            ax[0][0].set_title('Retention time deviation')
            ax[0][0].legend()
            ax[0][1].hist(mzDevRelList, bins = 100, color = 'grey')
            ax[0][1].set_title('Mz deviation [ppm]')
            ax[1][0].scatter(rtPosList, rtDevList, color = 'grey', alpha = 0.1, label='observed')
            ax[1][0].plot(timePoints,corrValues, color='red', alpha=0.5, label='correction function')
            ax[1][0].set_title('Retention time deviation over time')
            ax[1][0].legend()
            ax[1][1].scatter(rtPosList, mzDevRelList, color = 'grey', alpha = 0.1)
            ax[1][1].set_title('Mz deviation over time')
            ax[2][0].scatter(rtPosList, corrDevArr, color = 'grey', alpha = 0.1)
            ax[2][0].set_title('Aligned retention time deviation over time')
            plt.show()

        featureArrays = featureContainer.getArrays(['rt'], specfiles=specfile, sort='rt', filterAttribute=None)
        featureArrays['corrRt'] = featureArrays['rt'] - dataFit.corrArray(featureArrays['rt'])
        for featureId, corrRt, rt in zip(featureArrays['containerId'], featureArrays['corrRt'], featureArrays['rt']):
            featureContainer.index[featureId].rt = corrRt


def groupFeatures(featureContainer, specfiles=None, featureFilterAttribute='isAnnotated', massTolerance=3,
                  toleranceUnit='ppm', rtTolerance=10, massKey='mh', rtKey='rt', largerBetter=False,
                  matchCharge=False, labelDescriptor=LabelDescriptor(), targetChargeStates=(2, 3, 4, 5, 6)
                  ):
    """Function to group :class:`FeatureItem` according to mass, retention time, charge and label state.

    Feature groups can be used for labelfree quantification and quantification by labeling with stable isotopes.
    If multiple specfiles are specified, the :class:`FeatureItem` annotation information is shared between the
    specfiles. Feature groups are only annotated with the best scoring identification of all grouped :class:`FeatureItem`.
    Before grouping :class:`FeatureItem` of multiple specfiles, a retention time alignment should be performed with
    e.g. :func:`rtCalibration`.

    :ivar featureContainer: :class:`FeatureContainer` contains :class:`FeatureItem` which should be grouped
    :ivar specfiles: Group only items of :attr:`FeatureContainer.container[specfile]`,
    if None all specfiles present in featureContainer are processed
    :type specfiles: str, list or None
    :ivar featureFilterAttribute: Attribute for :param:`FeatureContainer.getItems(filterAttribute=featureFilterAttribute)`,
    used to filter :class:`FeatureItem` elements that should be used as seeds for feature grouping ('isValid', 'isMatched', 'isAnnotated')
    'isValid' all valid features are used for grouping, 'isMatched' all features with matched MSn items are used for grouping or
    'isAnnotated' to use only features with an peptide id for grouping.
    :ivar massTolerance: is used to calculate the mass window to match :class:`FeatureItem`
    :ivar toleranceUnit: defines how the massTolerance is applied to the mass value, 'ppm' * (1 +/- tolerance*1E-6) or 'da': +/- value
    :ivar rtTolerance: is used to calculate the retention time window to match :class:`FeatureItem` (is used as + / - rt tolerance)
    :ivar massKey: mass attribute key in :attr:`FeatureItem.__dict__`
    :ivar rtKey: retention time attribute key in :attr:`FeatureItem.__dict__`
    :ivar largerBetter: boolean, True if a larger PSM score is better,
    used to annotate feature group with the peptide sequence from the best scoring :class:`SpectrumIdentificationItem`
    :ivar matchCharge: If True only the charge from :attr:`FeatureItem.charge` is used for feature grouping, else all charge states
    from :ivar:`targetChargeStates` are used
    :ivar labelDescriptor: :class:`LabelDescriptor` describes the label setup of an experiment, must only be specified if
    peptides are labeled with stable isotopes.
    :ivar targetChargeStates: list of charge states which are used for feature Grouping if :ivar:`matchCharge` is False

    See also :class:`FeatureGroupContainer` and :class:`FeatureGroupItem`

    NOTE: TODO -> feature groups still can contain multiple features per charge/specfile/label position.
                  such entries should be removed or the group tagged as self.isValid = False
    """
    specfiles = featureContainer.specfiles if specfiles is None else aux.toList(specfiles)

    matchedFeatures = set()

    groupContainer = FeatureGroupContainer(specfiles, labelDescriptor=labelDescriptor)

    featureArrays = featureContainer.getArrays([rtKey, massKey, 'charge'], sort=massKey, specfiles=specfiles)
    for feature in featureContainer.getItems(filterAttribute=featureFilterAttribute, specfiles=specfiles):
        if feature.containerId in matchedFeatures:
            continue
        if matchCharge:
            targetChargeStates = [feature.charge]

        groupItem = FeatureGroupItem()
        matching = 'featureSeed'
        labelMassDiff = float()

        while matching is not False:
            if toleranceUnit.lower() == 'ppm':
                limMassUp = getattr(feature, massKey) * (1 + massTolerance * 1e-6) + labelMassDiff
                limMassDown = getattr(feature, massKey) * (1 - massTolerance * 1e-6) + labelMassDiff
            elif toleranceUnit.lower() == 'da':
                limMassUp = getattr(feature, massKey) + massTolerance + labelMassDiff
                limMassDown = getattr(feature, massKey) - massTolerance + labelMassDiff

            limRtUp = getattr(feature, rtKey) + rtTolerance
            limRtDown = getattr(feature, rtKey) - rtTolerance

            arrayPosLeft = bisect.bisect_left(featureArrays[massKey], limMassDown)
            arrayPosRight = bisect.bisect_right(featureArrays[massKey], limMassUp)

            matchMask = ((featureArrays[rtKey][arrayPosLeft:arrayPosRight] <= limRtUp) &
                         (featureArrays[rtKey][arrayPosLeft:arrayPosRight] >= limRtDown)
                         )

            chargeIdMatches = ddict(list)

            for targetCharge in targetChargeStates:
                chargMask = matchMask & (featureArrays['charge'][arrayPosLeft:arrayPosRight] == targetCharge)
                for containerId in featureArrays['containerId'][arrayPosLeft:arrayPosRight][chargMask]:
                    groupItem.featureIds.append(containerId)
                    chargeIdMatches[targetCharge].append(containerId)

            if matching == 'featureSeed':
                scorePepList = list()
                for containerId in groupItem.featureIds:
                    if featureContainer[containerId].isAnnotated:
                        scorePepList.append((featureContainer[containerId].score, featureContainer[containerId].peptide))
                scorePepList.sort(reverse=largerBetter)

                try:
                    peptide = scorePepList[0][1]
                    labelState = returnLabelState(peptide, labelDescriptor, labelSymbols=None)
                except IndexError:
                    labelState = -1

                for charge, idMatches in chargeIdMatches.items():
                    if charge not in groupItem.matchMatrix:
                        groupItem.matchMatrix[charge] = numpy.empty((len(specfiles), len(labelDescriptor.labels)+1), dtype='object')
                        groupItem.intensityMatrix[charge] = numpy.empty((len(specfiles), len(labelDescriptor.labels)+1), dtype='object')

                    for idMatch in idMatches:
                        specfilePos = groupContainer.specfilePositions[idMatch[0]]
                        groupItem.matchMatrix[charge][specfilePos][labelState] = idMatch
                        groupItem.intensityMatrix[charge][specfilePos][labelState] = featureContainer[idMatch].intensity

                #If the labelState is not -1, define the next labelState and expected labelMassDiff
                if labelState < 0:
                    matching = False
                else:
                    matching = 'labels'
                    labelMassDifferences = returnLabelStateMassDifferences(peptide, labelDescriptor)
                    labelState = labelMassDifferences.keys()[0]
                    labelMassDiff = labelMassDifferences[labelState]
                    del(labelMassDifferences[labelState])
            else:
                for charge, idMatches in chargeIdMatches.items():
                    if charge not in groupItem.matchMatrix:
                        groupItem.matchMatrix[charge] = numpy.empty((len(specfiles), len(labelDescriptor.labels)+1), dtype='object')
                        groupItem.intensityMatrix[charge] = numpy.empty((len(specfiles), len(labelDescriptor.labels)+1), dtype='object')

                    for idMatch in idMatches:
                        specfilePos = groupContainer.specfilePositions[idMatch[0]]
                        groupItem.matchMatrix[charge][specfilePos][labelState] = idMatch
                        groupItem.intensityMatrix[charge][specfilePos][labelState] = featureContainer[idMatch].intensity

                if len(labelMassDifferences) == 0:
                    matching = False
                else:
                    labelState = labelMassDifferences.keys()[0]
                    labelMassDiff = labelMassDifferences[labelState]
                    del(labelMassDifferences[labelState])

        mhList = list()
        rtList = list()
        siIds = list()
        siiIds = list()
        scorePepList = list()
        for featureId in groupItem.featureIds:
            rtList.append(featureContainer[featureId].rt)
            mhList.append(featureContainer[featureId].mh)
            siIds.extend(featureContainer[featureId].siIds)
            siiIds.extend(featureContainer[featureId].siiIds)
            if featureContainer[featureId].isAnnotated:
                scorePepList.append((featureContainer[featureId].score,
                                     featureContainer[featureId].peptide,
                                     featureContainer[featureId].sequence))
            if featureId not in groupContainer.index:
                groupContainer.index[featureId] = list()
            groupContainer.index[featureId].append(groupItem)
        scorePepList.sort(reverse=largerBetter)

        groupItem.rt = numpy.average(rtList)
        groupItem.mh = numpy.average(mhList)
        groupItem.siIds = siIds
        groupItem.siiIds = siiIds
        if len(groupItem.siIds) > 0:
            groupItem.isMatched = True
        else:
            groupItem.isMatched = False
        if len(scorePepList) > 0:
            groupItem.score = scorePepList[0][0]
            groupItem.peptide = scorePepList[0][1]
            groupItem.sequence = scorePepList[0][2]
            groupItem.isAnnotated = True
        else:
            groupItem.isMatched = False
        groupItem.isValid = True
        groupContainer.container.append(groupItem)
        matchedFeatures.update(groupItem.featureIds)
    return groupContainer
