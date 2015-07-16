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
        for the other parameters see :class:`ItemContainer.getValidItems()`

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
            if key == 'containerId':
                emptyArray = numpy.empty(len(arrays[key]), dtype='object')
                emptyArray[:] = arrays[key]
                arrays[key] = emptyArray
            else:
                arrays[key] = numpy.array(arrays[key])

        return arrays

    def save(self, fileFolder, fileName):
        """Store a pickled version of self using  __class__.__name__ as file-appendix."""
        fileName = '.'.join((fileName, self.__class__.__name__))
        filePath = os.path.join(fileFolder, fileName).replace('\\', '/')
        with open(filePath, 'w') as openFile:
            pickle.dump(self, openFile)

    @classmethod
    def load(cls, fileFolder, fileName):
        """Load a pickled version of self using the __class__.__name__ as file-appendix."""
        fileName = '.'.join((fileName, cls.__name__))
        filePath = os.path.join(fileFolder, fileName).replace('\\', '/')
        with open(filePath, 'r') as openFile:
            return pickle.load(openFile)


class SpectrumItem(ContainerItem):
    """Representation of a spectrum."""
    def __init__(self, identifier, specfile):
        super(SpectrumItem, self).__init__(identifier, specfile)
        self.msLevel = None


class SiContainer(ItemContainer):
    """
    ItemContainer for mass spectrometry data (spectra) (for example MS1, MS2),
    SiContainer ... Spectrum Item Container.
    see also :class:`SiiContainer` (Spectrum Identification Item Container) which contains sequence data.

    :ivar ionLists: spectrum ion information, not loaded by default
    dict(containerId=dict(mz=nump.array([mass / charge, ...]), i=numpy.array([intensity, ...]))).
    """
    def __init__(self):
        super(SiContainer, self).__init__()
        self.ionLists = dict()

    def save(self, fileFolder, fileName, saveIonList=True):
        """Store a pickled version of the self, using the __class__.__name__ as file-appendix.

        Stores the ionList in a separate file with appendix '.ionList'.
        """
        ionLists = self.ionLists
        self.ionLists = dict()
        try:
            super(self.__class__, self).save(fileFolder, fileName)
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

            ionListFileName = '.'.join((fileName, 'ionList'))
            ionListFilePath = os.path.join(fileFolder, ionListFileName).replace('\\', '/')
            with open(ionListFilePath, 'wb') as openFile:
                numpy.save(openFile, ionArray)

    @classmethod
    def load(cls, fileFolder, fileName, importIonList=True):
        """Load a pickled version of the self, using the __class__.__name__ as file-appendix."""
        siContainer = super(cls, cls).load(fileFolder, fileName)

        ionListFileName = '.'.join((fileName, 'ionList'))
        ionListFilePath = os.path.join(fileFolder, ionListFileName).replace('\\', '/')
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


class SiiContainer(ItemContainer):
    """
    ItemContainer for msn spectrum identifications (Peptide Spectrum Matches),
    SiiContainer ... Spectrum Identification Item Container.
    see also :class:`SiContainer` (Spectrum Item Container) which contains spectrum data.
    """
    def __init__(self):
        super(SiiContainer, self).__init__()

    def addSiInfo(self, siContainer, specfiles=None, attributes=['obsMz', 'rt', 'charge']):
        """ Copy attributes into sii from the corresponding SpectrumItem in siContainer,
        if an attribute is not presend in the SpectrumItem the attribute value is set to None
        Attribute examples: 'obsMz', 'rt', 'charge', 'TIC', 'IIT', 'ms1Id'
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
                        tempPeptideMasses[peptide] = calcPeptideMass(peptide)
                    peptideMass = tempPeptideMasses[peptide]
                    if charge is not None:
                        sii.calcMz = calcMzFromMass(peptideMass, charge)
                    elif guessCharge:
                        guessedCharge = round(peptideMass / (sii.obsMz - aux.atomicMassProton), 0)
                        sii.calcMz = calcMzFromMass(peptideMass, guessedCharge)
                        sii.charge = guessedCharge
        del(tempPeptideMasses)


class FeatureItem(ContainerItem):
    """Representation of a peptide elution feature.

    ivar isMatched: None if unspecified, should be set to False on import, True if any Si or Sii elements could be matched
    ivar isAnnotated: None if unspecified, should be set to False on import, True if any Sii elements could be matched
    ivar siIds: containerId values of matched Si entries
    ivar siiIds: containerId values of matched Sii entries
    ivar peptide: peptide sequence of best scoring Sii match
    ivar sequence: plain amino acid sequence of best scoring Sii match, used to retrieve protein information
    ivar score: score of best scoring Sii match
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

    see also :class:`SiContainer` (Spectrum Item Container) which contains spectrum data.
    see also :class:`SiiContainer` (Spectrum Identification Item Container) which contains sequence data.
    """
    def __init__(self):
        super(FeatureContainer, self).__init__()


class Peptide(object):
    """Describes a peptide derived by one or more proteins.

    :ivar sequence: amino acid sequence of the peptide
    :ivar missedCleavage: number of missed cleavages, dependens on enzyme specificity
    :ivar proteinList: protein ids that generate this peptide under certain digest condition
    :ivar startPosDict: {proteinId:startPosition, ...} start position of peptide in protein
    :ivar endPosDict: {proteinId:startPosition, ...} end position of peptide in protein
    :ivar mass: peptide mass in Daltons
    :ivar length: number of amino acids
    """
    def __init__(self, sequence, mc=None):
        self.sequence = sequence
        self.missedCleavage = mc

        self.proteinList = list()
        self.startPosDict = dict()
        self.endPosDict = dict()

    @lazyAttribute
    def mass(self):
        return pyteomics.mass.calculate_mass(self.sequence, charge=0)

    @lazyAttribute
    def length(self):
        return len(self.sequence)


class PeptideContainer(object):
    """Container object for peptides :class:`Peptide`.

    :ivar peptides: {peptideSequence:Peptide(), peptideSequence:Peptide(), ...}
    """
    def __init__(self):
        self.peptides = dict()

    def __getitem__(self, peptide):
        """Return entry from self.peptides using peptide as key."""
        return self.peptides[peptide]


class Protein(object):
    """Describes a protein.

    ivar id: identifier of the protein eg. UniprotId
    ivar name: name of the protein
    ivar sequence: amino acid sequence of the protein
    ivar mass: protein mass in Daltons
    ivar length: number of amino acids
    """
    def __init__(self, sequence, identifier=str(), name=str()):
        self.id = identifier
        self.name = name
        self.sequence = sequence

    @lazyAttribute
    def mass(self):
        return pyteomics.mass.calculate_mass(self.sequence, charge=0)

    @lazyAttribute
    def length(self):
        return len(self.sequence)


class ProteinContainer(object):
    """Container object for proteins :class:`Protein`, protein entries can be accessd via id or name.

    ivar proteinIds: {proteinId:Protein(), proteinId:Protein()}
    ivar proteinNames: {proteinName:Protein(), proteinName:Protein()}
    """
    # A container for Protein or ProteinEvidence objects
    def __init__(self):
        self.proteinIds = dict()
        self.proteinNames = dict()

    def __getitem__(self, key):
        """Uses key to return protein entries :class:`Protein`.

        :ivar key: either a protein id or a protein name
        """
        if key in self.proteinIds:
            return self.proteinIds[key]
        elif key in self.proteinNames:
            return self.proteinNames[key]
        else:
            raise KeyError(key)


class PeptideEvidence(Peptide):
    """Summarizes all the evidence (:class:`SpectrumIdentificationItem`) for a certain peptide.

    ivar peptide: amino acid sequence of the peptide including modifications
    ivar sequence: amino acid sequence of the peptide, corresponds to peptideRef of mzidentml files
    ivar bestId: containerId of best scoring Sii item
    ivar siiIds: containerIds of all added Sii items
    ivar score: best score of all added Sii items
    ivar scores: scores of all added Sii items

    see also :class:`Peptide`
    """
    def __init__(self, peptide, sequence=None):
        sequence = removeModifications(peptide) if sequence is None else sequence
        super(PeptideEvidence, self).__init__(sequence)
        del(self.startPosDict)
        del(self.endPosDict)
        del(self.missedCleavage)

        self.peptide = peptide
        self.sequence = sequence
        self.bestId = tuple()
        self.siiIds = list()
        self.score = float()
        self.scores = list()

    @lazyAttribute
    def mass(self):
        return calcPeptideMass(peptide)


class PeptideEvidenceContainer(PeptideContainer):
    """Container for peptide evidence class:`PeptideEvidence`
    :ivar peptides: {peptide:PeptideEvidence(), peptide:PeptideEvidence(), ...}
    :ivar siiContainer: SiiContainer() which is used to generate the PeptideEvidence
    :ivar scoreKey: SpectrumIdentificationItem attribute which is used to find the best scoring item
    :ivar largerBetter: True if a larger value of the scoreKey attribute means a better score
    :ivar modified: True if modified peptides are treated as unique entries,
    set False to use only the amino acid sequence of a peptide
    see also :class:`PeptideContainer`
    """
    def __init__(self, siiContainer, scoreKey='qValue', largerBetter=False, modified=False):
        super(PeptideEvidenceContainer, self).__init__()
        self.siiContainer = siiContainer

        self._scoreKey = scoreKey
        self._largerBetter = largerBetter
        self._modified = modified

        self._generatePeptideEvidence()

    def _generatePeptideEvidence(self):
        if self._modified:
            peptideKey = 'peptide'
        else:
            peptideKey = 'sequence'

        for sii in self.siiContainer.getItems(sort=self._scoreKey, reverse=self._largerBetter):
            peptide = getattr(sii, peptideKey)
            siiScore = getattr(sii, self._scoreKey)
            if peptide in self.peptides:
                self.peptides[peptide].siiIds.append(sii.containerId)
                self.peptides[peptide].scores.append(siiScore)
            else:
                peptideEvidence = PeptideEvidence(peptide, sequence=sii.sequence)
                peptideEvidence.bestId = sii.containerId
                peptideEvidence.siiIds.append(sii.containerId)
                peptideEvidence.score = siiScore
                peptideEvidence.scores.append(siiScore)
                self.peptides[peptide] = peptideEvidence


class ProteinEvidence(Protein):
    """ Summarizes all the PeptideEvidence information for a certain protein
    see also :class:`Protein`
    :ivar id: amino acid sequence of the peptide including modifications
    :ivar uniquePeptides: amino acid sequence of the peptide, corresponds to peptideRef of mzidentml files
    :ivar sharedPeptides: containerId of best scoring Sii item
    :ivar uniquePsmCount: containerIds of all added Sii items
    :ivar sharedPsmCount: best score of all added Sii items
    :ivar isValid: should evaluate to True or False, None if unspecified - used to filter data.
    see also :class:`Protein`
    """
    def __init__(self, identifier, sequence=str(), name=None):
        super(ProteinEvidence, self).__init__(sequence, identifier=identifier, name=name)
        self.uniquePeptides = list()
        self.sharedPeptides = list()
        self.uniquePsmCount = int()
        self.sharedPsmCount = int()
        self.isValid = None

    @lazyAttribute
    def coverage(self):
        """Calculate the number of identified amino acids by unique peptides"""
        raise NotImplementedError


class ProteinEvidenceContainer(ProteinContainer):
    """ Container for protein evidence :class:`ProteinEvidence`.

    :ivar proteinIds: {proteinId:ProteinEvidence(), proteinId:ProteinEvidence()}
    :ivar proteinNames: {proteinName:ProteinEvidence(), proteinName:ProteinEvidence()}
    :ivar peptideEvidenceContainer: class:`PeptideEvidenceContainer` used to generate protein evidence
    :ivar proteindb: fasta representation of proteins :class:`ProteinContainer`
    :ivar peptidedb: fasta representation of peptides :class:`PeptideContainer`
    """
    def __init__(self, peptideEvidenceContainer, proteindb, peptidedb):
        super(ProteinEvidenceContainer, self).__init__()
        del(self.proteinIds)
        del(self.proteinNames)
        self.proteins = dict()
        self.peptideEvidenceContainer = peptideEvidenceContainer
        self.peptidedb = peptidedb
        self.proteindb = proteindb
        self.validProteinList = list()

        self._generateProteinEvidence()

    def __getitem__(self, key):
        """Uses key to return protein evidence entries :class:`ProteinEvidence`.

        :ivar key: proteinId
        """
        return self.proteins[key]

    def _generateProteinEvidence(self):
        for peptide, _peptideEvidence in self.peptideEvidenceContainer.peptides.items():
            isUnique = self.peptidedb[_peptideEvidence.sequence].unique
            proteins = set(self.peptidedb[_peptideEvidence.sequence].proteinList)
            for protein in proteins:
                if protein in self.proteins:
                    _proteinEvidence = self.proteins[protein]
                else:
                    _proteinEvidence = ProteinEvidence(proteindb[protein].id, sequence=proteindb[protein].sequence, name=proteindb[protein].name)
                    self.proteins[protein] = _proteinEvidence

                if isUnique:
                    _proteinEvidence.uniquePeptides.append(peptide)
                    _proteinEvidence.uniquePsmCount += len(_peptideEvidence.siiIds)
                else:
                    _proteinEvidence.sharedPeptides.append(peptide)
                    _proteinEvidence.sharedPsmCount += len(_peptideEvidence.siiIds)

        for _proteinEvidence in self.proteins.values():
            if len(_proteinEvidence.uniquePeptides) > 0:
                _proteinEvidence.valid = True
                self.validProteinList.append(_proteinEvidence.id)
            else:
                _proteinEvidence.valid = False


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
        specfilePath = aux.searchFileLocation(specfile, 'SiContainer', fileDirectory)
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


def generateSiContainerFiles(fileDirectory):
    """Generate SiContainer and ionList files for all mzML files in the fileDirectory and its subfolders.

    see also :func:`removeSiContainerFiles` and :meth:`SiContainer.save`
    """
    for filePath in aux.matchingFilePaths('', fileDirectory, targetFileExtension='mzML', selector=lambda x: True):
        dotPosition = [x for x in aux.findAllSubstrings(filePath, '.')][-1]

        fileFolder = os.path.dirname(filePath)
        fileName = os.path.basename(filePath[:dotPosition])
        targetFilePath = '.'.join((filePath[:dotPosition], 'SiContainer'))
        if not os.path.isfile(targetFilePath):
            siContainer = SiContainer()
            importSpectrumItems(siContainer, filePath, fileName, importIonList=True)
            print('Saving SiContainer / ionList :', fileName)
            siContainer.save(fileFolder, fileName)


def removeSiContainerFiles(fileDirectory):
    """Remove all SiContainer and ionList files in the fileDirectory and its subfolders.

    see also :func:`generateSiContainerFiles` and :meth:`SiContainer.save`
    """
    for filePath in aux.matchingFilePaths('', fileDirectory, targetFileExtension='SiContainer', selector=lambda x: True):
        os.remove(filePath)
    for filePath in aux.matchingFilePaths('', fileDirectory, targetFileExtension='ionList', selector=lambda x: True):
        os.remove(filePath)


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
                       ('MS:1000927', ['value'])
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
                                                  ('charge', 'MS:1000041', 2)
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
        tsvreader = csv.reader(openFile, delimiter="\t")
        headerLine = tsvreader.next()
        headerDict = dict([ [y,x] for (x,y) in enumerate( headerLine ) ])
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
                featureItem.rtLow = min(rtArea)
                featureItem.rtHigh = max(rtArea)
                featureItem.charge = featureEntryDict['charge']
                featureItem.mz = featureEntryDict['mz']
                featureItem.mh = aux.returnMh(featureEntryDict['mz'], featureEntryDict['charge'])
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

def returnDigestedFasta(filePath, minLength=5, maxLength=40, missedCleavage=2,
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
    proteindb = ProteinContainer()
    peptidedb = PeptideContainer()

    for fastaEntry in fastaRead:
        protein = Protein(fastaEntry['sequence'], identifier=fastaEntry['sysName'],
                          name = fastaEntry['stdName']
                          )
        proteindb.proteinIds[fastaEntry['sysName']] = protein
        proteindb.proteinNames[fastaEntry['stdName']] = protein

        for unmodPeptide, info in digestInSilico(fastaEntry['sequence'], missedCleavage,
                                                  removeNtermM=True, minLength=minLength,
                                                  maxLength=maxLength
                                                  ):
            if ignoreIsoleucine:
                unmodPeptideNoIsoleucine = unmodPeptide.replace('I', 'L')
                if unmodPeptideNoIsoleucine in peptidedb.peptides:
                    currPeptide = peptidedb[unmodPeptideNoIsoleucine]
                else:
                    currPeptide = Peptide(unmodPeptideNoIsoleucine, mc=info['missedCleavage'])
                    peptidedb.peptides[unmodPeptideNoIsoleucine] = currPeptide

                if unmodPeptide not in peptidedb.peptides:
                    peptidedb.peptides[unmodPeptide] = currPeptide
            else:
                if unmodPeptide in peptidedb.peptides:
                    currPeptide = peptidedb[unmodPeptide]
                else:
                    currPeptide = Peptide(unmodPeptide, mc=info['missedCleavage'])
                    peptidedb.peptides[unmodPeptide] = currPeptide

            currPeptide.proteinList.append(fastaEntry['sysName'])
            currPeptide.startPosDict[fastaEntry['sysName']] = info['startPos']
            currPeptide.endPosDict[fastaEntry['sysName']] = info['endPos']

    for peptide in peptidedb.peptides.keys():
        numProteinMatches = len(peptidedb[peptide].proteinList)
        if numProteinMatches == 1:
            peptidedb[peptide].unique = True
        elif numProteinMatches > 1:
            peptidedb[peptide].unique = False
        else:
            print('No protein matches in peptidedb for peptide sequence: ', peptide)

    return proteindb, peptidedb


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

    :param proteinSequence: amino acid sequence of the protein to be digested
    :param missedCleavages: number of allowed missed cleavage sites
    :param removeNtermM: If True, consider peptides with the n-terminal methionine of the protein removed
    :param minLength: only yield peptides with length >= minLength
    :param maxLength: only yield peptides with length <= maxLength

    NOTE: at the moment it only works for trypsin (K/R) and c-terminal cleavage
    """
    # Return in silico digested peptides, peptide start position, peptide end position
    # Peptide position start at 1 and end at len(proteinSequence)
    passFilter = lambda seq: (len(seq) >= minLength and len(seq) <= maxLength)

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
            sequence = proteinSequence[startPos:endPos]
            if passFilter(sequence):
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
            sequence = proteinSequence[startPos:endPos]
            if passFilter(sequence):
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
                sequence = proteinSequence[startPos:endPos]
                if passFilter(sequence):
                    info = dict()
                    info['startPos'] = startPos+1
                    info['endPos'] = endPos
                    info['missedCleavage'] = missedCleavage
                    yield sequence, info
        lastCleavagePos += 1


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
        unimodSymbol = '[UNIMOD:' + unimodNumber + ']'
        numMod = peptide.count(unimodSymbol)
        unmodPeptide = unmodPeptide.replace(unimodSymbol, '')
        additionalModMass += unimodMass * numMod

    if unmodPeptide.find('[') != -1:
        raise Exception()

    unmodPeptideMass = pyteomics.mass.calculate_mass(unmodPeptide, charge=0)
    modPeptideMass = unmodPeptideMass + additionalModMass
    return modPeptideMass


def calcMzFromMass(mass, charge):
    """Calculate the mz value of a peptide from its mass and charge.

    :type mass: float
    :type charge: int
    """
    mz = (mass + (aux.atomicMassProton * charge) ) / charge
    return mz


def calcMassFromMz(mz, charge):
    """Calculate the mass of a peptide from its mz and charge.

    :type mz: float
    :type charge: int
    """
    mass = (mz - aux.atomicMassProton) * charge
    return mass


def removeModifications(peptide):
    """Removes all modifications from a peptide string

    :ivar peptide: peptide sequence, modifications have to be written in the format "[modificationName]"

    :type peptide: str
    """
    unmodPeptide = peptide
    if unmodPeptide.find('.') != -1:
        unmodPeptide = unmodPeptide.split('.')[1]
    while unmodPeptide.find('[') != -1:
        unmodPeptide = unmodPeptide.split('[', 1)[0] + unmodPeptide.split(']', 1)[1]
    return unmodPeptide


def returnModPositions(peptide, indexStart=1, removeModString='UNIMOD:'):
    """Determines the amino acid positions of all present modifications.

    :ivar peptide: peptide sequence, modifications have to be written in the format "[modificationName]"
    :ivar indexStart: returned amino acids positions of the peptide start with this number (1st amino acid position = indexStart)
    :ivar removeModString: string to remove from the returned modification name

    :return: {modificationName:[position1, position2, ...], ...}

    TEST:
    peptide = 'GFHIHEFGDATN[UNIMOD:7]GC[UNIMOD:4]VSAGPHFN[UNIMOD:7]PFKK'
    returnModPositions(peptide) == {'4': [14], '7': [12, 22]}
    """
    unidmodPositionDict = dict()
    while peptide.find('[') != -1:
        currModification = peptide.split('[')[1].split(']')[0]
        currPosition = peptide.find('[') - 1
        if currPosition == -1: # move n-terminal modifications to first position
            currPosition = 0
        currPosition += indexStart

        peptide = peptide.replace('['+currModification+']', '', 1)

        currModification = currModification.replace(removeModString, '')
        unidmodPositionDict.setdefault(currModification,list())
        unidmodPositionDict[currModification].append(currPosition)
    return unidmodPositionDict


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
            #TODO REMOVE specId = tuple(specId)
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
                    #TODO: REMOVE featureId = tuple(featureId)
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
