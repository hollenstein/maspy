from __future__ import print_function

import bisect
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

from pymzml import run
from pyteomics import mass

import msfunctions.common as cFU


class ContainerItem(object):
    """
    :ivar containerId: used to look up item in ItemContainer.index
    :ivar id: identifier in original file
    :ivar specfile: spectrum filename
    :ivar isValid: should evaluate to True or False, None if unspecified - used to filter data.
    """
    def __init__(self, identifier, specfile):
        self.containerId  = (specfile, identifier)
        self.id = identifier
        self.specfile = specfile
        self.isValid = None


class ItemContainer(object):
    """
    :ivar index: used to look up items with their containerId
    :ivar container: {specfile:[ContainerItem(), ContainerItem(), ...]}
    :ivar specfiles: list of filenames of spectrum files
    """
    def __init__(self):
        self.index = dict()
        self.container = dict()
        self.specfiles = list()

    def __getitem__(self, key):
        """
        Return an item from index, using the containerId
        """
        return self.index[key]

    def getItems(self, specfiles=None, sort=None, reverse=False, filterAttribute='isValid', filterTargetValue=True, selector=None):
        """
        :param specfiles: filenames of spectrum files - return only items from those files. (str or [str, str...])
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

        specfiles = self.specfiles if specfiles == None else cFU.toList(specfiles)
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
        """
        Return a condensed array of data selected from ContainerItems for faster data processing.

        :param attributes: list of item attributes that should be written to the returned array.

        for the other parameters see getValidItems.

        :returns: dict(key1 from keylist: numpy.array, key2 from keylist: numpy.array, ..., indexPos: numpy.array, id: numpy.array, specfile: numpy.array),
        i.e. returns the columns of the table specified by the list of keys, each numpy.array has the dimensions Nx1. If a value is not present, None, is substituted.
        """
        items = self.getItems(specfiles=specfiles, sort=sort, reverse=reverse,
                              filterAttribute=filterAttribute, filterTargetValue=filterTargetValue,
                              selector=selector
                              )
        arrays = dict()
        attributes = set(['containerId', 'id', 'specfile'] + cFU.toList(attributes))
        for key in attributes:
            arrays[key] = list()

        for item in items:
            for key in attributes:
                arrays[key].append(getattr(item, key, None))

        for key in arrays.keys():
            arrays[key] = numpy.array(arrays[key])

        return arrays

    def save(self, fileFolder, fileName):
        """Store a pickled version of the self, using the __class__.__name__ as file-appendix."""
        fileName = '.'.join((fileName, self.__class__.__name__))
        filePath = os.path.join(fileFolder, fileName).replace('\\', '/')
        with open(filePath, 'w') as openFile:
            pickle.dump(self, openFile)

    @classmethod
    def load(cls, fileFolder, fileName):
        """Load a pickled version of the self, using the __class__.__name__ as file-appendix."""
        fileName = '.'.join((fileName, cls.__name__))
        filePath = os.path.join(fileFolder, fileName).replace('\\', '/')
        with open(filePath, 'r') as openFile:
            return pickle.load(openFile)


class SpectrumItem(ContainerItem):
    """Representation of a spectrum"""
    def __init__(self, identifier, specfile):
        super(SpectrumItem, self).__init__(identifier, specfile)
        self.msLevel = None


class SiContainer(ItemContainer):
    """
    ItemContainer for mass spectrometry data (spectra) (for example MS1, MS2),
    SiContainer ... Spectrum Item Container.
    see also 'class::SiiContainer' (Spectrum Identification Item Container) which contains sequence data.

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
    def load(cls, fileFolder, fileName, loadIonList=True):
        """Load a pickled version of the self, using the __class__.__name__ as file-appendix."""
        _siContainer = super(cls, cls).load(fileFolder, fileName)

        ionListFileName = '.'.join((fileName, 'ionList'))
        ionListFilePath = os.path.join(fileFolder, ionListFileName).replace('\\', '/')
        if loadIonList and os.path.isfile(ionListFilePath):
            importedArray = numpy.load(ionListFilePath)
            for key, mzList, iList in itertools.izip(importedArray[0], importedArray[1], importedArray[2]):
                _siContainer.ionLists[key] = dict()
                _siContainer.ionLists[key]['mz'] = mzList
                _siContainer.ionLists[key]['i'] = iList
        return _siContainer


def importSpectrumItems(_siContainer, specfilePath, specfile, msLevel=[1, 2], importIonList=False, mgfType=None):
    """ Import spectra from mzml or mgf files
    :param _siContainer: Spectra are added to to this instance of SiContainer()
    :param specfilePath: spectrum filename (.mzml or .mgf)
    :param specfile: name to represent this file in the SiContainer. Each filename
    can only occure once, therefore importing the same filename again is prevented
    :param msLevel: msLevels to load
    :param importIonList: bool whether the ion spectra should be loaded
    :param mgfType: if the file is of type '.mgf', and the mgf was generated by pParse set to mgfType="pParse" because of header information ambiguity.
    """
    if not os.path.isfile(specfilePath):
        print('File does not exits:', specfilePath)
    elif not specfilePath.lower().endswith('.mzml') and not specfilePath.lower().endswith('.mgf'):
        print('File is not an "mzml" or "mgf" file:', specfilePath)
    else:
        if specfile not in _siContainer.specfiles:
            _siContainer.specfiles.append(specfile)
            _siContainer.container[specfile] = list()
            if specfilePath.lower().endswith('.mzml'):
                msrun = pymzmlReadMzml(specfilePath)
                importMzmlSpectrumItems(_siContainer, msrun, specfile, msLevel=[1, 2], importIonList=importIonList)
            elif specfilePath.lower().endswith('.mgf'):
                importMgfSpectrumItems(_siContainer, specfilePath, specfile, importIonList=importIonList, mgfType=mgfType)
        else:
            print(specfile, 'is already present in the SiContainer, import interrupted.')


def importMzmlSpectrumItems(_siContainer, msrun, specfile, msLevel=[1, 2], importIonList=False):
    """load mzml spectrum items.

    Mostly used as a private function by importSpectrumItems.

    :param msrun: pymzml.msrun class instance of mzml file containing parameters to load using pymzml.
    :param specfile: filename (keyword-name not path) of spectrum file

    see also importSpectrumItems.
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
            _siContainer.container[specfile].append(si)
            _siContainer.index[si.containerId] = si

            if importIonList:
                _siContainer.ionLists[si.containerId] = dict()
                _siContainer.ionLists[si.containerId]['mz'] = numpy.array(spectrum.mz, dtype='float64')
                _siContainer.ionLists[si.containerId]['i'] = numpy.array(spectrum.i, dtype='float64')

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


def importMgfSpectrumItems(_siContainer, specfilePath, specfile, importIonList=False, mgfType=None):
    """Load mgf file spectrum items.

    Mostly used as a private function by importSpectrumItems.

    :param specfilePath: actual path to file.
    :param specfile: filename (keyword-name not path) of spectrum file
    see also importSpectrumItems.
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
                _siContainer.container[specfile].append(si)
                _siContainer.index[si.containerId] = si

                if importIonList:
                    mzList, iList = _splitMgfIonList(ionList)
                    _siContainer.ionLists[si.containerId] = dict()
                    _siContainer.ionLists[si.containerId]['mz'] = numpy.array(mzList, dtype='float64')
                    _siContainer.ionLists[si.containerId]['i'] = numpy.array(iList, dtype='float64')


class SpectrumIdentificationItem(ContainerItem):
    """Representation of a MS2 sequence annotation (Peptide Spectrum Match)"""
    def __init__(self, identifier, specfile):
        super(SpectrumIdentificationItem, self).__init__(identifier, specfile)


class SiiContainer(ItemContainer):
    def __init__(self):
        super(SiiContainer, self).__init__()

    def addSiInfo(self, siContainer, specfiles=None, attributes=['obsMz', 'rt', 'charge']):
        """ Copy attributes into sii from the corresponding SpectrumItem in siContainer,
        if an attribute is not presend in the SpectrumItem the attribute value is set to None
        Attribute examples: 'obsMz', 'rt', 'charge', 'TIC', 'IIT', 'ms1Id'
        """
        specfiles = self.specfiles if specfiles == None else cFU.toList(specfiles)

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


    def calcMz(self, specfiles=None, guessCharge=True):
        # Guess charge uses the calculated mass and the observed m/z value to calculate the charge
        specfiles = self.specfiles if specfiles is None else cFU.toList(specfiles)
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
                        guessedCharge = round(peptideMass / (sii.obsMz - cFU.atomicMassProton), 0)
                        sii.calcMz = calcMzFromMass(peptideMass, guessedCharge)
                        sii.charge = guessedCharge
        del(tempPeptideMasses)


def importPsmResults(_siiContainer, psmArrays, specfile, psmType='percolator', qValue=0.01):
    if specfile not in _siiContainer.container:
        _siiContainer.container[specfile] = list()
    if specfile not in _siiContainer.specfiles:
        _siiContainer.specfiles.append(specfile)

    if psmType == 'percolator':
        _importFromPercolatorArray(_siiContainer, psmArrays, specfile, qValueCutOff=qValue)

def _importFromPercolatorArray(_siiContainer, psmArrays, specfile, qValueCutOff=None):
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

        if sii.containerId in _siiContainer.index:
            siiList = _siiContainer.index[sii.containerId]
            sii.rank = len(siiList) + 1
        else:
            sii.rank = 1
            _siiContainer.index[sii.containerId] = list()

        if sii.rank == 1:
            if qValueCutOff is not None:
                if sii.qValue <= qValueCutOff:
                    sii.isValid = True
            else:
                sii.isValid = True

        _siiContainer.index[sii.containerId].append(sii)
        _siiContainer.container[specfile].append(sii)


