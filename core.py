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
import pyteomics

import pyms.auxiliary as aux


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
        attributes = set(['containerId', 'id', 'specfile'] + aux.toList(attributes))
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


class SpectrumIdentificationItem(ContainerItem):
    """Representation of a msn sequence annotation (Peptide Spectrum Match)"""
    def __init__(self, identifier, specfile):
        super(SpectrumIdentificationItem, self).__init__(identifier, specfile)


class SiiContainer(ItemContainer):
    """
    ItemContainer for msn spectrum identifications (Peptide Spectrum Matches),
    SiiContainer ... Spectrum Identification Item Container.
    see also 'class::SiContainer' (Spectrum Item Container) which contains spectrum data.
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
    """Representation of a peptide elution feature
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
    """
    ItemContainer for peptide elution features,
    see also 'class::SiContainer' (Spectrum Item Container) which contains spectrum data.
    see also 'class::SiiContainer' (Spectrum Identification Item Container) which contains sequence data.
    """
    def __init__(self):
        super(FeatureContainer, self).__init__()


###############################################
### Import functions for various file types ###
###############################################
def pymzmlReadMzml(mzmlPath):
    """Auxiliary function, specifies extra accesions to read from an mzml file, returns a pymzml run object and """
    extraAccessions = [('MS:1000827', ['value']),
                       ('MS:1000828', ['value']),
                       ('MS:1000829', ['value']),
                       ('MS:1000744', ['value']),
                       ('MS:1000016', ['value']),
                       ('MS:1000927', ['value'])
                       ]
    return pymzml.run.Reader(mzmlPath, extraAccessions=extraAccessions)


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


def importPsmResults(_siiContainer, fileLocation, specfile, psmType='percolator', psmEngine='comet', qValue=0.01):
    """
    Function to control the import of PSM results into a SiiContainer()
    See also _importFromPercolatorArray()
    """
    if specfile not in _siiContainer.container:
        _siiContainer.container[specfile] = list()
    if specfile not in _siiContainer.specfiles:
        _siiContainer.specfiles.append(specfile)

    if psmType == 'percolator':
        _psmArrays = _importPercolatorResults(fileLocation, psmEngine=psmEngine)
        _importFromPercolatorArrays(_siiContainer, _psmArrays, specfile, qValueCutOff=qValue)


def _importFromPercolatorArrays(_siiContainer, psmArrays, specfile, qValueCutOff=None):
    """
    Write Spectrum Identification Items into the siiContainer
    See also: _importPercolatorResults()
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


def _importPercolatorResults(fileLocation, psmEngine=None):
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


def importPeptideFeatures(_featureContainer, filelocation, specfile):
    """ Import peptide features from a featureXml file (eg. generated by OPENMS featureFinderCentroided)
    :param _featureContainer: Spectra are added to to this FeatureContainer()
    :param filelocation: file path of the file to import
    :param specfile: name to represent this file in the FeatureContainer. Each filename
    can only occure once, therefore importing the same filename again is prevented
    """
    if not os.path.isfile(filelocation):
        print('File does not exits:', filelocation)
    elif not filelocation.lower().endswith('.featurexml'):
        print('File is not a featurexml file:', filelocation)
    else:
        if specfile not in _featureContainer.specfiles:
            _featureContainer.specfiles.append(specfile)
            _featureContainer.container[specfile] = list()
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

                _featureContainer.index[featureItem.containerId] = featureItem
                _featureContainer.container[specfile].append(featureItem)
        else:
            print(specfile, 'is already present in the SiContainer, import interrupted.')


def _importFeatureXml(fileLocation):
    """Reads a featureXml file and returns {featureKey1: {attribute1:value1, attribute2:value2, ...}, ...}"""
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


################################################
### Functions to work with peptide sequences ###
################################################
def calcPeptideMass(peptide):
    """Calculate the mass of a peptide, modifications have to be in unimod format [UNIMOD:x]
    and 'x' has to be present in aux.unimodToMassDict

    :type peptide: str
    """
    unimodMassDict = aux.unimodToMassDict

    additionalModMass = float()
    unmodPeptide = peptide
    for unimodNumber, unimodMass in unimodMassDict.items():
        unimodSymbol = '[UNIMOD:' + unimodNumber + ']'
        numMod = peptide.count(unimodSymbol)
        unmodPeptide = unmodPeptide.replace(unimodSymbol, '')
        additionalModMass += unimodMass * numMod

    if unmodPeptide.find('UNIMOD') != -1:
        raise Exception()

    unmodPeptideMass = pyteomics.mass.calculate_mass(unmodPeptide, charge=0)
    modPeptideMass = unmodPeptideMass + additionalModMass
    return modPeptideMass


def calcMzFromMass(mass, charge):
    """Calculate the mz value of a peptide from mass and charge.

    :type mass: float
    :type charge: int
    """
    mz = (mass + (aux.atomicMassProton * charge) ) / charge
    return mz


def calcMassFromMz(mz, charge):
    """Calculate the mass of a peptide from mz and charge.

    :type mz: float
    :type charge: int
    """
    mass = (mz - aux.atomicMassProton) * charge
    return mass


def removeModifications(peptide):
    """Removes all '[x]' tags from a peptide; x can be a string of any length

    :type peptide: str
    """
    unmodPeptide = peptide
    if unmodPeptide.find('.') != -1:
        unmodPeptide = unmodPeptide.split('.')[1]
    while unmodPeptide.find('[') != -1:
        unmodPeptide = unmodPeptide.split('[', 1)[0] + unmodPeptide.split(']', 1)[1]
    return unmodPeptide


def returnModPositions(peptide, indexStart=1, removeModString='UNIMOD:'):
    """ returns a dictionary: key = modification, value = list of positions, positions start at var indexStart
    #test:
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
