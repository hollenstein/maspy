from __future__ import print_function, division
from future.utils import viewkeys, viewvalues, viewitems, listvalues, listitems

try: # python 2.7
    import cPickle as pickle
except ImportError: # python 3.x series
    import pickle

import bisect
from collections import defaultdict as ddict
import functools
import io
import operator
import os
from matplotlib import pyplot as plt
import numpy

import maspy.auxiliary as aux
import maspy.core
import maspy.sil

class FeatureGroupItem(object):
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
    """
    def __init__(self):
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


class FeatureGroupContainer(object):
    """ItemContainer for peptide feature groups :class`FeatureGroupItem`.

    :ivar container: Storage list of :class:`FeatureGroupItem`
    :ivar index: Use :attr:`FeatureItem.containerId` to which :class:`FeatureGroupItem` the feature was grouped
    :ivar labelDescriptor: :class:`maspy.sil.LabelDescriptor` describes the label setup of an experiment
    :ivar specfiles: List of keywords (filenames) representing files
    :ivar specfilePositions: {specfile:arrayPosition, ...}
    arrayPosition respresents the array position of a specfile in :attr:`FeatureGroupItem.matchMatrix`
    """
    def __init__(self, specfiles, labelDescriptor=None):
        self.index = dict()
        self.specfiles = specfiles
        self.container = list()
        self.labelDescriptor = maspy.sil.LabelDescriptor() if labelDescriptor is None else labelDescriptor

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

        :param report: 'lfq' or 'sil'

        :returns: dict(key1 from keylist: numpy.array, key2 from keylist: numpy.array, ..., indexPos: numpy.array, id: numpy.array, specfile: numpy.array),
        i.e. returns the columns of the table specified by the list of keys, each numpy.array has the dimensions Nx1. If a value is not present, None, is substituted.
        """
        #TODO: return string is not clear, docstring is not clear
        #TODO: explain lfq/sil report type

        items = self.getItems(sort=sort, reverse=reverse, filterAttribute=filterAttribute,
                              filterTargetValue=filterTargetValue, selector=selector
                              )

        arrays = dict([(key, []) for key in attributes])

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
            for labelState in list(viewkeys(self.labelDescriptor.labels)) + [-1]:
                labelAttributeName = ' '.join(('label:', str(labelState)))
                arrays[labelAttributeName] = list()
                reportAttributes.append(labelAttributeName)

        if report == 'sil':
            for item in items:
                for charge in viewkeys(item.intensityMatrix):
                    for specfile, specfilePosition in viewitems(self.specfilePositions):
                        for key in attributes:
                            arrays[key].append(getattr(item, key, None))
                        arrays['charge'].append(charge)
                        arrays['specfile'].append(specfile)
                        for labelState in list(viewkeys(self.labelDescriptor.labels)) + [-1]:
                            labelAttributeName = ' '.join(('label:', str(labelState)))
                            arrays[labelAttributeName].append(item.intensityMatrix[charge][specfilePosition, labelState])
        elif report == 'lfq':
            for item in items:
                for charge in viewkeys(item.intensityMatrix):
                    for labelState in list(viewkeys(self.labelDescriptor.labels)) + [-1]:
                        for key in attributes:
                            arrays[key].append(getattr(item, key, None))
                        arrays['charge'].append(charge)
                        arrays['labelState'].append(labelState)
                        for specfile, specfilePosition in viewitems(self.specfilePositions):
                            arrays[specfile].append(item.intensityMatrix[charge][specfilePosition, labelState])
        else:
            raise Exception('report must be either "lfq" or "sil", not '+report)

        for key in list(viewkeys(arrays)):
            if key in reportAttributes:
                arrays[key] = numpy.array(arrays[key], dtype=numpy.float64)
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
        #TODO: replace by jjh safe replace method
        filename = '.'.join((filename, self.__class__.__name__.lower()))
        filepath = aux.joinpath(filefolder, filename)
        with io.open(filepath, 'wb') as openFile:
            pickle.dump(self, openFile)

    @classmethod
    def load(cls, filefolder, filename):
        """Load a pickled version of self using __class__.__name__.lower() as file-appendix.

        :ivar filefolder: folder where the container has been saved
        :ivar filename: filename of the stores container, without file appendix
        """
        classInstance = cls._load(filefolder, filename)
        classInstance.index = dict()
        for item in classInstance.container:
            for featureId in item.featureIds:
                if featureId not in classInstance.index:
                    classInstance.index[featureId] = list()
                classInstance.index[featureId].append(item)
        return classInstance

    @classmethod
    def _load(cls, filefolder, filename):
        filename = '.'.join((filename, cls.__name__.lower()))
        filepath = aux.joinpath(filefolder, filename)
        with io.open(filepath, 'rb') as openFile:
            return pickle.load(openFile)

    def __str__(self):
        numSpecfiles = len(self.specfiles)
        numItems = len(self.container)

        output = [str(self.__class__)]
        output.append(' '.join(['Containing', str(numSpecfiles), 'specfiles and', str(numItems), 'items.']))
        output.extend(self.specfiles)
        return '\n'.join(output)

    def __repr__(self):
        return self.__str__()


# --- Functions to work with FeatureItems --- #
def matchToFeatures(featureContainer, specContainer, specfiles=None, fMassKey='mz', sMassKey='obsMz', isotopeErrorList=(0),
                    precursorTolerance=5, toleranceUnit='ppm', rtExpansionUp=0.10, rtExpansionDown=0.05, matchCharge=True,
                    scoreKey='pep', largerBetter=False):
    """Annotate :class:`FeatureItem` in :class:`FeatureContainer` by matching :class:`SpectrumItem` (Si) or :class:`SpectrumIdentificationItem` (Sii).

    :ivar featureContainer: :class:`maspy.core.FeatureItem` of this instance of :class:`maspy.core.FeatureContainer` are annotated
    :type specContainer: :class:`maspy.core.SiContainer` or :class:`maspy.core.SiiContainer`
    :ivar specfiles: Annotate only items of :attr:`FeatureContainer.container[specfile]`,
    if None all specfiles present in featureContainer and specContainer are processed
    :type specfiles: str, list or None
    :ivar fMassKey: mass attribute key in :attr:`FeatureItem.__dict__`
    :ivar sMassKey: mass attribute key in :attr:`SpectrumItem.__dict__` or :attr:`SpectrumIdentificationItem.__dict__` (eg 'obsMz', 'calcMz')
    :ivar isotopeErrorList: allowed isotope errors relative to the spectrum mass. Eg. [0, 1] if no feature has been matched with isotope error 0,
    the spectrum mass is increased by the mass difference of carbon isotopes 12 and 13 and matched again. The different isotope error values are
    tested in the specified order therefore 0 should normally be the 1st value of the tuple
    :type isotopeErrorList: list or tuple of int
    :ivar precursorTolerance: is the largest allowed mass deviation of Si or Sii to :class:`FeatureItem`
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
    #TODO: this function is nested and should maybe be rewritten
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
        specArrays = specContainer.getArrays([sMassKey, 'rt', 'charge', 'msLevel'], specfiles=specfile)
        featureArrays = featureContainer.getArrays(['rtHigh', 'rtLow', 'charge', fMassKey],
                                                   specfiles=specfile, sort=fMassKey
                                                   )
        featureArrays['rtHighExpanded'] = featureArrays['rtHigh'] + (featureArrays['rtHigh'] - featureArrays['rtLow']) * rtExpansionUp
        featureArrays['rtLowExpanded'] = featureArrays['rtLow'] - (featureArrays['rtHigh'] - featureArrays['rtLow']) * rtExpansionDown

        specFeatureDict = dict() ## key = scanNr, value = set(featureKeys)
        featureSpecDict = dict() ## key = featureKey, value = set(scanNrs)

        for specPos, specId in enumerate(specArrays['containerId']):
            specZ = specArrays['charge'][specPos]
            if specZ is None:
                continue
            specMass = specArrays[sMassKey][specPos]
            specRt = specArrays['rt'][specPos]

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
        stats['totalSpectra'] = len(specArrays['containerId'][(specArrays['msLevel'] != 1)])
        stats['matchedSpectra'] = len(specFeatureDict)
        stats['relMatchedSpectra'] = round(1.*stats['matchedSpectra']/stats['totalSpectra'], 3)

        print('------', specfile, '------')
        print('Annotated features:\t\t\t', stats['matchedFeatures'], '/', stats['totalFeatures'], '=', stats['relMatchedFeatures'], '%')
        print('Spectra matched to features:\t\t', stats['matchedSpectra'], '/', stats['totalSpectra'], '=', stats['relMatchedSpectra'], '%')
        if multiMatchCounter != 0:
                print('Discarded because of multiple matches:\t', multiMatchCounter)
        if isotopeErrorMatchCounter != 0:
                print('Isotope error matched spectra:\t\t', isotopeErrorMatchCounter)

        #annotate feature with sii information (peptide, sequence, score)
        if isinstance(specContainer, maspy.core.SiiContainer):
            for featureId in viewkeys(featureSpecDict):
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


def rtCalibration(featureContainer, allowedRtDev=60, allowedMzDev=2.5, showPlots=False, reference=None, specfiles=None):
    """Performs a retention time calibration between :class:`FeatureItem` of multiple specfiles.

    :ivar featureContainer: Perform alignment on :class:`FeatureItem` in :attr:`FeatureContainer.specfiles`
    :ivar allowedRtDev: maxium retention time difference of two features in two runs to be matched
    :ivar allowedMzDev: maxium relative m/z difference (in ppm) of two features in two runs to be matched
    :ivar showPlots: boolean, True if a plot should be generated which shows to results of the calibration
    :ivar reference: Can be used to specifically specify a reference specfile
    :ivar specfiles: Limit alignment to those specfiles in the featureContainer
    """
    #TODO: long function, maybe split into subfunctions
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
            for key in list(viewkeys(featureArrays)):
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

        for featurePos in range(len(featureArrays[mzKey])):
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

        splineInitialKnots = int(max(rtPosList) - min(rtPosList))
        dataFit = aux.DataFit(rtDevList, rtPosList)
        dataFit.splineInitialKnots = splineInitialKnots
        dataFit.splineTerminalExpansion = 0.2
        dataFit.processInput(dataAveraging='median', windowSize=10)
        dataFit.generateSplines()

        if showPlots:
            corrDevArr = rtDevList - dataFit.corrArray(rtPosList)
            timePoints = [min(rtPosList) + x for x in range(int(max(rtPosList)-min(rtPosList)))]
            corrValues  = dataFit.corrArray(timePoints)
            fig, ax = plt.subplots(3, 2, sharex=False, sharey=False, figsize=(20, 18))
            fig.suptitle(specfile)
            ax[0][0].hist(rtDevList, bins=100, color='grey', alpha=0.5, label='observed')
            ax[0][0].hist(corrDevArr, bins=100, color='red', alpha=0.5, label='corrected')
            ax[0][0].set_title('Retention time deviation')
            ax[0][0].legend()
            ax[0][1].hist(mzDevRelList, bins=100, color='grey')
            ax[0][1].set_title('Mz deviation [ppm]')
            ax[1][0].scatter(rtPosList, rtDevList, color='grey', alpha=0.1, label='observed')
            ax[1][0].plot(timePoints,corrValues, color='red', alpha=0.5, label='correction function')
            ax[1][0].set_title('Retention time deviation over time')
            ax[1][0].legend()
            ax[1][1].scatter(rtPosList, mzDevRelList, color='grey', alpha=0.1)
            ax[1][1].set_title('Mz deviation over time')
            ax[2][0].scatter(rtPosList, corrDevArr, color='grey', alpha=0.1)
            ax[2][0].set_title('Aligned retention time deviation over time')
            plt.show()

        featureArrays = featureContainer.getArrays(['rt'], specfiles=specfile, sort='rt', filterAttribute=None)
        featureArrays['corrRt'] = featureArrays['rt'] - dataFit.corrArray(featureArrays['rt'])
        for featureId, corrRt, rt in zip(featureArrays['containerId'], featureArrays['corrRt'], featureArrays['rt']):
            featureContainer.index[featureId].rt = corrRt


def groupFeatures(featureContainer, specfiles=None, featureFilterAttribute='isAnnotated', massTolerance=3,
                  toleranceUnit='ppm', rtTolerance=10, massKey='mh', rtKey='rt', largerBetter=False,
                  matchCharge=False, labelDescriptor=maspy.sil.LabelDescriptor(), targetChargeStates=(2, 3, 4, 5, 6)
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
    :ivar labelDescriptor: :class:`maspy.sil.LabelDescriptor` describes the label setup of an experiment, must only be specified if
    peptides are labeled with stable isotopes.
    :ivar targetChargeStates: list of charge states which are used for feature Grouping if :ivar:`matchCharge` is False

    See also :class:`FeatureGroupContainer` and :class:`FeatureGroupItem`
    """
    #TODO: -> feature groups still can contain multiple features per charge/specfile/label position.
    # such entries should be removed or the group tagged as self.isValid = False
    #TODO: function nested and too long, split in subfunctions
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
                    labelState = maspy.sil.returnLabelState(peptide, labelDescriptor, labelSymbols=None)
                except IndexError:
                    labelState = -1

                for charge, idMatches in viewitems(chargeIdMatches):
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
                    labelMassDifferences = maspy.sil.returnLabelStateMassDifferences(peptide, labelDescriptor)
                    labelState = list(viewkeys(labelMassDifferences))[0] #TODO: this is now an unsorted list of keys, if they need to be sorted use sort()
                    labelMassDiff = labelMassDifferences[labelState]
                    del(labelMassDifferences[labelState])
            else:
                for charge, idMatches in viewitems(chargeIdMatches):
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
                    labelState = list(viewkeys(labelMassDifferences))[0]
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
            groupItem.isAnnotated = False
        groupItem.isValid = True
        groupContainer.container.append(groupItem)
        matchedFeatures.update(groupItem.featureIds)
    return groupContainer
