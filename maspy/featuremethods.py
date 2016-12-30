"""
#TODO: module description
"""

#  Copyright 2015-2017 David M. Hollenstein, Jakob J. Hollenstein
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.

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

import bisect
from collections import defaultdict as ddict
import functools
import io
import operator
import os

#TODO: remove matplotlib from this modul, return a report object
from matplotlib import pyplot as plt
import numpy

import maspy.auxiliary as aux
import maspy.core
import maspy.sil


# --- Functions to work with FeatureItems --- #
def matchToFeatures(fiContainer, specContainer, specfiles=None, fMassKey='mz',
                    sMassKey='obsMz', isotopeErrorList=(0),
                    precursorTolerance=5, toleranceUnit='ppm',
                    rtExpansionUp=0.10, rtExpansionDown=0.05, matchCharge=True,
                    scoreKey='pep', largerBetter=False):
    """Annotate :class:`Fi <maspy.core.Fi>` (Feature items) by matching
    :class:`Si <maspy.core.Si>` (Spectrum items) or :class:`Sii
    <maspy.core.Sii>` (Spectrum identification items).

    :param fiContainer: :class:`maspy.core.FeatureContainer`, contains ``Fi``.
    :param specContainer: :class:`maspy.core.MsrunContainer` or
        :class:`maspy.core.SiiContainer`, contains ``Si`` or ``Sii``.
    :param specfiles: filenames of ms-run files, if specified consider only
        items from those files
    :type specfiles: str, list or None
    :param fMassKey: mass attribute key in :attr:`Fi.__dict__`
    :param sMassKey: mass attribute key in :attr:`Si.__dict__` or
        :attr:`Sii.__dict__` (eg 'obsMz', 'excMz')
    :param isotopeErrorList: allowed isotope errors relative to the spectrum
        mass, for example "0" or "1". If no feature has been matched with
        isotope error 0, the spectrum mass is increased by the mass difference
        of carbon isotopes 12 and 13 and matched again. The different isotope
        error values are tested in the specified order therefore "0" should
        normally be the first value of the list.
    :type isotopeErrorList: list or tuple of int
    :param precursorTolerance: the largest allowed mass deviation of ``Si`` or
        ``Sii`` relative to ``Fi``
    :param toleranceUnit: defines how the ``precursorTolerance`` is applied to
        the mass value of ``Fi``. ``"ppm": mass * (1 +/- tolerance*1E-6)`` or
        ``"da": mass +/- value``
    :param rtExpansionUp: relative upper expansion of ``Fi`` retention time
        area. ``limitHigh = Fi.rtHigh + (Fi.rtHigh - Fi.rtLow) * rtExpansionUp``
    :param rtExpansionDown: relative lower expansion of ``Fi`` retention time
        area. ``limitLow = Fi.rtLow - (Fi.rtHigh - Fi.rtLow) * rtExpansionDown``
    :param matchCharge: bool, True if ``Fi`` and ``Si`` or ``Sii`` must have the
        same ``charge`` state to be matched.
    :param scoreKey: ``Sii`` attribute name used for scoring the identification
        reliability
    :param largerBetter: bool, True if higher score value means a better
        identification reliability

    .. note:
        Concerning the feature retention area expansion. If ``Si`` or ``Sii`` is
        matched to multiple ``Fi`` the rt expansion is removed and the matching
        is repeated.

    .. note:
        If the ``specContainer`` is a ``SiiContainer`` then matched ``Fi`` are
        annotated with :attr:`Sii.peptide`, if multiple ``Sii`` are matched to
        ``Fi`` the one with the best score is used.

    #TODO: this function is nested pretty badly and should maybe be rewritten
    #TODO: replace tolerance unit "ppm" by tolerance mode "relative" and change
        repsective calculations
    """
    isotopeErrorList = aux.toList(isotopeErrorList)

    if specContainer.__class__.__name__ == 'MsrunContainer':
        listKeySpecIds = 'siIds'
    else:
        listKeySpecIds = 'siiIds'
    specContainerSpecfiles = [_ for _ in viewkeys(specContainer.info)]

    if specfiles is not None:
        specfiles = aux.toList(specfiles)
    else:
        specfiles = [_ for _ in viewkeys(fiContainer.info)]
    specfiles = list(set(specfiles).intersection(set(specContainerSpecfiles)))

    for specfile in specfiles:
        multiMatchCounter = int()
        isotopeErrorMatchCounter = int()
        specArrays = specContainer.getArrays([sMassKey, 'rt', 'charge',
                                              'msLevel'], specfiles=specfile
                                              )
        featureArrays = fiContainer.getArrays(['rtHigh', 'rtLow', 'charge',
                                               fMassKey], specfiles=specfile,
                                               sort=fMassKey
                                              )
        featureArrays['rtHighExpanded'] = (featureArrays['rtHigh'] +
                                           (featureArrays['rtHigh'] -
                                            featureArrays['rtLow']) *
                                           rtExpansionUp
                                           )
        featureArrays['rtLowExpanded'] = (featureArrays['rtLow'] -
                                          (featureArrays['rtHigh'] -
                                           featureArrays['rtLow']) *
                                          rtExpansionDown
                                          )

        specFeatureDict = dict() ## key = scanNr, value = set(featureKeys)
        featureSpecDict = dict() ## key = featureKey, value = set(scanNrs)

        for specPos, specId in enumerate(specArrays['id']):
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
                    specMassHigh = ((specMass + isotopeError * 1.003355 / specZ)
                                    * (1 + precursorTolerance*1E-6)
                                    )
                    specMassLow = ((specMass + isotopeError * 1.003355 / specZ)
                                   * (1 - precursorTolerance*1E-6)
                                   )
                elif toleranceUnit.lower() == 'da':
                    specMassHigh = ((specMass + isotopeError * 1.003355 / specZ)
                                    + precursorTolerance
                                    )
                    specMassLow  = ((specMass + isotopeError * 1.003355 / specZ)
                                    - precursorTolerance
                                    )

                posL = bisect.bisect_left(featureArrays[fMassKey],
                                          specMassLow
                                          )
                posR = bisect.bisect_right(featureArrays[fMassKey],
                                           specMassHigh
                                           )

                if matchCharge:
                    chargeMask = (featureArrays['charge'][posL:posR] == specZ)

                fRtHighKey = 'rtHighExpanded'
                fRtLowKey = 'rtLowExpanded'
                for fRtHighKey, fRtLowKey in [('rtHighExpanded',
                                               'rtLowExpanded'),
                                              ('rtHigh', 'rtLow')
                                              ]:
                    rtMask = ((featureArrays[fRtLowKey][posL:posR] <= specRt) &
                              (featureArrays[fRtHighKey][posL:posR] >= specRt)
                              )
                    if matchCharge:
                        matchedFeatureIds = featureArrays['id'][posL:posR][rtMask & chargeMask]
                    else:
                        matchedFeatureIds = featureArrays['id'][posL:posR][rtMask]

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
                    getattr(fiContainer.container[specfile][featureId],
                            listKeySpecIds
                            ).append(specId)
                    fiContainer.container[specfile][featureId].isMatched = True
                    specFeatureDict[specId] = featureId
                    featureSpecDict[featureId] = specId

        stats = dict()
        stats['totalFeatures'] = len(featureArrays['id'])
        stats['matchedFeatures'] = len(featureSpecDict)
        stats['relMatchedFeatures'] = round(100*stats['matchedFeatures']/stats['totalFeatures'], 1)
        stats['totalSpectra'] = len(specArrays['id'][(specArrays['msLevel'] != 1)])
        stats['matchedSpectra'] = len(specFeatureDict)
        stats['relMatchedSpectra'] = round(100*stats['matchedSpectra']/stats['totalSpectra'], 1)

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
                for specId in fiContainer.container[specfile][featureId].siiIds:
                    _sii = specContainer.getValidItem(specfile, specId)
                    score = getattr(_sii, scoreKey)
                    peptide = _sii.peptide
                    sequence = _sii.sequence
                    matches.append([score, peptide, sequence])
                matches.sort(reverse=largerBetter)

                fiContainer.container[specfile][featureId].isAnnotated = True
                fiContainer.container[specfile][featureId].score = matches[0][0]
                fiContainer.container[specfile][featureId].peptide = matches[0][1]
                fiContainer.container[specfile][featureId].sequence = matches[0][2]


def rtCalibration(fiContainer, allowedRtDev=60, allowedMzDev=2.5,
                  reference=None, specfiles=None, showPlots=False,
                  plotDir=None, minIntensity=1e5):
    """Performs a retention time calibration between :class:`FeatureItem` of multiple specfiles.

    :ivar fiContainer: Perform alignment on :class:`FeatureItem` in :attr:`FeatureContainer.specfiles`
    :ivar allowedRtDev: maxium retention time difference of two features in two runs to be matched
    :ivar allowedMzDev: maxium relative m/z difference (in ppm) of two features in two runs to be matched
    :ivar showPlots: boolean, True if a plot should be generated which shows to results of the calibration
    :ivar plotDir: if not None and showPlots is True, the plots are saved to
        this location.
    :ivar reference: Can be used to specifically specify a reference specfile
    :ivar specfiles: Limit alignment to those specfiles in the fiContainer
    :ivar minIntensity: consider only features with an intensity above this value
    """
    #TODO: long function, maybe split into subfunctions
    specfiles = [_ for _ in viewkeys(fiContainer.info)] if specfiles is None else specfiles
    matchCharge = True

    refMzKey = 'mz'
    mzKey = 'mz'

    if reference is not None:
        if reference in specfiles:
            specfiles = [reference] + list(set(specfiles).difference(set([reference])))
        else:
            print('Specified reference specfile not present, using reference: ', specfiles[0])

    for featureItem in fiContainer.getItems(specfiles=specfiles):
        if not hasattr(featureItem, 'obsRt'):
            setattr(featureItem, 'obsRt', featureItem.rt)

    referenceArrays = None
    for specfile in specfiles:
        featureArrays = fiContainer.getArrays(['rt', 'charge', 'mz', 'intensity'],
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
            ax[0][0].set_xlim(allowedRtDev*-1, allowedRtDev)
            ax[0][1].hist(mzDevRelList, bins=100, color='grey')
            ax[0][1].set_title('Mz deviation [ppm]')
            ax[1][0].scatter(rtPosList, rtDevList, color='grey', alpha=0.1, label='observed')
            ax[1][0].plot(timePoints,corrValues, color='red', alpha=0.5, label='correction function')
            ax[1][0].set_title('Retention time deviation over time')
            ax[1][0].legend()
            ax[1][0].set_ylim(allowedRtDev*-1, allowedRtDev)
            ax[1][1].scatter(rtPosList, mzDevRelList, color='grey', alpha=0.1)
            ax[1][1].set_title('Mz deviation over time')
            ax[1][1].set_ylim(allowedMzDev*-1, allowedMzDev)
            ax[2][0].scatter(rtPosList, corrDevArr, color='grey', alpha=0.1)
            ax[2][0].set_title('Aligned retention time deviation over time')
            ax[2][0].set_ylim(allowedRtDev*-1, allowedRtDev)
            if plotDir is not None:
                plotloc = aux.joinpath(plotDir, specfile+'.rtAlign.png')
                fig.savefig(plotloc)
            else:
                fig.show()

        featureArrays = fiContainer.getArrays(['rt'], specfiles=specfile, sort='rt')
        featureArrays['corrRt'] = featureArrays['rt'] - dataFit.corrArray(featureArrays['rt'])
        for featureId, corrRt, rt in zip(featureArrays['id'], featureArrays['corrRt'], featureArrays['rt']):
            fiContainer.container[specfile][featureId].rt = corrRt


##TODO: Code is deprecated, new classes are currently located in maspy.featuregrouping
#class FeatureGroupItem(object):
#    """Representation of a group of :class:`FeatureItem`.
#
#    :ivar isMatched: False by default, True if any :class:`FeatureItem` in the group are matched.
#    :ivar isAnnotated: False by default, True if any :class:`FeatureItem` in the group are annotated.
#    :ivar siIds: containerId values of matched Si entries
#    :ivar siiIds: containerId values of matched Sii entries
#    :ivar featureIds: containerId values of :class:`FeatureItem` in the feature group
#    :ivar peptide: peptide sequence of best scoring Sii match
#    :ivar sequence: plain amino acid sequence of best scoring Sii match, used to retrieve protein information
#    :ivar score: score of best scoring Sii match
#    :ivar matchMatrix: structured representation of :attr:`FeatureItem.containerId` in the feature group.
#    :ivar intensityMatrix: similar to :attr:`matchMatrix` but contains :attr:`FeatureItem.intensity` values.
#    {chargeState: 2d numpy.array with specfiles as 1st dimension and labelState as 2nd dimension}
#    """
#    def __init__(self):
#        self.isMatched = None
#        self.isAnnotated = None
#        self.siIds = list()
#        self.siiIds = list()
#        self.featureIds = list()
#        self.peptide = None
#        self.sequence = None
#        self.score = None
#        self.matchMatrix = dict()
#        self.intensityMatrix = dict()
#
#
#class FeatureGroupContainer(object):
#    """ItemContainer for peptide feature groups :class`FeatureGroupItem`.
#
#    :ivar container: Storage list of :class:`FeatureGroupItem`
#    :ivar index: Use :attr:`FeatureItem.containerId` to which :class:`FeatureGroupItem` the feature was grouped
#    :ivar labelDescriptor: :class:`maspy.sil.LabelDescriptor` describes the label setup of an experiment
#    :ivar specfiles: List of keywords (filenames) representing files
#    :ivar specfilePositions: {specfile:arrayPosition, ...}
#    arrayPosition respresents the array position of a specfile in :attr:`FeatureGroupItem.matchMatrix`
#    """
#    def __init__(self, specfiles, labelDescriptor=None):
#        self.container = dict()
#        self.labelDescriptor = maspy.sil.LabelDescriptor() if labelDescriptor is None else labelDescriptor
#        self._index = 0
#
#        self.info = dict()
#        for position, specfile in enumerate(specfiles):
#            self.info[specfile] = {'matrixPosition': position}
#
#    def getItems(self, specfiles=None, sort=False, reverse=False, selector=lambda fgi: True):
#        """Generator that yields filtered and/or sorted :class:`Si` objects from :instance:`self.sic`
#
#        :param specfiles: filenames of msrun files - if specified return only items from those files
#        :type specfiles: str or [str, str, ...]
#        :param sort: if "sort" is specified the returned list of items is sorted according to the :class:`Si`
#        attribute specified by "sort", if the attribute is not present the item is skipped.
#        :param reverse: boolean to reverse sort order
#        :param selector: a function which is called with each :class:`Si` item and returns
#        True (include item) or False (discard item). If not specified all items are returned
#        """
#        specfiles = [_ for _ in viewkeys(self.info)] if specfiles is None else aux.toList(specfiles)
#        return _getItems(self.container, specfiles, sort, reverse, selector)
#
#    def getArrays(self, report='lfq', attr=None, specfiles=None, sort=False, reverse=False, selector=lambda si: True, defaultValue=None):
#        """Return a condensed array of data selected from :class:`Si` objects of :instance:`self.sic`
#        for fast and convenient data processing.
#
#        :param attr: list of :class:`Si` item attributes that should be added to the returned array.
#        If an attribute is not present the "defaultValue" is added instead. The attributes "id" and "specfile"
#        are always included, in combination they serve as a unique id.
#        :param specfiles: filenames of msrun files - if specified return only items from those files
#        :type specfiles: str or [str, str, ...]
#        :param sort: if "sort" is specified the returned list of items is sorted according to the :class:`Si`
#        attribute specified by "sort", if the attribute is not present the item is skipped.
#        :param reverse: boolean to reverse sort order
#        :param selector: a function which is called with each :class:`Si` item and returns
#        True (include item) or False (discard item). If not specified all items are returned
#
#        return {'attribute1': numpy.array(), 'attribute1': numpy.array(), ...}
#        """
#        attr = attr if attr is not None else []
#        attr = set(['id', 'specfile'] + aux.toList(attr))
#        specfiles = [_ for _ in viewkeys(self.info)] if specfiles is None else aux.toList(specfiles)
#
#        arrays = arrays = dict([(key, []) for key in attr])
#        reportAttributes = list()
#        if report == 'lfq':
#            arrays['charge'] = list()
#            arrays['labelState'] = list()
#            for specfile in self.specfiles:
#                arrays[specfile] = list()
#                reportAttributes.append(specfile)
#        elif report == 'sil':
#            arrays['charge'] = list()
#            arrays['specfile'] = list()
#            for labelState in list(viewkeys(self.labelDescriptor.labels)) + [-1]:
#                labelAttributeName = ' '.join(('label:', str(labelState)))
#                arrays[labelAttributeName] = list()
#                reportAttributes.append(labelAttributeName)
#
#        if report == 'sil':
#            for item in _getItems(self.container, specfiles, sort, reverse, selector):
#                for charge in viewkeys(item.intensityMatrix):
#                    for specfile in specfiles:
#                        specfilePosition = self.info[specfile]['matrixPosition']
#                        for key in attributes:
#                            arrays[key].append(getattr(item, key, None))
#                        arrays['charge'].append(charge)
#                        arrays['specfile'].append(specfile)
#                        for labelState in list(viewkeys(self.labelDescriptor.labels)) + [-1]:
#                            labelAttributeName = ' '.join(('label:', str(labelState)))
#                            arrays[labelAttributeName].append(item.intensityMatrix[charge][specfilePosition, labelState])
#        elif report == 'lfq':
#            for item in _getItems(self.container, specfiles, sort, reverse, selector):
#                for charge in viewkeys(item.intensityMatrix):
#                    for labelState in list(viewkeys(self.labelDescriptor.labels)) + [-1]:
#                        for key in attributes:
#                            arrays[key].append(getattr(item, key, None))
#                        arrays['charge'].append(charge)
#                        arrays['labelState'].append(labelState)
#                        for specfile in specfiles:
#                            specfilePosition = self.info[specfile]['matrixPosition']
#                            arrays[specfile].append(item.intensityMatrix[charge][specfilePosition, labelState])
#        else:
#            raise Exception('report must be either "lfq" or "sil", not '+report)##
#
#        for key in  [_ for _ in viewkeys(arrays)]:
#            if key in reportAttributes:
#                arrays[key] = numpy.array(arrays[key], dtype=numpy.float64)
#            else:
#                arrays[key] = numpy.array(arrays[key])
#        return arrays