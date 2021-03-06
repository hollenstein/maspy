"""
Work in progress: This module contains functions to perform feature grouping
between multiple ms-runs/specfiles and the class Fgi and FgiContainer to
represent and store such groups.
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

import io
import json
import zipfile

import numpy
import scipy

import maspy.auxiliary as aux
from maspy.core import _getItems
from maspy.core import _getArrays


# --- Fgi (feature group item) and FgiContainer classes --- #
class Fgi(object):
    """ #TODO: docstring
    """
    def __init__(self, identifier):
        self.id = identifier
        self.isValid = None
        self.charge = None
        self.isAnnotated = None

        self.specfiles = None
        self.featureIds = None

        self.clusterType = None
        self.intensities = None
        self.coordinates = None

    def _reprJSON(self):
        """Returns a JSON serializable represenation of a ``Fgi`` class
        instance. Use :func:`maspy.core.Fgi._fromJSON()` to generate a new
        ``Fgi`` instance from the return value.

        :returns: a JSON serializable python object
        """
        return {'__Fgi__': self.__dict__}

    @classmethod
    def _fromJSON(cls, jsonobject):
        """Generates a new instance of :class:`maspy.core.Fgi` from a decoded
        JSON object (as generated by :func:`maspy.core.Fgi._reprJSON()`).

        :param jsonobject: decoded JSON object

        :returns: a new instance of :class:`Fgi`
        """
        newInstance = cls(None)
        newInstance.__dict__.update(jsonobject)
        #newInstance.featureIds = numpy.array(newInstance.featureIds)
        #newInstance.specfiles = numpy.array(newInstance.specfiles)

        return newInstance

    @staticmethod
    def jsonHook(encoded):
        """Custom JSON decoder that allows construction of a new ``Fgi``
        instance from a decoded JSON object.

        :param encoded: a JSON decoded object literal (a dict)

        :returns: "encoded" or :class:`Fgi`
        """
        if '__Fgi__' in encoded:
            return Fgi._fromJSON(encoded['__Fgi__'])
        else:
            return encoded


class FgiContainer(object):
    """ #TODO: docstring
    """
    def __init__(self, specfiles):
        self.container = dict()
        self.info = dict()
        self._nextFgiId = 0
        self._matrixTemplate = list()
        for position, specfile in enumerate(specfiles):
            self.info[specfile] = {'matrixPosition': position}
            self._matrixTemplate.append(specfile)

    def _getNextFgiId(self):
        fgiId = str(self._nextFgiId)
        self._nextFgiId += 1
        return fgiId

    def getItem(self, identifier):
        """ #TODO: docstring
        """
        return self.container[identifier]

    def getArrays(self, attr=None, sort=False, reverse=False,
                  selector=None, defaultValue=None, report='lfq'):
        """ #TODO: docstring
        """
        selector = (lambda fgi: fgi.isValid) if selector is None else selector
        attr = attr if attr is not None else []
        attr = set(['id', 'intensities'] + aux.toList(attr))
        items = self.getItems(sort, reverse, selector)
        arrays = _getArrays(items, attr, defaultValue)

        for specfile in self._matrixTemplate:
            arrays[specfile] = list()
        for intensities in arrays['intensities']:
            for specfile, intensitiy in zip(self._matrixTemplate, intensities):
                arrays[specfile].append(intensitiy)
        for specfile in self._matrixTemplate:
            arrays[specfile] = numpy.array(arrays[specfile],
                                           dtype=numpy.float64
                                           )
        del arrays['intensities']

        return arrays

    def getItems(self, sort=False, reverse=False,  selector=None):
        """ #TODO: docstring
        """
        selector = (lambda fgi: fgi.isValid) if selector is None else selector
        _container = {'_': self.container}
        return _getItems(_container, '_', sort, reverse, selector)

    def save(self, path, name, compress=True):
        """Writes the specified specfiles to ``fgic`` files on the hard disk.

        .. note::
            If ``.save()`` is called and no ``fgic`` files are present in the
            specified path new files are generated, otherwise old files are
            replaced.

        :param path: filedirectory to which the ``fgic`` file is written.
        :param name: filename, without file extension
        :param compress: bool, True to use zip file compression
        """
        with aux.PartiallySafeReplace() as msr:
            filename = name + '.fgic'
            filepath = aux.joinpath(path, filename)
            with msr.open(filepath) as openfile:
                self._writeContainer(openfile, compress)

    def _writeContainer(self, filelike, compress):
        """Writes the ``self.container`` entry to the ``fgic`` format.

        :param filelike: path to a file (str) or a file-like object
        :param compress: bool, True to use zip file compression

        .. note::
            In addition it could also dump the ``self.info`` entry to the
            zipfile with the filename ``info``, but this is not used at the
            moment. For details see :func:`maspy.auxiliary.writeJsonZipfile()`
        """
        aux.writeJsonZipfile(filelike, self.container, compress=compress)
        zipcomp = zipfile.ZIP_DEFLATED if compress else zipfile.ZIP_STORED
        with zipfile.ZipFile(filelike, 'a', allowZip64=True) as containerFile:
            infodata = self.info
            infodata['_matrixTemplate'] = self._matrixTemplate
            containerFile.writestr('info', json.dumps(infodata, zipcomp))

    def load(self, path, name):
        """Imports the specified ``fgic`` file from the hard disk.

        :param path: filedirectory to which the ``fgic`` file is written.
        :param name: filename, without file extension
        """

        filename = name + '.fgic'
        filepath = aux.joinpath(path, filename)
        with zipfile.ZipFile(filepath, 'r') as containerZip:
            #Convert the zipfile data into a str object, necessary since
            #containerZip.read() returns a bytes object.
            jsonString = io.TextIOWrapper(containerZip.open('data'),
                                          encoding='utf-8'
                                          ).read()
            infoString = io.TextIOWrapper(containerZip.open('info'),
                                          encoding='utf-8'
                                          ).read()
        self.container = json.loads(jsonString, object_hook=Fgi.jsonHook)
        self.info.update(json.loads(infoString))
        self._matrixTemplate = self.info['_matrixTemplate']
        del self.info['_matrixTemplate']

    def updateIntensities(self, fiContainer, iKey='intensity'):
        """ #TODO: docstring
        :param fiContainer:
        :param iKey: Attribute name of :class:`Fi` that contains the feature
            intensity or an abundance measure. Default "intensity"
        """
        for fgi in listvalues(self.container):
            intensities = list()
            specfileIds = {i: j for i, j in zip(fgi.specfiles, fgi.featureIds)}
            for specfile in self._matrixTemplate:
                if specfile in specfileIds:
                    fi = fiContainer.getItem(specfile, specfileIds[specfile])
                    intensities.append(getattr(fi, iKey))
                else:
                    intensities.append(None)
            fgi.intensities = intensities


def updateFgiAnnotationFromFi(fgiContainer, fiContainer, largerBetter):
    """ #TODO: docstring

    :param fgiContainer:
    :param fiContainer:
    :param largerBetter:
    """
    for fgi in listvalues(fgiContainer.container):
        annotations = list()
        for specfile, fiId in zip(fgi.specfiles, fgi.featureIds):
            fi = fiContainer.getItem(specfile, fiId)
            if not fi.isAnnotated:
                continue
            annotations.append([fi.score, fi.peptide, fi.sequence])
        annotations.sort(reverse=largerBetter)
        if len(annotations) > 0:
            fgi.isAnnotated = True
            fgi.score = annotations[0][0]
            fgi.peptide = annotations[0][1]
            fgi.sequence = annotations[0][2]
        else:
            fgi.isAnnotated = False


# --- Functions to do a continuity grouping on feature arrays --- #
def continuityGrouping(values, limit):
    """ #TODO docstring

    :param values: ``numpy.array`` containg ``int`` or ``float``, must be sorted
    :param limit: the maximal difference between two values, if this number is
        exceeded a new group is generated

    :returns: a list containing array start and end positions of continuous
        groups
    """
    lastValue = values[0]
    lastPos = 0
    groupStartPos = 0

    groupPos = list()
    for currPos, currValue in enumerate(values):
        if currValue - lastValue > limit:
            groupPos.append((groupStartPos, lastPos))
            groupStartPos = currPos
        lastPos = currPos
        lastValue = currValue
    groupPos.append((groupStartPos, lastPos))
    return groupPos


def massTimeContinuityGroups(arrays, mKey, tKey, mLimit, tLimit):
    """ #TODO docstring

    :param arrays: a dictionary containing ``numpy.arrays``, must be sorted
        according to the "mKey" (mass key) value.
    :param mKey: "arrays" key that contains the mass ``numpy.array``
    :param tKey: "arrays" key that contains the time ``numpy.array``
    :param mLimit: maximal mass difference for separating continuity groups
    :param tLimit: maximal time difference for separating continuity groups

    :returns: a list containing array positions of continuous groups."""
    arrayPositions = numpy.array(range(listvalues(arrays)[0].size))

    finalGroupPositions = list()
    for start, end in continuityGrouping(arrays[mKey], mLimit):
        if start == end:
            finalGroupPositions.append(arrayPositions[start:end+1])
            continue

        #Perform time continuity grouping on the mass continuity groups
        preSelectionT = arrays[tKey][start:end+1]
        preSelectionM = arrays[mKey][start:end+1]
        preSelectionPositions = arrayPositions[start:end+1]
        _sort = numpy.argsort(preSelectionT)
        preGroups = continuityGrouping(preSelectionT[_sort], tLimit)

        #Perform a second round of mass continuity grouping
        finalGroupPrePos = list()
        for _start, _end in preGroups:
            preGroupPos = sorted(_sort[_start:_end+1])
            secGroups = continuityGrouping(preSelectionM[preGroupPos], mLimit)
            for fStart, fEnd in secGroups:
                finalGroupPrePos.append(preGroupPos[fStart:fEnd+1])

        #Add the final group positions
        for _pos in finalGroupPrePos:
            finalGroupPositions.append(preSelectionPositions[_pos])

    return finalGroupPositions


def getContGroupArrays(arrays, groupPositions, arrayKeys=None):
    """Convinience function to generate a subset of arrays from specified array
    positions.

    :param arrays: a dictionary containing ``numpy.arrays``
    :param groupPositions: arrays positions that should be included in the
        subset of arrays
    :param arrayKeys: a list of "arrays" keys that should be included in the
        subset of arrays, if None all keys are selected

    :returns: a dictionary containing ``numpy.arrays``
    """
    if arrayKeys is None:
        arrayKeys = list(viewkeys(arrays))
    matchingArrays = dict()
    for key in arrayKeys:
        matchingArrays[key] = arrays[key][groupPositions]
    return matchingArrays


def log2RelativeMassLimit(relativeMassLimit):
    """Converts a relative mass limit to a corresponding value, which can be
    used as an absolute limit for log2 transformed mass values.

    :param relativeMassLimit: float, for example the commonly used value 10ppm
        (10*1e-6)

    :returns: numpy.log2(1 + relativeMassLimit)
    """
    return numpy.log2(1 + relativeMassLimit)


# --- Performing a proximity clustering of featrues --- #
def calcDistMatchArr(matchArr, tKey, mKey):
    """Calculate the euclidean distance of all array positions in "matchArr".

    :param matchArr: a dictionary of ``numpy.arrays`` containing at least two
        entries that are treated as cartesian coordinates.
    :param tKey: #TODO: docstring
    :param mKey: #TODO: docstring

    :returns: #TODO: docstring

            {'eucDist': numpy.array([eucDistance, eucDistance, ...]),
             'posPairs': numpy.array([[pos1, pos2], [pos1, pos2], ...])
             }
    """
    #Calculate all sorted list of all eucledian feature distances
    matchArrSize = listvalues(matchArr)[0].size

    distInfo = {'posPairs': list(), 'eucDist': list()}
    _matrix = numpy.swapaxes(numpy.array([matchArr[tKey], matchArr[mKey]]), 0, 1)

    for pos1 in range(matchArrSize-1):
        for pos2 in range(pos1+1, matchArrSize):
            distInfo['posPairs'].append((pos1, pos2))
    distInfo['posPairs'] = numpy.array(distInfo['posPairs'])
    distInfo['eucDist'] = scipy.spatial.distance.pdist(_matrix)

    distSort = numpy.argsort(distInfo['eucDist'])
    for key in list(viewkeys(distInfo)):
        distInfo[key] = distInfo[key][distSort]

    return distInfo


def proximityGrouping(matchArr, distInfo, distLimit, categoryKey):
    """ #TODO: docstring. Group according to the distance value provided by
        ``distInfo['eucDist']`` with the limitation that each ... category value
        can occur only once per group.

    :param matchArr: #TODO: docstring
    :param distInfo: #TODO: docstring, must be sorted, provide keys "posPairs"
        and "eucDist". As generated by :func:`calcDistMatchArr()`
    :param distLimit: #TODO: docstring
    :param categoryKey: #TODO: docstring

    :returns: #TODO: docstring
    """
    #Group fi according to their proximity
    matchArrSize = listvalues(matchArr)[0].size

    linkageGroups = {p: [p] for p in range(matchArrSize)}
    posToGroup = {p: p for p in range(matchArrSize)}
    groupCategories = {p: set([s]) for p, s in zip(range(matchArrSize),
                                                  matchArr[categoryKey]
                                                  )
                      }
    for (pos1, pos2), dist in zip(distInfo['posPairs'], distInfo['eucDist']):
        if dist > distLimit:
            break

        id1 = posToGroup[pos1]
        id2 = posToGroup[pos2]
        if groupCategories[id1].intersection(groupCategories[id2]):
            continue

        linkageGroups[id1].extend(linkageGroups[id2])
        groupCategories[id1].update(groupCategories[id2])
        for _pos in linkageGroups[id2]:
            posToGroup[_pos] = id1
        del linkageGroups[id2]
        del groupCategories[id2]

    return linkageGroups


def fiGroupFromLinkageGroup(matchArr, arrPos, groupId, timeKey, massKey):
    """ #TODO: docstring
    """
    fgi = Fgi(groupId)
    matchArr['isAnnotated'][arrPos]

    minT = numpy.min(matchArr[timeKey][arrPos])
    maxT = numpy.max(matchArr[timeKey][arrPos])
    minM = numpy.min(matchArr[massKey][arrPos])
    maxM = numpy.max(matchArr[massKey][arrPos])

    fgi.isValid = True
    fgi.specfiles = matchArr['specfile'][arrPos].tolist()
    fgi.featureIds = matchArr['id'][arrPos].tolist()
    fgi.isAnnotated = numpy.any(matchArr['isAnnotated'][arrPos])
    fgi.coordinates = ((minT, maxT), (minM, maxM))
    #fgi.clusterType = clusterType

    return fgi


def generateFeatureGroups(fgiContainer, linkageGroups, matchArr, timeKey,
                          massKey, logMassKey, massScalingFactor):
    """ #TODO: docstring

    :param fgiContainer:
    :param linkageGroups:

    :returns: a list of ids of the newly generated :class:`Fgi`
    """
    #Generate feature groups from the linked features
    newFgiIds = list()
    for linkageGroup in viewvalues(linkageGroups):
        fgiId = fgiContainer._getNextFgiId()
        fgi = fiGroupFromLinkageGroup(matchArr, linkageGroup, fgiId,
                                      timeKey, massKey
                                      )
        fgiContainer.container[fgiId] = fgi
        fgi.metrics = clusterMetrics(matchArr[timeKey][linkageGroup],
                                     matchArr[logMassKey][linkageGroup],
                                     massScalingFactor=massScalingFactor
                                     )
        fgi.rt = fgi.metrics['meanTime']
        fgi.mz = fgi.metrics['meanMass']
        newFgiIds.append(fgiId)
    return newFgiIds


def findFgiOverlaps(fgiContainer, fgiIds):
    """

    :param fgiContainer:
    :param fgiIds:

    :return: bool, True if any two groups overlap else False
    """
    #Use the fgi min/max values of mass and time to define overlapping groups
    overlap = False
    for id1Pos in range(len(fgiIds)-1):
        fgi1 = fgiContainer.container[fgiIds[id1Pos]]
        for id2Pos in range(id1Pos+1, len(fgiIds)):
            fgi2 = fgiContainer.container[fgiIds[id2Pos]]

            c1 = fgi1.coordinates
            c2 = fgi2.coordinates
            if rangeOverlap(c1[0][0], c1[0][1], c2[0][0], c2[0][1]) and\
                    rangeOverlap(c1[1][0], c1[1][1], c2[1][0], c2[1][1]):
                fgi1.isValid = False
                fgi1.clusterType = 'tangled'
                fgi2.isValid = False
                fgi2.clusterType = 'tangled'
                overlap = True
    return overlap


def rangeOverlap(a_min, a_max, b_min, b_max):
    """Neither range is completely greater than the other
    #TODO docstring
    """
    return (a_min <= b_max) and (b_min <= a_max)


def clusterMetrics(timeValues, massValues, massScalingFactor=1):
    """ #TODO: docstring
    """
    metrics = dict()
    metrics['meanTime'] = numpy.mean(timeValues)
    metrics['meanMass'] = numpy.mean(massValues)
    metrics['devTime'] = timeValues - metrics['meanTime']
    metrics['devMass'] = massValues - metrics['meanMass']
    #metrics['devMass'] = (1-metrics['meanMass']/massValues)
    metrics['spreadTime'] = numpy.max(timeValues) - numpy.min(timeValues)
    metrics['spreadMass'] = numpy.max(massValues) - numpy.min(massValues)
    #metrics['spreadMass'] = (1-numpy.min(massValues) / numpy.max(massValues))
    metrics['devEuc'] = numpy.sqrt(numpy.power(metrics['devTime'], 2) +
                                   numpy.power(metrics['devMass']*massScalingFactor, 2)
                                   )
    metrics['meanEuc'] = numpy.mean(metrics['devEuc'])
    metrics['devTime'] = metrics['devTime'].tolist()
    metrics['devMass'] = metrics['devMass'].tolist()
    metrics['devEuc'] = metrics['devEuc'].tolist()
    return metrics


def lfqFeatureGrouping(fiContainer, timeLimit=40, massLimit=10*1e-6,
                       eucLimit=None, timeKey='rt', massKey='mz',
                       massScalingFactor=None, categoryKey='specfile',
                       charges=None, matchArraySelector=None, specfiles=None):
    """ #TODO: docstring

    :param fiContainer: #TODO: docstring
    :param timeLimit: #TODO: docstring
    :param massLimit: #TODO: docstring
    :param eucLimit: #TODO: docstring
    :param timeKey: #TODO: docstring
    :param massKey: #TODO: docstring
    :param massScalingFactor: #TODO: docstring
    :param categoryKey: #TODO: docstring
    :param charges: #TODO: docstring
    :param matchArraySelector: #TODO: docstring
    :param specfiles: limit grouping to these specfiles

    :returns: #TODO docstring, :class:`FgiContainer`
    """
    # --- perform the whole feature grouping process --- #
    targetChargeStates = range(1, 6) if charges is None else charges
    if matchArraySelector is None:
        matchArraySelector = lambda arr: numpy.any(arr['isAnnotated'])
    if massScalingFactor is None:
        massScalingFactor = timeLimit / massLimit
    if eucLimit is None:
        eucLimit = timeLimit
    if specfiles is None:
        specfiles = sorted(viewkeys(fiContainer.info))

    #'massToleranceMode': 'relative'
    #'timeToleranceMode': 'absolute'
    fgiContainer = FgiContainer(specfiles)

    logMassLimit = log2RelativeMassLimit(massLimit)
    logMassKey = 'logMass'

    logToleranceFactor = massLimit / log2RelativeMassLimit(massLimit)
    logMassScalingFactor = massScalingFactor * logToleranceFactor
    """ Note: because "a" is similar to "b"
    a = (1- 400 / 400.001) * massScalingFactor
    b = (numpy.log2(400.001) - numpy.log2(400)) * logMassScalingFactor
    """

    fiArrayKeys = [massKey, timeKey, 'isAnnotated', 'isMatched']

    for _charge in targetChargeStates:
        # - Prepare feature arrays - #
        fiSelector = lambda fi: fi.charge == _charge and fi.isValid
        fiArrays = fiContainer.getArrays(fiArrayKeys, specfiles, sort=massKey,
                                         selector=fiSelector)
        fiArrays['logMass'] = numpy.log2(fiArrays[massKey])
        if listvalues(fiArrays)[0].size == 0:
            continue

        # - group features which are in close mass and time proximity - #
        continuousGroups = massTimeContinuityGroups(fiArrays, logMassKey,
                                                    timeKey, logMassLimit,
                                                    timeLimit
                                                    )

        # - perform proximity grouping - #
        matchArrayKeys = list(viewkeys(fiArrays))

        for groupId in range(len(continuousGroups)):
            #Grab the arrays of the current feature continuity group
            groupPositions = continuousGroups[groupId]
            matchArr = getContGroupArrays(fiArrays, groupPositions,
                                          matchArrayKeys
                                          )
            if not matchArraySelector(matchArr):
                continue

            #Calculate a sorted list of all euclidean feature distances
            matchArr['mNorm'] = matchArr[logMassKey] * logMassScalingFactor
            distInfo = calcDistMatchArr(matchArr, timeKey, 'mNorm')

            #Group fi according to their proximity
            linkageGroups = proximityGrouping(matchArr, distInfo, eucLimit,
                                              categoryKey
                                              )

            #Generate feature groups from the linked features
            fgiIds = generateFeatureGroups(fgiContainer, linkageGroups,
                                           matchArr, timeKey, massKey,
                                           logMassKey, logMassScalingFactor
                                           )

            #Set charge manually
            for fgiId in fgiIds:
                fgiContainer.container[fgiId].charge = _charge

            #Mark overlapping groups as not valid (fgi.isValid = False)
            fgiDoOverlap = findFgiOverlaps(fgiContainer, fgiIds)

        #Add feature intensities to the feature groups
        fgiContainer.updateIntensities(fiContainer)
    return fgiContainer
