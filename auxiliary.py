from __future__ import print_function

from collections import defaultdict as ddict
import functools
import math
import operator
import os
import random

import numpy
from scipy.interpolate import LSQUnivariateSpline

from maspy.mit_stats import runningMode as runningMode
from maspy.mit_stats import runningMean as runningMean
from maspy.mit_stats import runningMedian as runningMedian


def lazyAttribute(fn):
    attributeName = '_lazy_' + fn.__name__
    @property
    def _lazyAttribute(self):
        if not hasattr(self, attributeName):
            setattr(self, attributeName, fn(self))
        return getattr(self, attributeName)
    return _lazyAttribute


def searchFileLocation(targetFileName, targetFileExtension, rootDirectory, recursive=True):
    """Search for file with specified file extension in all subfolders of specified rootDirectory, returns first matching instance.

    :type targetFileName: str
    :type rootDirectory: str
    :type targetFileExtension: str
    :ivar recursive: bool, specify whether subdirectories should be searched
    """
    expectedFileName = targetFileName.split('.')[0] + '.' + targetFileExtension
    targetFilePath = None

    if recursive:
        for dirpath, dirnames, filenames in os.walk(rootDirectory):
            for filename in filenames:
                if filename == expectedFileName:
                    targetFilePath = joinpath(dirpath, filename)
                    break
            if targetFilePath is not None:
                break
    else:
        for filename in os.listdir(rootDirectory):
            filePath = joinpath(rootDirectory, filename)
            if not os.path.isfile(filePath):
                continue
            if filename == expectedFileName:
                targetFilePath = filePath
                break

    return targetFilePath


def matchingFilePaths(targetfilename, directory, targetFileExtension=None, selector=None):
    """Search for files in all subfolders of specified directory, return filepaths of all matching instances.
    :param targetfilename: filename to search for, only the string before the last '.' is used for filename matching
    :param directory: search directory, including all subdirectories
    :param targetFileExtension: string after the last '.' in the filename, has to be identical if specified.
                                '.' in targetFileExtension are ignored, thus '.txt' is equal to 'txt'.

    :param selector: a function which is called with the value of targetfilename
    and has to return True (include value) or False (discard value).
    If no selector is specified, equality to targetfilename is used.
    """
    targetFilePaths = list()

    targetfilename = os.path.splitext(targetfilename)[0]
    targetFileExtension = targetFileExtension.replace('.', '')
    matchExtensions = False if targetFileExtension is None else True
    if selector is None:
        selector = functools.partial(operator.eq, targetfilename)

    for dirpath, dirnames, filenames in os.walk(directory):
        for filename in filenames:
            filenameNoextension = os.path.splitext(filename)[0]

            if selector(filenameNoextension):
                if matchExtensions:
                    if not filename.endswith('.' + targetFileExtension):
                        continue
                targetFilePaths.append(joinpath(dirpath, filename))
    return targetFilePaths


def findAllSubstrings(a_str, sub):
    start = 0
    while True:
        start = a_str.find(sub, start)
        if start == -1: return
        yield start
        start +=  1# use start += 1 to find overlapping matches, otherwise len(sub)


def toList(variable, types=(str, int, float)):
    """ Converts variable of type string, int, float to a list: variable -> [variable]

    :type variable: (str, int, float, others)
    """
    if isinstance(variable, types):
        return [variable]
    else:
        return variable


def joinpath(path, *paths):
    """Join two or more pathname components, inserting "/" as needed and replacing all "\\" by "/". """
    return os.path.join(path, *paths).replace('\\','/')


class Factorial():
    def __getitem__(self, n):
        try:
            return getattr(self, str(n))
        except AttributeError:
            setattr(self, str(n), math.factorial(int(n)))
            return getattr(self, str(n))


class DataFit(object):
    def __init__(self, dependentVarInput, independentVarInput):
        self.dependentVarInput = dependentVarInput
        self.independentVarInput = independentVarInput

        self.splines = None
        self.splineCycles = 10
        self.splineMinKnotPoins = 10
        self.splineOrder = 2
        self.splineInitialKnots = 200
        self.splineSubsetPercentage = 0.4
        self.splineTerminalExpansion = 0.1

        self.dependentVar = list()
        self.independentVar = list()

    def processInput(self, dataAveraging=False, windowSize=None):
        self.dependentVar = numpy.array(self.dependentVarInput, dtype=numpy.float64)
        self.independentVar = numpy.array(self.independentVarInput, dtype=numpy.float64)

        sortMask = self.independentVar.argsort()
        self.dependentVar = self.dependentVar[sortMask]
        self.independentVar = self.independentVar[sortMask]

        if dataAveraging:
            averagedData = averagingData(self.dependentVar, windowSize=windowSize, averagingType=dataAveraging)
            averagedData = numpy.array(averagedData, dtype=numpy.float64)

            missingNumHigh = numpy.floor((self.independentVar.size - averagedData.size) / 2)
            missingNumLow = (self.independentVar.size - averagedData.size) - missingNumHigh

            self.dependentVar = averagedData
            self.independentVar = self.independentVar[missingNumLow:-missingNumHigh]

    def generateSplines(self):
        self.splines = returnSplineList(self.dependentVar, self.independentVar, subsetPercentage=self.splineSubsetPercentage,
                                        cycles=self.splineCycles, minKnotPoints=self.splineMinKnotPoins, initialKnots=self.splineInitialKnots,
                                        splineOrder=self.splineOrder, terminalExpansion=self.splineTerminalExpansion
                                        )

    def __getitem__(self, value):
        returnValue = numpy.mean([numpy.nan_to_num(currSpline(value)) for currSpline in self.splines])
        return returnValue

    def corrArray(self, inputArray):
        outputArray = numpy.vstack([numpy.nan_to_num(currSpline(inputArray)) for currSpline in self.splines]).mean(axis=0)
        return outputArray


def averagingData(array, windowSize=None, averagingType='median'):
    assert averagingType in ['median', 'mean']
    if averagingType == 'median':
        if windowSize is None:
            windowSize = int(float(len(array)) / (50)) if int(float(len(array) / (50))) > 100 else 100
        averagedData = runningMedian(array, windowSize)
    elif averagingType == 'mean':
        if windowSize is None:
            windowSize = int(float(len(array)) / (50)) if int(float(len(array) / (50))) > 100 else 100
        averagedData = runningMean(array, None, windowSize)
    return averagedData


def returnSplineList(dependentVar, independentVar, subsetPercentage=0.4, cycles=10,
                     minKnotPoints=10, initialKnots=200, splineOrder=2, terminalExpansion=0.1
                     ):
    """Expects sorted arrays.

    :ivar terminalExpansion: expand subsets on both sides

    return TODO: what is returned?
    """
    expansions = ddict(list)
    expansionArea = (independentVar[-1] - independentVar[0]) * terminalExpansion
    #adds 100 data points at both ends of the dependent and independent array
    for i in range(100):
        expansions['indUp'].append(independentVar[-1] + expansionArea/100*i)
        expansions['indDown'].append(independentVar[0] - expansionArea/100*(100-i+1))
        expansions['depUp'].append(dependentVar[-1])
        expansions['depDown'].append(dependentVar[0])

    dependentVar = numpy.array(expansions['depDown'] + list(dependentVar) + expansions['depUp'], dtype=numpy.float64)
    independentVar = numpy.array(expansions['indDown'] + list(independentVar) + expansions['indUp'], dtype=numpy.float64)

    splineList = list()
    for cycle in range(cycles):
        subset = sorted(random.sample(xrange(len(dependentVar)), int(len(dependentVar)*subsetPercentage )) )
        terminalExpansion

        dependentSubset = dependentVar[subset]
        independentSubset = independentVar[subset]

        minIndVar = independentSubset[minKnotPoints]
        maxIndVar = independentSubset[-minKnotPoints]

        knots = [float(i) * (maxIndVar-minIndVar) / initialKnots + minIndVar for i in range(1, initialKnots)]
        ## remove knots with less then minKnotPoints data points  ##
        lastKnot = knots[0]
        newKnotList = [lastKnot]
        for knotPos in range(1,len(knots)):
            nextKnot = knots[knotPos]
            numHits = len(independentSubset[(independentSubset >= lastKnot) & (independentSubset <= nextKnot)] )
            if numHits >= minKnotPoints:
                newKnotList.append( nextKnot )
                lastKnot = nextKnot
        knots = newKnotList

        spline = LSQUnivariateSpline(independentSubset, dependentSubset, knots, k=splineOrder)
        splineList.append(spline)
    return splineList
