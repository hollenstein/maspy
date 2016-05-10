from __future__ import print_function, division, unicode_literals
from future.utils import viewkeys, viewvalues, viewitems, listvalues, listitems

try: # python 2.7
    from itertools import izip as zip
except ImportError: # python 3.x series
    pass
###############################################################################

try:
   basestring
except NameError:
   basestring=str #necessary for type comparison

from collections import defaultdict as ddict
import contextlib
import functools
import io
import json
import math
import operator
import os
import random
import re
import shutil
import tempfile
import zipfile

import numpy
from scipy.interpolate import LSQUnivariateSpline

from maspy.mit_stats import runningMode as runningMode
from maspy.mit_stats import runningMean as runningMean
from maspy.mit_stats import runningMedian as runningMedian


class PartiallySafeReplace(object):
    #TODO: Docstring missing
    def __init__(self):
        pass

    def __enter__(self):
        self._files = {}
        return self

    @contextlib.contextmanager
    def open(self, filepath, mode='w+b'):
        #TODO: check if file already in self._files
        #Check if the filepath can be accessed and is writable before creating the tempfile
        if not _isFileAccessible(filepath):
            raise IOError('File %s is not writtable' %(filepath, ))

        tempfilepath = None
        with tempfile.NamedTemporaryFile(delete=False, mode=mode) as tmpf:
            tempfilepath = tmpf.name
            yield tmpf
        self._files[filepath] = tempfilepath

    def __exit__(self, x, y, z):
        #Check if all filepaths can be accessed and are writable before moving the tempfiles
        for filepath in self._files:
            if not _isFileAccessible(filepath):
                raise IOError('File %s is not writtable' %(filepath, ))
        for filepath, tempfilepath in viewitems(self._files):
            #Note: here unhandled exceptions may still occur because of race conditions, messing things up.
            shutil.move(tempfilepath, filepath)


@contextlib.contextmanager
def openSafeReplace(filepath, mode='w+b'):
    #TODO: Docstring missing
    tempfileName = None
    #Check if the filepath can be accessed and is writable before creating the tempfile
    if not _isFileAccessible(filepath):
        raise IOError('File %s is not writtable' %(filepath, ))
    with tempfile.NamedTemporaryFile(delete=False, mode=mode) as tmpf:
        tempfileName = tmpf.name
        yield tmpf
    #Check if the filepath can be accessed and is writable before moving the tempfile
    if not _isFileAccessible(filepath):
        raise IOError('File %s is not writtable' %(filepath, ))
    #Note: here unhandled exceptions may still occur because of race conditions, messing things up.
    shutil.move(tempfileName, filepath)


def _isFileAccessible(filepath):
    """Returns True if the specified filepath is writable."""
    directory = os.path.dirname(filepath)
    if not os.access(directory, os.W_OK):
        #Return False if directory does not exist or is not writable
        return False
    if os.path.exists(filepath):
        if not os.access(filepath, os.W_OK):
            #Return False if file is not writable
            return False
        try:
            openfile = os.open(filepath, os.O_WRONLY)
            os.close(openfile)
        except IOError:
            #Return False if file is locked
            return False
    #Return True if file is writtable
    return True


# --- Json serialization helper functions --- #
class MaspyJsonEncoder(json.JSONEncoder):
    """ Extension of the json.JSONEncoder to serialze Maspy classes
    Maspy classes need to define a _reprJSON() method, which returns a json
    serializable object.

    return calls the obj._reprJSON() method if it is defined, else calls
    json.JSONEncoder.default(obj)
    """
    def default(self, obj):
        if hasattr(obj, '_reprJSON'):
            return obj._reprJSON()
        #Let the base class default method raise the TypeError
        return json.JSONEncoder.default(self, obj)


def writeJsonZipfile(filelike, data, compress=True, mode='w', name='data'):
    """ Serializes the objects contained in data to a JSON formated string and writes
    it to a zipfile.

    :param filelike: can be either a path to a file (a string) or a file-like object
    :param data: object that should be converted to a JSON formated string
    objects and types in data must be supported by the json.JSONEncoder or
    have the method _reprJSON() defined.
    :param compress: boolean, True if the zipfile should be compressed
    :param mode: 'w' to truncate and write a new file, or 'a' to append to an existing file
    :param name: the file name that will be given to the JSON output in the archive
    """
    zipcomp = zipfile.ZIP_DEFLATED if compress else zipfile.ZIP_STORED
    with zipfile.ZipFile(filelike, mode, allowZip64=True) as containerFile:
        containerFile.writestr(name, json.dumps(data, cls=MaspyJsonEncoder), zipcomp)


def writeBinaryItemContainer(filelike, binaryItemContainer, compress=True):
    """Serializes the binaryItems contained in binaryItemContainer and writes them
    into a zipfile archive.

    Examples of binaryItem classes are maspy.core.Ci() and maspy.core.Sai()

    A binaryItem class has to define the function _reprJSON() which returns
    a JSON formated string representation of the class instance. In addition
    it has to contain an attribute :attr:`arrays`, a dictionary which values
    are numpy arrays. These numpy arrays are serialized to bytes and written
    to the "binarydata" file in the archive. see :func:`_dumpArrayDictToFile`

    The JSON formated string representation of the binaryItems, together with
    the metadata, necessary to restore serialized numpy arrays, is written
    to the "metadata" file in the archive:
    [[serialized binaryItem, [metadata of a numpy array, ...]], ... ]

    Use the method :func:`loadBinaryItemContainer()` to restore a
    binaryItemContainer from a zipfile.

    :param filelike: can be either a path to a file (a string) or a file-like object
    :param binaryItemContainer: a dictionary containing binaryItems.
    """
    allMetadata = dict()
    binarydatafile = io.BytesIO()
    #Note: It would be possible to sort the items here
    for index, binaryItem in enumerate(viewvalues(binaryItemContainer)):
        metadataList = _dumpArrayDictToFile(binarydatafile, binaryItem.arrays)
        allMetadata[index] = [binaryItem._reprJSON(), metadataList]

    #TODO: Is seek here still necessary?
    binarydatafile.seek(0)

    zipcomp = zipfile.ZIP_DEFLATED if compress else zipfile.ZIP_STORED
    with zipfile.ZipFile(filelike, 'w', allowZip64=True) as containerFile:
        containerFile.writestr('metadata', json.dumps(allMetadata, cls=MaspyJsonEncoder), zipcomp)
        containerFile.writestr('binarydata', binarydatafile.getvalue(), zipcomp)


def _dumpArrayDictToFile(filelike, arrayDict):
    """Function to serialize and write numpy arrays contained in a dictionary
    to a file.

    :param filelike: can be a file or a file-like object that provide the methods
    write() and tell()
    :param arrayDict: a dictionary which values are numpy arrays, these numpy
    arrays are serialized to bytes and written to the filelike.

    returns a list of metadata dictionaries, a metadata dictionary contains
    information necessary to restore the array from the file and the corresponding
    key from the arrayDict as 'arrayKey'.

    see :func:`_dumpArrayToFile` and :func:`_dumpNdarrayToFile`
    """
    metadataList = list()
    for arrayKey in sorted(arrayDict):
        array = arrayDict[arrayKey]
        if array.ndim == 1:
            metadata = _dumpArrayToFile(filelike, array)
        else:
            metadata = _dumpNdarrayToFile(filelike, array)
        metadata['arrayKey'] = arrayKey
        metadataList.append(metadata)
    return metadataList


def _dumpArrayToFile(filelike, array):
    """Serializes a 1-dimensional numpy array to bytes, writes the bytes to the filelike object
    and returns a dictionary with metadata, necessary to restore the array from the file.

    :param filelike: can be a file or a file-like object that provide the methods
        write() and tell()
    :param array: a 1-dimensional numpy array

    returns a dictionary with metadata: {'start': start position in the file, 'end': end position in the file,
                                         'size': size of the array, 'dtype': numpy data type of the array}
    """
    bytedata = array.tobytes('C')
    start = filelike.tell()
    end = start + len(bytedata)
    metadata = {'start': start, 'end': end, 'size': array.size, 'dtype': array.dtype.name}
    filelike.write(bytedata)
    return metadata


def _dumpNdarrayToFile(filelike, ndarray):
    """Serializes a N-dimensional numpy array to bytes, writes the bytes to the filelike object
    and returns a dictionary with metadata, necessary to restore the array from the file.

    :param filelike: can be a file or a file-like object that provide the methods
        write() and tell()
    :param array: a N-dimensional numpy array

    returns a dictionary with metadata: {'start': start position in the file, 'end': end position in the file,
                                         'size': size of the array, 'dtype': numpy data type of the array,
                                         'shape': description of the array shape}
    """
    bytedata = array.tobytes('C')
    start = filelike.tell()
    end = start + len(bytedata)
    metadata = {'start': start, 'end': end, 'size': array.size, 'dtype': array.dtype.name, 'shape': array.shape}
    filelike.write(bytedata)
    return metadata

    """ Serializes the binaryItems contained in binaryItemContainer and writes them
    into a zipfile archive."""


def loadBinaryItemContainer(zippedfile, jsonHook):
    """Imports binaryItems from a zipfile generated by :func:`writeBinaryItemContainer`

    :param zipfile: can be either a path to a file (a string) or a file-like object
    :param jsonHook: a custom decoding function for JSON formated strings of the
    binaryItems stored in the zipfile.

    returns a dictionary containing binaryItems in the form: {binaryItem.id: binaryItem, ... }
    """
    binaryItemContainer = dict()
    with zipfile.ZipFile(zippedfile, 'r') as containerZip:
        #Convert the zipfile data into a str object, necessary since containerZip.read() returns a bytes object.
        metadataText = io.TextIOWrapper(containerZip.open('metadata'), encoding='utf-8').read()
        allMetadata = json.loads(metadataText, object_hook=jsonHook)
        metadataIndex = [str(_) for _ in sorted([int(i) for i in viewkeys(allMetadata)])]
        binarydataFile = containerZip.open('binarydata')
        for index in metadataIndex:
            binaryItem = allMetadata[index][0]
            for binaryMetadata in allMetadata[index][1]:
                arrayKey = binaryMetadata['arrayKey']
                rawdata = binarydataFile.read(binaryMetadata['end']-binaryMetadata['start'])
                array = _arrayFromBytes(rawdata, binaryMetadata)
                binaryItem.arrays[arrayKey] = array
            binaryItemContainer[binaryItem.id] = binaryItem
    return binaryItemContainer


def _arrayFromBytes(dataBytes, metadata):
    """Generates and returns a numpy array from raw data bytes.

    :param bytes: raw data bytes as generated by numpy.ndarray.tobytes()
    :param metadata: a dictionary containing the data type and optionally the
    shape parameter to reconstruct a numpy array from the raw data bytes.
    eg. {"dtype": "float64", "shape": (2, 3)}
    """
    array = numpy.frombuffer(dataBytes, dtype=numpy.typeDict[metadata['dtype']])
    if 'shape' in metadata:
        array = array.reshape(metadata['shape'])
    return array


# --- not yet named section --- #
def searchFileLocation(targetFileName, targetFileExtension, rootDirectory, recursive=True):
    """Search for files with specified file extension in all subfolders of specified rootDirectory, returns first matching instance.

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


def listFiletypes(targetfilename, directory):
    """Looks for all occurences of a specified filename in a directory
    and returns a list of all present file extensions of this filename.

    :param targetfilename: a filename without any extensions
    :param directory: only files present in this directory are compared
    to the targetfilename

    An extension is everything from the first dot to the end.
    In this case the filename is everything before the first dot.

    Only the first dot of the extension is omitted. eg "txt" but "txt.zip"

    :returns: a list of extensions
    """
    targetextensions = list()
    for filename in os.listdir(directory):
        if not os.path.isfile(joinpath(directory, filename)):
            continue
        splitname = filename.split('.')
        basename = splitname[0]
        extension = '.'.join(splitname[1:])
        if basename == targetfilename:
            targetextensions.append(extension)
    return targetextensions


def findAllSubstrings(string, substring):
    """ Returns a list of all substring starting positions in string
    or an empty list if substring is not present in string."""
    #TODO: solve with regex? what about '.': return [m.start() for m in re.finditer('(?='+substring+')', string)]
    start = 0
    positions = []
    while True:
        start = string.find(substring, start)
        if start == -1:
            break
        positions.append(start)
        start += 1 #to find overlapping matches
    return positions


def toList(variable, types=(basestring, int, float, )):
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


# --- not yet named section, involving all sorts of calculations --- #
class Memoize(object):
    """A general memoization class, specify a function when creating a
    new instance of the class. The functions return value is returned
    and stored in self.memo when the instance is called with an argument
    for the first time. Later calls with the same argument return the
    cached value, instead of calling the function again.
    """
    def __init__(self, function):
        self.function = function
        self.memo = {}
    def __call__(self, arg):
        if not arg in self.memo:
            self.memo[arg] = self.function(arg)
        return self.memo[arg]
    #For multiple arguments use instead:
    #def __call__(self, *args):
    #    if not args in self.memo:
    #        self.memo[args] = self.function(*args)
    #    return self.memo[args]


factorial = Memoize(lambda n: math.factorial(int(n)))
"""Returns the factorial of a number, the results of already
calculated numbers are stores in factorial.memo """


log10factorial = Memoize(lambda n: math.log10(math.factorial(int(n))))
"""Returns the log10 factorial of a number, the results of already
calculated numbers are stores in log10factorial.memo """


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
    if windowSize is None:
        windowSize = int(len(array) / 50) if int(len(array) / 50) > 100 else 100

    if averagingType == 'median':
        averagedData = runningMedian(array, windowSize)
    elif averagingType == 'mean':
        averagedData = runningMean(array, len(array), windowSize)
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
        subset = sorted(random.sample(range(len(dependentVar)), int(len(dependentVar) * subsetPercentage)))
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
            numHits = len(independentSubset[(independentSubset >= lastKnot) & (independentSubset <= nextKnot)])
            if numHits >= minKnotPoints:
                newKnotList.append(nextKnot)
                lastKnot = nextKnot
        knots = newKnotList

        spline = LSQUnivariateSpline(independentSubset, dependentSubset, knots, k=splineOrder)
        splineList.append(spline)
    return splineList


def tolerantArrayMatching(referenceArray, matchArray, matchTolerance=20, matchUnit='ppm'):
    #TODO: Docstring
    #TODO: change matchUnit to "absoulte", "relative" and remove the "*1e-6"
    """arrays must be sorted"""
    if matchUnit == 'ppm':
        lowerLimitMatchArray = matchArray * (1 - matchTolerance*1e-6)
        upperLimitMatchArray = matchArray * (1 + matchTolerance*1e-6)
    elif matchUnit == 'da':
        lowerLimitMatchArray = matchArray - matchTolerance
        upperLimitMatchArray = matchArray + matchTolerance
    else:
        raise Exception('wrong matchUnit type specified (da or ppm): ', matchUnit)

    lowerLimitMask = numpy.zeros_like(matchArray, dtype=int)
    upperLimitMask = numpy.zeros_like(matchArray, dtype=int)

    refPosLow = int()
    maxReferenceValue = referenceArray[-1]
    for matchPos, (lowerMatchEntry, upperMatchEntry) in enumerate(zip(lowerLimitMatchArray, upperLimitMatchArray)):
        if lowerMatchEntry < maxReferenceValue:
            while referenceArray[refPosLow] < lowerMatchEntry:
                refPosLow += 1
            refPosHigh = refPosLow
            #Try except statement because this case can only happen once at the end of the array
            try:
                while referenceArray[refPosHigh] <= upperMatchEntry:
                    refPosHigh += 1
            except IndexError:
                refPosHigh = len(referenceArray) - 1
            lowerLimitMask[matchPos] = refPosLow
            upperLimitMask[matchPos] = refPosHigh
        else:
            refPosLow = len(referenceArray) - 1
            refPosHigh = len(referenceArray) - 1
            lowerLimitMask[matchPos] = refPosLow
            upperLimitMask[matchPos] = refPosHigh
            break

    matchPos += 1
    lowerLimitMask[matchPos:len(matchArray)] = refPosLow
    upperLimitMask[matchPos:len(matchArray)] = refPosHigh

    return lowerLimitMask, upperLimitMask
