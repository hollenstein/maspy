"""
A collection of helper functions used in different modules. Functions deal for
example with saving files, encoding and decoding data in the json format,
filtering of numpy arrays, data fitting, etc.
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

try:
    basestring
except NameError:
    #necessary for type comparison
    basestring=str
################################################################################

import bisect
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

from maspy.mit_stats import runningMean as runningMean
from maspy.mit_stats import runningMedian as runningMedian


class PartiallySafeReplace(object):
    """Indirectly overwrite files by writing to temporary files and replacing
    them at once.

    This is a context manager. When the context is entered, subsequently opened
    files will actually open temporary files. Each time the same file-path is
    opened, the same temporary file will be used.

    When the context is closed, it will attempt to replace the original files
    with the content of the temporary files. Thus, several files can be prepared
    with less risk of data loss. Data loss is still possible if the replacement-
    operation fails (due to locking, not handled yet) or is interrupted.
    """
    def __init__(self):
        pass

    def __enter__(self):
        self._files = {}
        return self

    @contextlib.contextmanager
    def open(self, filepath, mode='w+b'):
        """Opens a file - will actually return a temporary file but replace the
        original file when the context is closed.
        """
        #Check if the filepath can be accessed and is writable before creating
        #the tempfile
        if not _isFileAccessible(filepath):
            raise IOError('File %s is not writable' % (filepath,))

        if filepath in self._files:
            with open(self._files[filepath], mode=mode) as tmpf:
                yield tmpf
        else:
            tempfilepath = None
            with tempfile.NamedTemporaryFile(delete=False, mode=mode) as tmpf:
                tempfilepath = tmpf.name
                yield tmpf
            self._files[filepath] = tempfilepath

    def __exit__(self, x, y, z):
        #Check if all filepaths can be accessed and are writable before moving
        #the tempfiles
        for filepath in self._files:
            if not _isFileAccessible(filepath):
                raise IOError('File %s is not writable' % (filepath, ))
        for filepath, tempfilepath in viewitems(self._files):
            #Note: here unhandled exceptions may still occur because of race
            #conditions, messing things up.
            shutil.move(tempfilepath, filepath)


@contextlib.contextmanager
def openSafeReplace(filepath, mode='w+b'):
    """Context manager to open a temporary file and replace the original file on
    closing.
    """
    tempfileName = None
    #Check if the filepath can be accessed and is writable before creating the
    #tempfile
    if not _isFileAccessible(filepath):
        raise IOError('File %s is not writtable' % (filepath, ))
    with tempfile.NamedTemporaryFile(delete=False, mode=mode) as tmpf:
        tempfileName = tmpf.name
        yield tmpf
    #Check if the filepath can be accessed and is writable before moving the
    #tempfile
    if not _isFileAccessible(filepath):
        raise IOError('File %s is not writtable' % (filepath, ))
    #Note: here unhandled exceptions may still occur because of race conditions,
    #messing things up.
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
    """Extension of the json.JSONEncoder to serialize MasPy classes.

    Maspy classes need to define a _reprJSON() method, which returns a json
    serializable object.
    """

    def default(self, obj):
        """
        :returns: obj._reprJSON() if it is defined, else
            json.JSONEncoder.default(obj)
        """
        if hasattr(obj, '_reprJSON'):
            return obj._reprJSON()
        #Let the base class default method raise the TypeError
        return json.JSONEncoder.default(self, obj)


def writeJsonZipfile(filelike, data, compress=True, mode='w', name='data'):
    """Serializes the objects contained in data to a JSON formated string and
    writes it to a zipfile.

    :param filelike: path to a file (str) or a file-like object
    :param data: object that should be converted to a JSON formated string.
        Objects and types in data must be supported by the json.JSONEncoder or
        have the method ``._reprJSON()`` defined.
    :param compress: bool, True to use zip file compression
    :param mode: 'w' to truncate and write a new file, or 'a' to append to an
        existing file
    :param name: the file name that will be given to the JSON output in the
        archive
    """
    zipcomp = zipfile.ZIP_DEFLATED if compress else zipfile.ZIP_STORED
    with zipfile.ZipFile(filelike, mode, allowZip64=True) as containerFile:
        containerFile.writestr(name, json.dumps(data, cls=MaspyJsonEncoder),
                               zipcomp
                               )


def writeBinaryItemContainer(filelike, binaryItemContainer, compress=True):
    """Serializes the binaryItems contained in binaryItemContainer and writes
    them into a zipfile archive.

    Examples of binaryItem classes are :class:`maspy.core.Ci` and
    :class:`maspy.core.Sai`. A binaryItem class has to define the function
    ``_reprJSON()`` which returns a JSON formated string representation of the
    class instance. In addition it has to contain an attribute ``.arrays``, a
    dictionary which values are ``numpy.array``, that are serialized to bytes
    and written to the ``binarydata`` file of the zip archive. See
    :func:`_dumpArrayDictToFile()`

    The JSON formated string representation of the binaryItems, together with
    the metadata, necessary to restore serialized numpy arrays, is written
    to the ``metadata`` file of the archive in this form:
    ``[[serialized binaryItem, [metadata of a numpy array, ...]], ...]``

    Use the method :func:`loadBinaryItemContainer()` to restore a
    binaryItemContainer from a zipfile.

    :param filelike: path to a file (str) or a file-like object
    :param binaryItemContainer: a dictionary containing binaryItems
    :param compress: bool, True to use zip file compression
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
        containerFile.writestr('metadata',
                               json.dumps(allMetadata, cls=MaspyJsonEncoder),
                               zipcomp
                               )
        containerFile.writestr('binarydata', binarydatafile.getvalue(), zipcomp)


def _dumpArrayDictToFile(filelike, arrayDict):
    """Function to serialize and write ``numpy.array`` contained in a dictionary
    to a file. See also :func:`_dumpArrayToFile` and :func:`_dumpNdarrayToFile`.


    :param filelike: can be a file or a file-like object that provides the
        methods ``.write()`` and ``.tell()``.
    :param arrayDict: a dictionary which values are ``numpy.array``, that are
        serialized to bytes and written to the filelike.

    :returns: a list of metadata dictionaries
        a metadata dictionary contains information necessary to restore the
        ``numpy.arrays`` from the file and the corresponding key from the
        arrayDict as 'arrayKey'.
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
    """Serializes a 1-dimensional ``numpy.array`` to bytes, writes the bytes to
    the filelike object and returns a dictionary with metadata, necessary to
    restore the ``numpy.array`` from the file.

    :param filelike: can be a file or a file-like object that provides the
        methods ``.write()`` and ``.tell()``.
    :param array: a 1-dimensional ``numpy.array``

    :returns: a metadata dictionary ::
        {'start': start position in the file, 'end': end position in the file,
         'size': size of the array, 'dtype': numpy data type of the array
         }
    """
    bytedata = array.tobytes('C')
    start = filelike.tell()
    end = start + len(bytedata)
    metadata = {'start': start, 'end': end, 'size': array.size,
                'dtype': array.dtype.name
                }
    filelike.write(bytedata)
    return metadata


def _dumpNdarrayToFile(filelike, ndarray):
    """Serializes an N-dimensional ``numpy.array`` to bytes, writes the bytes to
    the filelike object and returns a dictionary with metadata, necessary to
    restore the ``numpy.array`` from the file.

    :param filelike: can be a file or a file-like object that provides the
        methods ``.write()`` and ``.tell()``.
    :param ndarray: a N-dimensional ``numpy.array``

    :returns: a metadata dictionary ::
        {'start': start position in the file, 'end': end position in the file,
         'size': size of the array, 'dtype': numpy data type of the array,
         'shape': description of the array shape
         }
    """
    bytedata = ndarray.tobytes('C')
    start = filelike.tell()
    end = start + len(bytedata)
    metadata = {'start': start, 'end': end, 'size': ndarray.size,
                'dtype': ndarray.dtype.name, 'shape': ndarray.shape
                }
    filelike.write(bytedata)
    return metadata


def loadBinaryItemContainer(zippedfile, jsonHook):
    """Imports binaryItems from a zipfile generated by
    :func:`writeBinaryItemContainer`.

    :param zipfile: can be either a path to a file (a string) or a file-like
        object
    :param jsonHook: a custom decoding function for JSON formated strings of the
        binaryItems stored in the zipfile.

    :returns: a dictionary containing binaryItems
        ``{binaryItem.id: binaryItem, ... }``
    """
    binaryItemContainer = dict()
    with zipfile.ZipFile(zippedfile, 'r') as containerZip:
        #Convert the zipfile data into a str object, necessary since
        #containerZip.read() returns a bytes object.
        metadataText = io.TextIOWrapper(containerZip.open('metadata'),
                                        encoding='utf-8'
                                        ).read()
        allMetadata = json.loads(metadataText, object_hook=jsonHook)
        metadataIndex = [str(_) for _ in sorted([int(i) for i in
                                                 viewkeys(allMetadata)
                                                 ])
                         ]
        binarydataFile = containerZip.open('binarydata')
        for index in metadataIndex:
            binaryItem = allMetadata[index][0]
            for binaryMetadata in allMetadata[index][1]:
                arrayKey = binaryMetadata['arrayKey']
                rawdata = binarydataFile.read(binaryMetadata['end'] -
                                              binaryMetadata['start']
                                              )
                array = _arrayFromBytes(rawdata, binaryMetadata)
                binaryItem.arrays[arrayKey] = array
            binaryItemContainer[binaryItem.id] = binaryItem
    return binaryItemContainer


def _arrayFromBytes(dataBytes, metadata):
    """Generates and returns a numpy array from raw data bytes.

    :param bytes: raw data bytes as generated by ``numpy.ndarray.tobytes()``
    :param metadata: a dictionary containing the data type and optionally the
        shape parameter to reconstruct a ``numpy.array`` from the raw data
        bytes. ``{"dtype": "float64", "shape": (2, 3)}``

    :returns: ``numpy.array``
    """
    array = numpy.fromstring(dataBytes, dtype=numpy.typeDict[metadata['dtype']])
    if 'shape' in metadata:
        array = array.reshape(metadata['shape'])
    return array


# --- not yet named section --- #
def searchFileLocation(targetFileName, targetFileExtension, rootDirectory,
                       recursive=True):
    """Search for a filename with a specified file extension in all subfolders
    of specified rootDirectory, returns first matching instance.

    :param targetFileName: #TODO: docstring
    :type targetFileName: str
    :param rootDirectory: #TODO: docstring
    :type rootDirectory: str
    :param targetFileExtension: #TODO: docstring
    :type targetFileExtension: str
    :param recursive: bool, specify whether subdirectories should be searched

    :returns: a filepath (str) or None
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


def matchingFilePaths(targetfilename, directory, targetFileExtension=None,
                      selector=None):
    """Search for files in all subfolders of specified directory, return
    filepaths of all matching instances.

    :param targetfilename: filename to search for, only the string before the
        last "." is used for filename matching. Ignored if a selector function
        is specified.
    :param directory: search directory, including all subdirectories
    :param targetFileExtension: string after the last "." in the filename, has
        to be identical if specified. "." in targetFileExtension are ignored,
        thus ".txt" is treated equal to "txt".
    :param selector: a function which is called with the value of targetfilename
        and has to return True (include value) or False (discard value). If no
        selector is specified, equality to targetfilename is used.

    :returns: list of matching file paths (str)
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
    """Looks for all occurences of a specified filename in a directory and
    returns a list of all present file extensions of this filename.

    In this cas everything after the first dot is considered to be the file
    extension: ``"filename.txt" -> "txt"``, ``"filename.txt.zip" -> "txt.zip"``

    :param targetfilename: a filename without any extensions
    :param directory: only files present in this directory are compared
        to the targetfilename

    :returns: a list of file extensions (str)
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
    """ Returns a list of all substring starting positions in string or an empty
    list if substring is not present in string.

    :param string: a template string
    :param substring: a string, which is looked for in the ``string`` parameter.

    :returns: a list of substring starting positions in the template string
    """
    #TODO: solve with regex? what about '.':
    #return [m.start() for m in re.finditer('(?='+substring+')', string)]
    start = 0
    positions = []
    while True:
        start = string.find(substring, start)
        if start == -1:
            break
        positions.append(start)
        #+1 instead of +len(substring) to also find overlapping matches
        start += 1
    return positions


def toList(variable, types=(basestring, int, float, )):
    """Converts a variable of type string, int, float to a list, containing the
    variable as the only element.

    :param variable: any python object
    :type variable: (str, int, float, others)

    :returns: [variable] or variable
    """
    if isinstance(variable, types):
        return [variable]
    else:
        return variable


def joinpath(path, *paths):
    """Join two or more pathname components, inserting "/" as needed and
    replacing all "\\" by "/".

    :returns: str
    """
    return os.path.join(path, *paths).replace('\\','/')


# --- not yet named section, involving all sorts of calculations --- #
class Memoize(object):
    """A general memoization class, specify a function when creating a new
    instance of the class. The functions return value is returned and stored in
    ``self.memo`` when the instance is called with an argument for the first
    time. Later calls with the same argument return the cached value, instead of
    calling the function again.

    :ivar function: when ``Memoize`` is called this functions return value is
        returned.
    :ivar memo: a dictionary that records the ``function`` return values of
        already called variables.
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
"""Returns the factorial of a number, the results of already calculated numbers
are stored in factorial.memo """


log10factorial = Memoize(lambda n: math.log10(math.factorial(int(n))))
"""Returns the log10 factorial of a number, the results of already calculated
numbers are stored in log10factorial.memo """


# --- some convenience functions to work with arrays and others... --- #
def calcDeviationLimits(value, tolerance, mode):
    """Returns the upper and lower deviation limits for a value and a given
    tolerance, either as relative or a absolute difference.

    :param value: can be a single value or a list of values if a list of values
        is given, the minimal value will be used to calculate the lower limit
        and the maximum value to calculate the upper limit
    :param tolerance: a number used to calculate the limits
    :param mode: either ``absolute`` or ``relative``, specifies how the
        ``tolerance`` should be applied to the ``value``.
    """
    values = toList(value)
    if mode == 'relative':
        lowerLimit = min(values) * (1 - tolerance)
        upperLimit = max(values) * (1 + tolerance)
    elif mode == 'absolute':
        lowerLimit = min(values) - tolerance
        upperLimit = max(values) + tolerance
    else:
        raise Exception('mode %s not specified' %(filepath, ))
    return lowerLimit, upperLimit


def returnArrayFilters(arr1, arr2, limitsArr1, limitsArr2):
    """#TODO: docstring

    :param arr1: #TODO: docstring
    :param arr2: #TODO: docstring
    :param limitsArr1: #TODO: docstring
    :param limitsArr2: #TODO: docstring

    :returns: #TODO: docstring
    """
    posL = bisect.bisect_left(arr1, limitsArr1[0])
    posR = bisect.bisect_right(arr1, limitsArr1[1])
    matchMask = ((arr2[posL:posR] <= limitsArr2[1]) &
                 (arr2[posL:posR] >= limitsArr2[0])
                 )
    return posL, posR, matchMask


def applyArrayFilters(array, posL, posR, matchMask):
    """#TODO: docstring

    :param array: #TODO: docstring
    :param posL: #TODO: docstring
    :param posR: #TODO: docstring
    :param matchMask: #TODO: docstring

    :returns: ``numpy.array``, a subset of the input ``array``.
    """
    return numpy.compress(matchMask, array[posL:posR], axis=0)


# --- data fitting section --- #
class DataFit(object):
    """#TODO: docstring

    :param splines: #TODO: docstring
    :param splineCycles: #TODO: docstring
    :param splineMinKnotPoins: #TODO: docstring
    :param splineOrder: #TODO: docstring
    :param splineInitialKnots: #TODO: docstring
    :param splineSubsetPercentage: #TODO: docstring
    :param splineTerminalExpansion: #TODO: docstring
    :param dependentVar: #TODO: docstring
    :param independentVar: #TODO: docstring
    """
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
        """ #TODO: docstring

        :param dataAveraging: #TODO: docstring
        :param windowSize: #TODO: docstring
        """
        self.dependentVar = numpy.array(self.dependentVarInput,
                                        dtype=numpy.float64
                                        )
        self.independentVar = numpy.array(self.independentVarInput,
                                          dtype=numpy.float64
                                          )

        sortMask = self.independentVar.argsort()
        self.dependentVar = self.dependentVar[sortMask]
        self.independentVar = self.independentVar[sortMask]

        if dataAveraging:
            averagedData = averagingData(self.dependentVar,
                                         windowSize=windowSize,
                                         averagingType=dataAveraging
                                         )
            averagedData = numpy.array(averagedData, dtype=numpy.float64)

            missingNumHigh = numpy.floor((self.independentVar.size
                                          - averagedData.size
                                          ) / 2
                                         )
            missingNumLow = ((self.independentVar.size - averagedData.size)
                             - missingNumHigh
                             )

            self.dependentVar = averagedData
            self.independentVar = self.independentVar[missingNumLow:
                                                      -missingNumHigh]

    def generateSplines(self):
        """#TODO: docstring
        """
        _ = returnSplineList(self.dependentVar, self.independentVar,
                             subsetPercentage=self.splineSubsetPercentage,
                             cycles=self.splineCycles,
                             minKnotPoints=self.splineMinKnotPoins,
                             initialKnots=self.splineInitialKnots,
                             splineOrder=self.splineOrder,
                             terminalExpansion=self.splineTerminalExpansion
                             )
        self.splines = _

    def __getitem__(self, value):
        """#TODO: docstring

        :param value: #TODO: docstring

        :returns: #TODO docstring
        """
        returnValue = numpy.mean([numpy.nan_to_num(currSpline(value)) for
                                  currSpline in self.splines
                                  ])
        return returnValue

    def corrArray(self, inputArray):
        """#TODO: docstring

        :param inputArray: #TODO: docstring

        :returns: #TODO docstring
        """
        outputArray = numpy.vstack([numpy.nan_to_num(currSpline(inputArray))
                                    for currSpline in self.splines
                                    ]).mean(axis=0)
        return outputArray


def averagingData(array, windowSize=None, averagingType='median'):
    """#TODO: docstring

    :param array: #TODO: docstring
    :param windowSize: #TODO: docstring
    :param averagingType: "median" or "mean"

    :returns: #TODO: docstring
    """
    assert averagingType in ['median', 'mean']
    if windowSize is None:
        windowSize = int(len(array) / 50) if int(len(array) / 50) > 100 else 100

    if averagingType == 'median':
        averagedData = runningMedian(array, windowSize)
    elif averagingType == 'mean':
        averagedData = runningMean(array, len(array), windowSize)
    return averagedData


def returnSplineList(dependentVar, independentVar, subsetPercentage=0.4,
                     cycles=10, minKnotPoints=10, initialKnots=200,
                     splineOrder=2, terminalExpansion=0.1
                     ):
    """ #TODO: docstring

    Note: Expects sorted arrays.

    :param dependentVar: #TODO: docstring
    :param independentVar: #TODO: docstring
    :param subsetPercentage: #TODO: docstring
    :param cycles: #TODO: docstring
    :param minKnotPoints: #TODO: docstring
    :param initialKnots: #TODO: docstring
    :param splineOrder: #TODO: docstring
    :param terminalExpansion: expand subsets on both sides

    :returns: #TODO: docstring
    """
    expansions = ddict(list)
    expansionArea = (independentVar[-1] - independentVar[0]) * terminalExpansion
    #adds 100 data points at both ends of the dependent and independent array
    for i in range(100):
        expansions['indUp'].append(independentVar[-1] + expansionArea/100*i)
        expansions['indDown'].append(independentVar[0] -
                                     expansionArea/100*(100-i+1)
                                     )
        expansions['depUp'].append(dependentVar[-1])
        expansions['depDown'].append(dependentVar[0])

    dependentVar = numpy.array(expansions['depDown'] + list(dependentVar) +
                               expansions['depUp'], dtype=numpy.float64
                               )
    independentVar = numpy.array(expansions['indDown'] + list(independentVar) +
                                 expansions['indUp'], dtype=numpy.float64
                                 )

    splineList = list()
    for cycle in range(cycles):
        subset = sorted(random.sample(range(len(dependentVar)),
                                      int(len(dependentVar) * subsetPercentage)
                                      )
                        )
        terminalExpansion

        dependentSubset = dependentVar[subset]
        independentSubset = independentVar[subset]

        minIndVar = independentSubset[minKnotPoints]
        maxIndVar = independentSubset[-minKnotPoints]

        knots = [float(i) * (maxIndVar-minIndVar) / initialKnots + minIndVar
                 for i in range(1, initialKnots)
                 ]
        ## remove knots with less then minKnotPoints data points  ##
        lastKnot = knots[0]
        newKnotList = [lastKnot]
        for knotPos in range(1,len(knots)):
            nextKnot = knots[knotPos]
            numHits = (len(independentSubset[(independentSubset >= lastKnot) &
                       (independentSubset <= nextKnot)])
                       )
            if numHits >= minKnotPoints:
                newKnotList.append(nextKnot)
                lastKnot = nextKnot
        knots = newKnotList

        spline = LSQUnivariateSpline(independentSubset, dependentSubset, knots,
                                     k=splineOrder)
        splineList.append(spline)
    return splineList


def tolerantArrayMatching(referenceArray, matchArray, matchTolerance,
                          matchUnit):
    """#TODO: docstring
    Note: arrays must be sorted

    :param referenceArray: #TODO: docstring
    :param matchArray: #TODO: docstring
    :param matchTolerance: #TODO: docstring
    :param matchUnit: #TODO: docstring

    :returns: #TODO: docstring

    #TODO: change matchUnit to "absolute", "relative" and remove the "*1e-6"
    """
    if matchUnit == 'ppm':
        lowLimMatchArr = matchArray * (1 - matchTolerance*1e-6)
        uppLimMatchArr = matchArray * (1 + matchTolerance*1e-6)
    elif matchUnit == 'da':
        lowLimMatchArr = matchArray - matchTolerance
        uppLimMatchArr = matchArray + matchTolerance
    else:
        raise Exception('wrong matchUnit type specified (da or ppm): ',
                        matchUnit)

    lowerLimitMask = numpy.zeros_like(matchArray, dtype=int)
    upperLimitMask = numpy.zeros_like(matchArray, dtype=int)

    refPosLow = int()
    maxReferenceValue = referenceArray[-1]
    for matchPos, (lowerMatch, upperMatch) in enumerate(zip(lowLimMatchArr,
                                                            uppLimMatchArr
                                                            )
                                                        ):
        if lowerMatch < maxReferenceValue:
            while referenceArray[refPosLow] < lowerMatch:
                refPosLow += 1
            refPosHigh = refPosLow
            #Try except statement because this case can only happen once at the
            #end of the array
            try:
                while referenceArray[refPosHigh] <= upperMatch:
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
