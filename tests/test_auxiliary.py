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

import os
import sys

import unittest
try:
    import unittest.mock as MOCK
except ImportError:
    import mock as MOCK

sys.path.append(os.path.abspath('..'))
import maspy.auxiliary as MODULE

class TestOtherAuxiliaryMethods(unittest.TestCase):
    def test_findAllSubstrings(self):
        self.assertEqual(MODULE.findAllSubstrings('AASDITSLYK', 'S'), [2, 6])
        self.assertEqual(MODULE.findAllSubstrings('AADISPLYSKABSPR', 'SP'), [4, 12])
        self.assertEqual(MODULE.findAllSubstrings('AAADISPLYSKAABSPR', 'AA'), [0, 1, 11])
        self.assertEqual(MODULE.findAllSubstrings('folder/filname.ext', '/'), [6])
        self.assertEqual(MODULE.findAllSubstrings('folder/filname.ext', '.'), [14])

    def test_toList(self):
        self.assertEqual(MODULE.toList((1, 2, 3, 'A')), (1, 2, 3, 'A'))
        self.assertEqual(MODULE.toList('A'), ['A'])
        self.assertEqual(MODULE.toList(123), [123])

    def test_joinpath(self):
        self.assertEqual(MODULE.joinpath('C:/basedir', 'adir', 'afile.ext'), 'C:/basedir/adir/afile.ext')

    def test_Memoize(self):
        memoize = MODULE.Memoize(lambda n: n)
        memoize(1); memoize(2); memoize(3)
        self.assertDictEqual(memoize.memo, {1: 1, 2: 2, 3: 3})
        self.assertEqual(memoize(1), 1)

        self.assertEqual(MODULE.factorial(10), 3628800)
        self.assertDictEqual(MODULE.factorial.memo, {10: 3628800})
        self.assertEqual(round(MODULE.log10factorial(3), 6), 0.778151)



class TestLookUpFileLocationMethods(unittest.TestCase):
    def mock_os_walk(self):
        filenames = ['test.xml', 'test.txt', 'notest.csv', 'test.txt.zip', 'nofile']
        dirnames = []
        dirpath = 'mockdir'
        return [(dirpath, dirnames, filenames)]
    def mock_os_listdir(self):
        return ['test.xml', 'test.txt', 'notest.csv', 'test.txt.zip', 'nofile']
    def mock_os_path_isfile_true(self):
        return True
    def mock_os_path_isfile_false(self):
        return False

    @MOCK.patch('os.walk', mock_os_walk)
    @MOCK.patch('os.listdir', mock_os_listdir)
    @MOCK.patch('os.path.isfile', mock_os_path_isfile_true)
    def test_searchFileLocation_isfile_true(self):
        self.assertIsNone(MODULE.searchFileLocation('test.tsv', 'csv', 'mockdir', recursive=False))
        self.assertEqual(MODULE.searchFileLocation('test.tsv', 'xml', 'mockdir', recursive=False), 'mockdir/test.xml')
        self.assertEqual(MODULE.searchFileLocation('test.tsv', 'xml', 'mockdir', recursive=True), 'mockdir/test.xml')

    @MOCK.patch('os.listdir', mock_os_listdir)
    @MOCK.patch('os.path.isfile', mock_os_path_isfile_false)
    def test_searchFileLocation_isfile_false(self):
        self.assertIsNone(MODULE.searchFileLocation('test.tsv', 'csv', 'mockdir', recursive=False))
        self.assertIsNone(MODULE.searchFileLocation('test.tsv', 'xml', 'mockdir', recursive=False))


    def test_matchingFilePaths(self):
        pass
        #MODULE.matchingFilePaths

    @MOCK.patch('os.listdir', mock_os_listdir)
    @MOCK.patch('os.path.isfile', mock_os_path_isfile_true)
    def test_listFiletypes_isfile_true(self):
        self.assertEqual(MODULE.listFiletypes('test', 'mockdir'), ['xml', 'txt', 'txt.zip'])

    @MOCK.patch('os.listdir', mock_os_listdir)
    @MOCK.patch('os.path.isfile', mock_os_path_isfile_false)
    def test_listFiletypes_isfile_false(self):
        self.assertEqual(MODULE.listFiletypes('test', 'mockdir'), [])


class TestRunningAverageMethods(unittest.TestCase):
    def test_runningAverageMethods(self):
        testData = [0, 1, 2, 3, 3, 6, 7, 10, 4, 5.5]
        N = len(testData)
        self.assertEqual(MODULE.runningMean(testData, N, 5), [1.8, 3.0, 4.2, 5.8, 6.0, 6.5])
        self.assertEqual(MODULE.runningMedian(testData, 5), [2, 3, 3, 6, 6, 6])


#PartiallySafeReplace
#openSafeReplace
#_isFileAccessible

#MaspyJsonEncoder
#writeJsonZipfile
#writeBinaryItemContainer
#_dumpArrayDictToFile
#_dumpArrayToFile
#_dumpNdarrayToFile
#loadBinaryItemContainer
#_arrayFromBytes

#class DataFit
#def averagingData
#def returnSplineList

#tolerantArrayMatching
