import sys
import mock
import unittest

import maspy.auxiliary as module

class TestOtherAuxiliaryMethods(unittest.TestCase):
    def test_findAllSubstrings(self):
        self.assertEqual(module.findAllSubstrings('AASDITSLYK', 'S'), [2, 6])
        self.assertEqual(module.findAllSubstrings('AADISPLYSKABSPR', 'SP'), [4, 12])
        self.assertEqual(module.findAllSubstrings('AAADISPLYSKAABSPR', 'AA'), [0, 1, 11])
        self.assertEqual(module.findAllSubstrings('folder/filname.ext', '/'), [6])
        self.assertEqual(module.findAllSubstrings('folder/filname.ext', '.'), [14])

    def test_toList(self):
        self.assertEqual(module.toList((1, 2, 3, 'A')), (1, 2, 3, 'A'))
        self.assertEqual(module.toList('A'), ['A'])
        self.assertEqual(module.toList(123), [123])

    def test_joinpath(self):
        self.assertEqual(module.joinpath('C:/basedir', 'adir', 'afile.ext'), 'C:/basedir/adir/afile.ext')

    def test_Memoize(self):
        memoize = module.Memoize(lambda n: n)
        memoize(1); memoize(2); memoize(3)
        self.assertDictEqual(memoize.memo, {1: 1, 2: 2, 3: 3})
        self.assertEqual(memoize(1), 1)

        self.assertEqual(module.factorial(10), 3628800)
        self.assertDictEqual(module.factorial.memo, {10: 3628800})
        self.assertEqual(round(module.log10factorial(3), 6), 0.778151)



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

    @mock.patch('os.walk', mock_os_walk)
    @mock.patch('os.listdir', mock_os_listdir)
    @mock.patch('os.path.isfile', mock_os_path_isfile_true)
    def test_searchFileLocation_isfile_true(self):
        self.assertIsNone(module.searchFileLocation('test.tsv', 'csv', 'mockdir', recursive=False))
        self.assertEqual(module.searchFileLocation('test.tsv', 'xml', 'mockdir', recursive=False), 'mockdir/test.xml')
        self.assertEqual(module.searchFileLocation('test.tsv', 'xml', 'mockdir', recursive=True), 'mockdir/test.xml')

    @mock.patch('os.listdir', mock_os_listdir)
    @mock.patch('os.path.isfile', mock_os_path_isfile_false)
    def test_searchFileLocation_isfile_false(self):
        self.assertIsNone(module.searchFileLocation('test.tsv', 'csv', 'mockdir', recursive=False))
        self.assertIsNone(module.searchFileLocation('test.tsv', 'xml', 'mockdir', recursive=False))


    def test_matchingFilePaths(self):
        pass
        #module.matchingFilePaths

    @mock.patch('os.listdir', mock_os_listdir)
    @mock.patch('os.path.isfile', mock_os_path_isfile_true)
    def test_listFiletypes_isfile_true(self):
        self.assertEqual(module.listFiletypes('test', 'mockdir'), ['xml', 'txt', 'txt.zip'])

    @mock.patch('os.listdir', mock_os_listdir)
    @mock.patch('os.path.isfile', mock_os_path_isfile_false)
    def test_listFiletypes_isfile_false(self):
        self.assertEqual(module.listFiletypes('test', 'mockdir'), [])


class TestRunningAverageMethods(unittest.TestCase):
    def test_runningAverageMethods(self):
        testData = [0, 1, 2, 3, 3, 6, 7, 10, 4, 5.5]
        N = len(testData)
        self.assertEqual(module.runningMean(testData, N, 5), [1.8, 3.0, 4.2, 5.8, 6.0, 6.5])
        self.assertEqual(module.runningMode(testData, N, 5), [3, 3, 3, 3, 10, 5.5])
        self.assertEqual(module.runningMedian(testData, 5),[2, 3, 3, 6, 6, 6])


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
