import sys
sys.path.append('D:/Dropbox/python')
import unittest
import maspy.auxiliary as module

class TestOtherAuxiliaryMethods(unittest.TestCase):
    def test_findAllSubstrings(self):
        self.assertEqual(list(module.findAllSubstrings('AASDITSLYK', 'S')), [2, 6])
        self.assertEqual(list(module.findAllSubstrings('AADISPLYSKABSPR', 'SP')), [4, 12])

    def test_toList(self):
        self.assertEqual(module.toList((1, 2, 3, 'A')), (1, 2, 3, 'A'))
        self.assertEqual(module.toList('A'), ['A'])
        self.assertEqual(module.toList(123), [123])

    def test_joinpath(self):
        self.assertEqual(module.joinpath('C:/basedir', 'adir', 'afile.ext'), 'C:/basedir/adir/afile.ext')

    def test_Factorial(self):
        factorial = module.Factorial()
        self.assertEqual(factorial[1], 1)
        self.assertEqual(factorial[2], 2)
        self.assertEqual(factorial[3], 6)
        self.assertEqual(factorial[12], 479001600)
        self.assertEqual(factorial.__dict__, {'1': 1, '12': 479001600, '2': 2, '3': 6})


class TestLookUpFileLocationMethods(unittest.TestCase):
    def test_searchFileLocation(self):
        pass
        #module.searchFileLocation

    def test_matchingFilePaths(self):
        pass
        #module.matchingFilePaths


#class DataFit
#def averagingData
#def returnSplineList
