import sys
sys.path.append('D:/Dropbox/python')
import unittest
import pyms.auxiliary as aux


class TestOtherAuxiliaryMethods(unittest.TestCase):
    def test_findAllSubstrings(self):
        self.assertEqual(list(aux.findAllSubstrings('AASDITSLYK', 'S')), [2, 6])
        self.assertEqual(list(aux.findAllSubstrings('AADISPLYSKABSPR', 'SP')), [4, 12])

    def test_toList(self):
        self.assertEqual(aux.toList((1, 2, 3, 'A')), (1, 2, 3, 'A'))
        self.assertEqual(aux.toList('A'), ['A'])
        self.assertEqual(aux.toList(123), [123])

    def test_Factorial(self):
        factorial = aux.Factorial()
        self.assertEqual(factorial[1], 1)
        self.assertEqual(factorial[2], 2)
        self.assertEqual(factorial[3], 6)
        self.assertEqual(factorial[12], 479001600)
        self.assertEqual(factorial.__dict__, {'1': 1, '12': 479001600, '2': 2, '3': 6})


class TestLookUpFileLocationMethods(unittest.TestCase):
    def test_searchFileLocation(self):
        pass
        #aux.searchFileLocation

    def test_matchingFilePaths(self):
        pass
        #aux.matchingFilePaths


#class DataFit
#def averagingData
#def returnSplineList
