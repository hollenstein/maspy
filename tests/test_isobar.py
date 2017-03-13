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

import numpy
import unittest

sys.path.append(os.path.abspath('..'))
import maspy.isobar as MODULE


class TestIsobaricTag(unittest.TestCase):
    def test_InitClass(self):
        reporterMz = [100, 101, 102, 103]
        impurityMatrix = [
            [0.00, 0.99, 0.01],
            [0.00, 0.99, 0.01],
            [0.00, 0.00, 0.00],
            [0.00, 0.90, 0.10]
            ]

        ionMzList = [80, 99.9995, 100.0003, 100.05, 101, 102, 102.1, 103]
        ionIntensityList = [0, 0, 0.99, 0, 1, 0.01, 0, 0.9]
        ionArrays = {'mz': numpy.array(ionMzList),
                     'i': numpy.array(ionIntensityList)
                     }

        isobar = MODULE.IsobaricTag('name', reporterMz, impurityMatrix, 1, 1)

        #Test matrix processing
        processedMatrix = numpy.array([
            [0.99, 0.00, 0.00, 0.00],
            [0.01, 0.99, 0.00, 0.00],
            [0.00, 0.01, 0.00, 0.00],
            [0.00, 0.00, 0.00, 0.90]
            ])
        numpy.testing.assert_almost_equal(isobar._processedMatrix,
                                          processedMatrix)

        #Test intensity correction
        intensities = numpy.array([0.99, 1, 0.01, 0.90])
        expectedIntensities = numpy.array([1., 1., 0., 1.])
        corrIntensities = isobar.corrImpurities(intensities)
        numpy.testing.assert_almost_equal(corrIntensities, expectedIntensities)


class TestIsobarMethods(unittest.TestCase):
    def test_extractReporterIons(self):
        ionMzList = [80, 99.9995, 100.0003, 100.05, 101, 102, 102.1, 103]
        ionIntensityList = [0, 0, 1, 0, 1, 1, 0, 1]
        ionArrays = {'mz': numpy.array(ionMzList),
                     'i': numpy.array(ionIntensityList)
                     }
        reporterMz = [100, 101, 102, 103]
        mzTolerance = 10e-6

        reporterArrays = MODULE._extractReporterIons(ionArrays, reporterMz,
                                                     mzTolerance)
        self.assertEqual(reporterArrays['mz'].size, len(reporterMz))
        self.assertTrue(numpy.all(reporterArrays['i'] == 1))

    def test_correctIsotopeImpurities(self):
        transposedMatrix = numpy.array([
            [0.99, 0.00, 0.00, 0.00],
            [0.01, 0.99, 0.00, 0.00],
            [0.00, 0.00, 0.00, 0.00],
            [0.00, 0.00, 0.00, 0.90]
            ])
        intensities = [0.99, 1, 0.01, 0.90]
        expectedIntensities = numpy.array([1., 1., 0., 1.])
        corrIntensities = MODULE._correctIsotopeImpurities(transposedMatrix,
                                                           intensities)
        self.assertEqual(corrIntensities.size, len(intensities))
        numpy.testing.assert_almost_equal(corrIntensities, expectedIntensities)

    def test_rearrangeItraq8plexMatrix(self):
        #TODO: deprecated
        impurityMatrix = [
            [-2, -1, 0, 1, 2],
            [-2, -1, 0, 1, 2],
            [-2, -1, 0, 1, 2],
            [-2, -1, 0, 1, 2],
            [-2, -1, 0, 1, 2],
            [-2, -1, 0, 1, 2],
            [-2, -1, 0, 1, 2],
            [-2, -1, 0, 1, 2]
        ]

        newMatrix = MODULE._rearrangeItraq8plexMatrix(impurityMatrix)
        self.assertEqual(impurityMatrix[6][3], 2)
        self.assertEqual(impurityMatrix[6][4], 0)
        self.assertEqual(impurityMatrix[7][1], -2)
        self.assertEqual(impurityMatrix[7][0], 0)

    def test_rearrangeTmt10plexMatrix(self):
        #TODO: not implemented yet
        impurityMatrix = [
            [-2, -1, 0, 1, 2],
            [-2, -1, 0, 1, 2],
            [-2, -1, 0, 1, 2],
            [-2, -1, 0, 1, 2],
            [-2, -1, 0, 1, 2],
            [-2, -1, 0, 1, 2],
            [-2, -1, 0, 1, 2],
            [-2, -1, 0, 1, 2],
            [-2, -1, 0, 1, 2],
            [-2, -1, 0, 1, 2]
        ]

        #newMatrix = MODULE._rearrangeTmt10plexMatrix(impurityMatrix)

    def test_normalizeImpurityMatrix(self):
        matrix = [
            [0.00, 92.87, 6.89],
            [0.94, 93.00, 5.90],
            [1.88, 93.12, 4.90],
            [2.82, 93.21, 3.90]
            ]
        normedMatrix = MODULE._normalizeImpurityMatrix(matrix)
        for row in normedMatrix:
            self.assertAlmostEqual(sum(row), 1, places=7)

    def test_padImpurityMatrix(self):
        matrix = [
            [0.00, 92.87, 6.89, 0.24],
            [0.94, 93.00, 5.90, 0.16],
            [1.88, 93.12, 4.90, 0.10],
            [2.82, 93.21, 3.90, 0.07],
            [3.77, 93.29, 2.88, 0.00],
            ]
        extendedMatrix = [
            [92.87, 6.89, 0.24, 0.00, 0.00],
            [0.94, 93.00, 5.90, 0.16, 0.00],
            [0.00, 1.88, 93.12, 4.90, 0.10],
            [0.00, 0.00, 2.82, 93.21, 3.90],
            [0.00, 0.00, 0.00, 3.77, 93.29],
            ]

        resultMatrix = MODULE._padImpurityMatrix(matrix, 1, 2)

        self.assertEqual(len(resultMatrix), len(extendedMatrix))
        self.assertEqual(len(resultMatrix[0]), len(resultMatrix))

        for expectedLine, resultLine in zip(extendedMatrix, resultMatrix):
            self.assertEqual(expectedLine, resultLine)

    def test_transposeMatrix(self):
        matrix = [
            [0.99, 0.01, 0.00, 0.00],
            [0.00, 0.99, 0.00, 0.00],
            [0.00, 0.00, 0.00, 0.00],
            [0.00, 0.00, 0.00, 0.90]
            ]
        transposedMatrix = numpy.array([
            [0.99, 0.00, 0.00, 0.00],
            [0.01, 0.99, 0.00, 0.00],
            [0.00, 0.00, 0.00, 0.00],
            [0.00, 0.00, 0.00, 0.90]
            ])
        resultMatrix = MODULE._transposeMatrix(matrix)
        self.assertTrue(numpy.array_equal(resultMatrix, transposedMatrix))

