"""
The module "isobar" is still in development. The interface of high and low level
functions is not yet stable!

This module provides methods for working with isobaric tag labeling
strategies, providing access to reporter intensities for quantification.

Example:
    isobaricTag = IsobaricTag()
    isobaricTag.setReporterMz(reporterMz)
    isobaricTag.setImpurityMatrix(impurityMatrix)

    mzTolerance = 20e-6

    selector = lambda si: si.msLevel > 1
    for si in msrunContainer.getItems(selector=selector):
        sai = msrunContainer.saic[si.specfile][si.id]
        reporterIons = _extractReporterIons(sai.arrays, isobaricTag.reporterMz,
                                            mzTolerance)
        reporterIons['corrI'] = isobar.correctIsotopeImpurities(reporterIons['i'])
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

import warnings
warnings.warn('Module maspy.isobar is currently in development.',
              ImportWarning)

import bisect
import itertools
import numpy

import scipy.optimize

import maspy.auxiliary as AUX


class IsobaricTag(object):
    """Representation of an isobaric labeling strategy.

    :ivar reagentName: Name to identify isobaric labeling reagent. For example
        "tmt6plex" or "itraq4plex.
    :ivar reporterMz: a list of reporter ion mz values
    :ivar impurityMatrix: a matrix (2d nested list) that describes reporter ion
        isotope impurities. Each isobaric channel must be present as a row. The
        number of rows must be equal to the number of "reporterMz" values.
    """
    def __init__(self, reagentName):
        self.reagentName = reagentName
        self.reporterMz = None
        self.impurityMatrix = None
        self._matrixPreChannels = None
        self._matrixPostChannels = None
        self._processedMatrix = None

    def setReporterMz(self, reporterMz):
        """Define a list of expected reporter ion mz values.

        :param reporterMz: a list of reporter mz values
        """
        self.reporterMz = reporterMz

    def setImpurityMatrix(self, impurityMatrix, preChannels=2, postChannels=2):
        """Add and process an isotope impurity matrix.

        The standard impurity matrix of an isobaric labeling reagent provided by
        the manufacturer contains five columns with a nominal mass shift of
        (-2, -1, 0, +1, +2), hence the default values of 2 for the arguments
        "preChannels" and "postChannels".

        :params impurityMatrix: a matrix (2d nested list) that describes
            reporter ion isotope impurities. Each isobaric channel must be
            present as a row.
        :params preChannels: number of matrix columns with a nominal mass shift
            < 0 (-1, -2,..) in respect to the reporter ion mz value.
        :params postChannels: number of matrix columns with a nominal mass shift
            > 0 (+1, +2,..) in respect to the reporter ion mz value.
        """
        self.impurityMatrix = impurityMatrix
        self._matrixPreChannels = preChannels
        self._matrixPostChannels = postChannels
        self._processImpurityMatrix()

    def _processImpurityMatrix(self):
        """Process the impurity matrix so that it can be used to correct
        observed reporter intensities.
        """
        processedMatrix = _normalizeImpurityMatrix(self.impurityMatrix)
        processedMatrix = _padImpurityMatrix(
            processedMatrix, self._matrixPreChannels, self._matrixPostChannels
            )
        processedMatrix = _transposeMatrix(processedMatrix)
        self._processedMatrix = processedMatrix

    def correctIsotopeImpurities(self, intensities):
        """Corrects observed reporter ion intensities for isotope impurities.

        Because of isotope impurities of the isobaric tagging reagents, the
        individual observed reporter ion intensities have some contributions
        from adjcent reporter ions. To account for these isotope overlaps the
        observed reporter ion intensities can be corrected using an isotope
        impurity matrix, which is supplied by the manufacturer and can vary
        between reagent batches.

        :param intensities: a list or numpy array of observed reporter ion
            intensities.
        :returns: a numpy array of reporter ion intensities corrected for
            isotope impurities.
        """
        return _correctIsotopeImpurities(self._processedMatrix, intensities)


def _extractReporterIons(ionArrays, reporterMz, mzTolerance):
    """Find and a list of reporter ions and return mz and intensity values.

    Expected reporter mz values are searched in "ionArray['mz']" and reported if
    the observed relative deviation is less than specified by "mzTolerance". In
    the case of multiple matches, the one with the minimal deviation is picked.
    If no matching entries are found numpy.nan is returned. The returned arrays
    are in the order of "reporterMz" values.

    :param ionArrays: a dictionary containing two numpy arrays of equal size,
        {"i": an array of ion intensities, "mz" an array of ion mz values}
    :param reporterMz: a list of reporter mz values
    :param mzTolerance: maximum allowed relative mz deviation
    :returns: {'mz': numpy.array(), 'i': numpy.array()}
    """
    reporterIons = {'mz': [], 'i': []}
    for reporterMzValue in reporterMz:
        limHi = reporterMzValue * (1+mzTolerance)
        limLo = reporterMzValue * (1-mzTolerance)
        loPos = bisect.bisect_left(ionArrays['mz'], limLo)
        upPos = bisect.bisect_right(ionArrays['mz'], limHi)

        matchingValues = ionArrays['mz'][loPos:upPos]
        if matchingValues.size == 0:
            reporterIons['i'].append(numpy.nan)
            reporterIons['mz'].append(numpy.nan)
        elif matchingValues.size == 1:
            reporterIons['i'].append(ionArrays['i'][loPos])
            reporterIons['mz'].append(ionArrays['mz'][loPos])
        else:
            mzDeviations = numpy.abs(matchingValues-reporterMzValue)
            minDeviationPos = numpy.argmin(mzDeviations)
            bestMatchArrayPos = range(loPos, upPos)[minDeviationPos]
            reporterIons['i'].append(ionArrays['i'][bestMatchArrayPos])
            reporterIons['mz'].append(ionArrays['mz'][bestMatchArrayPos])

    reporterIons['mz'] = numpy.array(reporterIons['mz'],
                                     dtype=ionArrays['mz'].dtype
                                     )
    reporterIons['i'] = numpy.array(reporterIons['i'],
                                    dtype=ionArrays['i'].dtype
                                    )

    return reporterIons


def _normalizeImpurityMatrix(matrix):
    """Normalize each row of the matrix that the sum of the row equals 1.

    :params matrix: a matrix (2d nested list) containing numbers, each isobaric
        channel must be present as a row.
    :returns: a matrix containing normalized values
    """
    newMatrix = list()
    for line in matrix:
        total = sum(line)
        if total != 0:
            newMatrix.append([i / total for i in line])
        else:
            newMatrix.append(line)
    return newMatrix


def _padImpurityMatrix(matrix, preChannels, postChannels):
    """Align the values of an isotope impurity matrix and fill up with 0.

    NOTE:
        The length of the rows in the "matrix" must be the sum of "preChannels"
        and "postChannels" + 1.

    :params matrix: a matrix (2d nested list) containing numbers, each isobaric
        channel must be present as a row.
    :params preChannels: number of matrix columns with a nominal mass shift < 0
        (-1, -2,..) in respect to the reporter ion mz value.
    :params postChannels: number of matrix columns with a nominal mass shift > 0
        (+1, +2,..) in respect to the reporter ion mz value.

    :returns: extended matrix, where the number of rows is unchanged but the
        length of each row is extend to the number of rows.
    """
    extendedMatrix = list()
    lastMatrixI = len(matrix)-1
    for i, line in enumerate(matrix):
        prePadding = itertools.repeat(0., i)
        postPadding = itertools.repeat(0., lastMatrixI-i)
        newLine = list(itertools.chain(prePadding, line, postPadding))
        extendedMatrix.append(newLine[preChannels:-postChannels])

    return extendedMatrix


def _transposeMatrix(matrix):
    """Convert an impurity matrix to a numpy array and invert columns and rows,
    so that each channel is present as a column.

    :params matrix: a matrix (2d nested list) containing numbers, each isobaric
        channel must be present as a row.
    :returns: transposed matrix
    """
    return numpy.array(matrix).transpose()


def _correctIsotopeImpurities(matrix, intensities):
    """Corrects observed reporter ion intensities for isotope impurities.

    :params matrix: a matrix (2d nested list) containing numbers, each isobaric
        channel must be present as a COLUMN. Use maspy.isobar._transposeMatrix()
        if channels are written in rows.
    :param intensities: a list or numpy array of observed reporter ion
        intensities.
    :returns: a numpy array of reporter ion intensities corrected for isotope
        impurities.
    """
    correctedIntensities, _ = scipy.optimize.nnls(matrix, intensities)
    return correctedIntensities
