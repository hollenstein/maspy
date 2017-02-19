"""
The module "isobar" is still in development. The interface of high and low level
functions is not yet stable!

This module provides methods for working with isobaric tag labeling
strategies, providing access to reporter intensities for quantification.

:Todos:
    Missing reporter ions
        At the moment missing reporter ion mz values are returned as numpy.nan
        but missing intensities as 0, this is because the linear regression
        function for normalizing isotope impurities can't handle nan.

    Simple normalization procedures
        Normalization procedures like equal median or summed intensity for each
        channel could be implemented here. However, these procedures are not
        necessarily unique for Isobaric labeling strategies. Therefore it makes
        more sense to generate a new module containing these functions. Either
        a module called something like "quant" or a separate module like
        "quantnorm".

    Extract reporter ions for one or multiple specfiles
        A function providing similar functionality as used in the example.

    Correct isotope impurities
        Check logic and correctness of iTRAQ8plex and TMT10plex matrix
        transformation. Implement TMT10plex matrix transformation.

    Unit tests for Individual Label strategies (iTRAQ4plex, TMT6plex, ...)

:Example:
    import numpy

    import maspy.reader
    import maspy.isobar

    mzmlPath = 'a_file_path/a_file_name.mzML'
    msrunContainer = maspy.reader.importMzml(mzmlPath)

    reporterMz = [126.127726, 127.124761, 128.134436,
                  129.131471, 130.141145, 131.138180]
    impurityMatrix = [
        [0.00, 0.00, 99.00, 1.00, 0.00],
        [0.00, 0.00, 99.00, 1.00, 0.00],
        [0.00, 0.00, 99.00, 1.00, 0.00],
        [0.00, 0.00, 99.00, 1.00, 0.00],
        [0.00, 0.00, 99.00, 1.00, 0.00],
        [0.00, 0.00, 99.00, 1.00, 0.00]
    ]

    tmt6plex = maspy.isobar.IsobaricTag('tmt6plex', reporterMz,
                                        impurityMatrix, 2, 2)
    mzTolerance = 20e-6
    selector = lambda si: si.msLevel > 1

    reporterArrays = {'mz': [], 'i': [], 'corrI': []}
    for si in msrunContainer.getItems(selector=selector):
        sai = msrunContainer.saic[si.specfile][si.id]
        reporterIons = maspy.isobar._extractReporterIons(
            sai.arrays, tmt6plex.reporterMz, mzTolerance
        )
        reporterIons['corrI'] = tmt6plex.correctIsotopeImpurities(reporterIons['i'])
        for key in reporterIons:
            reporterArrays[key].append(reporterIons[key])
    for key in list(reporterArrays):
        reporterArrays[key] = numpy.array(reporterArrays[key])

    #Plot mz deviation of observed reporter ions
    from matplotlib import pyplot as plt
    fig, ax = plt.subplots(6, figsize=(4, 9), sharex=True, sharey=True)
    bins = numpy.linspace(-10, 10, 40)
    for channel in range(6):
        exactMz = reporterMz[channel]
        mzDev = (1 - reporterArrays['mz'][:, channel] / exactMz) * 1e6
        m = numpy.isfinite(mzDev)
        ax[channel].hist(mzDev[m], bins=bins, color='grey')
        title = ''.join(['Channel ', str(channel+1), ': ', str(int(exactMz))])
        ax[channel].set_title(title)
    ax[channel].set_xlabel('delta m/z [ppm]')
    fig.tight_layout()
    fig.show()
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
        "TMT6plex-HCD" or "iTRAQ4plex".
    :ivar reporterMz: A list of reporter ion mz values observed after
        fragmentation.
    :ivar impurityMatrix: a matrix (2d nested list) that describes reporter ion
        isotope impurities. Each isobaric channel must be present as a row. The
        number of rows must be equal to the number of "reporterMz" values.
    :ivar matrixPreChannels: number of matrix columns with a nominal mass shift
        < 0 (-1, -2,..) in respect to the reporter ion mz value.
    :ivar matrixPostChannels: number of matrix columns with a nominal mass shift
        > 0 (+1, +2,..) in respect to the reporter ion mz value.
    :ivar labelReagents: #TODO: docstring
    """
    def __init__(self, reagentName, reporterMz, impurityMatrix,
                 matrixPreChannels, matrixPostChannels, labelReagents=None):
        self.reagentName = reagentName
        self.reporterMz = reporterMz
        self.impurityMatrix = impurityMatrix
        self.matrixPreChannels = matrixPreChannels
        self.matrixPostChannels = matrixPostChannels
        self._processedMatrix = self._processImpurityMatrix()
        if labelReagents is None:
            self.labelReagents = [
                'reporter'+str(i) for i in range(1, len(reporterMz)+1)
            ]
        else:
            self.labelReagents = labelReagents

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
        raise DeprecationWarning

    def _processImpurityMatrix(self):
        """Process the impurity matrix so that it can be used to correct
        observed reporter intensities.
        """
        processedMatrix = _normalizeImpurityMatrix(self.impurityMatrix)
        processedMatrix = _padImpurityMatrix(
            processedMatrix, self.matrixPreChannels, self.matrixPostChannels
            )
        processedMatrix = _transposeMatrix(processedMatrix)
        return processedMatrix

    def corrImpurities(self, intensities):
        """Corrects observed reporter ion intensities for isotope impurities.

        Because of isotope impurities of the isobaric tagging reagents, the
        individual observed reporter ion intensities have some contributions
        from adjcent reporter ions. To account for these isotope overlaps the
        observed reporter ion intensities can be corrected using an isotope
        impurity matrix, which is supplied by the manufacturer and can vary
        between reagent batches.

        :param intensities: numpy array of observed reporter ion intensities.
        :returns: a numpy array of reporter ion intensities corrected for
            isotope impuritiesself.
        """
        return _correctIsotopeImpurities(self._processedMatrix, intensities)


class Tmt2plex(IsobaricTag):
    """#TODO:

    :ivar reagentName: "TMT2plex", name of the isobaric labeling reagent
    :ivar reporterMz: A list of reporter ion mz values observed after
        fragmentation.
    :ivar impurityMatrix: #TODO
    :ivar labelReagents: A list of label reagent names that correspond to the
        reporter ions of "reporterMz".

    #TODO: Describe init here, use proper matrix, specify hcd or etd
    """
    _reporterMz = {'hcd': [126.127726, 127.131081],
                   'etd': [114.127725, 114.127725]
                   }
    _labelReagents = ['126', '127C']
    def __init__(self, impurityMatrix, fragmentation='hcd'):
        IsobaricTag.__init__(self, 'TMT2plex', self._reporterMz[fragmentation],
            impurityMatrix, 2, 2, labelReagents=self._labelReagents
        )

    def _processImpurityMatrix(self):
        raise NotImplementedError()


class Tmt6plex(IsobaricTag):
    """#TODO:

    :ivar reagentName: "TMT6plex", name of the isobaric labeling reagent
    :ivar reporterMz: A list of reporter ion mz values observed after
        fragmentation.
    :ivar impurityMatrix: #TODO
    :ivar labelReagents: A list of label reagent names that correspond to the
        reporter ions of "reporterMz".

    #TODO: Describe init here, use proper matrix, specify hcd or etd
    """
    _reporterMz = {'hcd': [126.127726, 127.124761, 128.134436,
                           129.131471, 130.141145, 131.138180],
                   'etd': [114.127725, 115.124760, 116.134433,
                           117.131468, 118.141141, 119.138176]
                   }
    _labelReagents = ['126',  '127N', '128C', '129N', '130C', '131']

    def __init__(self, impurityMatrix, fragmentation='hcd'):
        IsobaricTag.__init__(self, 'TMT6plex', self._reporterMz[fragmentation],
            impurityMatrix, 2, 2, labelReagents=self._labelReagents
        )

    def _processImpurityMatrix(self):
        raise NotImplementedError()


class Tmt10plex(IsobaricTag):
    """#TODO:

    :ivar reagentName: "TMT10plex", name of the isobaric labeling reagent
    :ivar reporterMz: A list of reporter ion mz values observed after
        fragmentation.
    :ivar impurityMatrix: #TODO
    :ivar labelReagents: A list of label reagent names that correspond to the
        reporter ions of "reporterMz".

    #TODO: Describe init here, use proper matrix, specify hcd or etd
    """
    #Values from the R package MSnbase
    _reporterMz = {'hcd': [126.127725, 127.124760, 127.131079, 128.128114,
                           128.134433, 129.131468, 129.137787, 130.134822,
                           130.141141, 131.138176],
                   'etd': [114.127725, 115.124760, 114.127725, 115.124760,
                           116.134433, 117.131468, 116.134433, 117.131468,
                           118.141141, 119.138176]
                   }
    _labelReagents = ['126', '127N', '127C', '128N', '128C',
                      '129N', '129C', '130N', '130C', '131']

    def __init__(self, impurityMatrix, fragmentation='hcd'):
        IsobaricTag.__init__(self, 'TMT10plex', self._reporterMz[fragmentation],
            impurityMatrix, 2, 2, labelReagents=self._labelReagents
        )

    def _processImpurityMatrix(self):
        raise NotImplementedError()


class Itraq4plex(IsobaricTag):
    """#TODO:

    :ivar reagentName: "iTRAQ4plex", name of the isobaric labeling reagent
    :ivar reporterMz: A list of reporter ion mz values observed after
        fragmentation.
    :ivar impurityMatrix: #TODO
    :ivar labelReagents: A list of label reagent names that correspond to the
        reporter ions of "reporterMz".

    #TODO: Describe init here
    """
    _labelReagents = ['114', '115', '116', '117']
    _reporterMz = [114.1112, 115.1083, 116.1116, 117.1150]
    _impurityMatrix =[
        [0.00, 1.00, 92.90, 5.90, 0.20],
        [0.00, 2.00, 92.30, 5.60, 0.10],
        [0.00, 3.00, 92.40, 4.50, 0.10],
        [0.10, 4.00, 92.30, 3.50, 0.10]
    ]

    def __init__(self):
        IsobaricTag.__init__(self, 'iTRAQ4plex', self._reporterMz,
            self._impurityMatrix, 2, 2, labelReagents=self._labelReagents)

    def _processImpurityMatrix(self):
        raise NotImplementedError()


class Itraq8plex(IsobaricTag):
    """#TODO:

    :ivar reagentName: "iTRAQ8plex", name of the isobaric labeling reagent
    :ivar reporterMz: A list of reporter ion mz values observed after
        fragmentation.
    :ivar impurityMatrix: #TODO
    :ivar labelReagents: A list of label reagent names that correspond to the
        reporter ions of "reporterMz".

    #TODO: Describe init here
    """
    _labelReagents = ['113', '114', '115', '116', '117', '118', '119', '121']
    _reporterMz = [113.10787, 114.11123, 115.10826, 116.11162,
                   117.11497, 118.11201, 119.11530, 121.12200]
    _impurityMatrix = [
        [0.00, 0.00, 92.87, 6.89, 0.24],
        [0.00, 0.94, 93.00, 5.90, 0.16],
        [0.00, 1.88, 93.12, 4.90, 0.10],
        [0.00, 2.82, 93.21, 3.90, 0.07],
        [0.06, 3.77, 93.29, 2.88, 0.00],
        [0.09, 4.71, 93.29, 1.91, 0.00],
        [0.14, 5.66, 93.33, 0.87, 0.00],
        [0.27, 7.44, 92.11, 0.18, 0.00]
    ]

    def __init__(self):
        IsobaricTag.__init__(self, 'iTRAQ8plex', self._reporterMz,
            self._impurityMatrix, 2, 2, labelReagents=self._labelReagents)

    def _processImpurityMatrix(self):
        raise NotImplementedError()
        processedMatrix = _normalizeImpurityMatrix(self.impurityMatrix)
        _rearrangeItraq8plexMatrix(processedMatrix)
        processedMatrix = _padImpurityMatrix(
            processedMatrix, self.matrixPreChannels, self.matrixPostChannels
            )
        processedMatrix = _transposeMatrix(processedMatrix)
        return processedMatrix


def _extractReporterIons(ionArrays, reporterMz, mzTolerance):
    """Find and a list of reporter ions and return mz and intensity values.

    Expected reporter mz values are searched in "ionArray['mz']" and reported if
    the observed relative deviation is less than specified by "mzTolerance". In
    the case of multiple matches, the one with the minimal deviation is picked.
    If no matching entries are found numpy.nan is returned for the mz value and
    an intensity of 0. The returned arrays are in the order of "reporterMz"
    values.

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
            #reporterIons['i'].append(numpy.nan)
            reporterIons['i'].append(0)
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


def _correctIsotopeImpurities(matrix, intensities):
    """Corrects observed reporter ion intensities for isotope impurities.

    :params matrix: a matrix (2d nested list) containing numbers, each isobaric
        channel must be present as a COLUMN. Use maspy.isobar._transposeMatrix()
        if channels are written in rows.
    :param intensities: numpy array of observed reporter ion intensities.
    :returns: a numpy array of reporter ion intensities corrected for isotope
        impurities.
    """
    correctedIntensities, _ = scipy.optimize.nnls(matrix, intensities)
    return correctedIntensities


def _rearrangeItraq8plexMatrix(matrix):
    """Rearranges an iTRAQ 8plex isotope impurity matrix, so that the impurities
    of 119 and 121 are positioned appropriately.

    #TODO: explain why

    :params matrix: a matrix (2d nested list) containing numbers, each isobaric
        channel must be present as a row.
    """
    matrix[6][3] = matrix[6][4]
    matrix[6][4] = 0
    matrix[7][1] = matrix[7][0]
    matrix[7][0] = 0


def _rearrangeTmt10plexMatrix(matrix):
    """Rearranges a TMT 10plex isotope impurity matrix, so that the impurities
    are positioned appropriately.

    #TODO: explain why

    :params matrix: a matrix (2d nested list) containing numbers, each isobaric
        channel must be present as a row.
    """
    #see https://tools.thermofisher.com/content/sfs/gallery/high/90110-004-TMT10-ratios.jpg
    raise NotImplementedError


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


