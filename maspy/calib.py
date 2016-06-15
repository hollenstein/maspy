""" The module ``maspy.calib`` contains methods for data calibration. This
comprises methods to aquire calibration data, to generate calibration functions
and to apply the calibration functions to items in one of the MasPy data
containers.."""

from __future__ import print_function, division, unicode_literals
from future.utils import viewkeys, viewvalues, viewitems, listvalues, listitems

try: # python 2.7
    from itertools import izip as zip
except ImportError: # python 3.x series
    pass
################################################################################
import bisect
from collections import defaultdict as ddict
import numpy


# --- Calibration of ion m/z values of MS1 spectra --- #
def aquireMs1CalibrationData(msrunContainer, specfile, siiArrays=None,
    lockMass=None, **kwargs):
    """Aquire mass error data, observed vs expected m/z, for calibration of MS1
    ion m/z values. Expected m/z values can be of ambient ions with known exact
    masses and from identified precursor masses of peptide spectrum matching
    results.

    :param msrunContainer: intance of :class:`maspy.core.MsrunContainer`
    :param specfile: filename of an ms-run file, used to extract mass error data
    :param siiArrays: optional, a dictionary of numpy.arrays
        Must provide the keys "obsMz" (observed precursor m/z), "excMz" (exact
        calculated m/z of the identified peptide) and "precursorId" (the scan id
        of the MS1 that preceeded the respective MSn scan).
    :param lockMass: None, True or a list of lock mass m/z values
        If ``None`` lock masses are not used to aquire calibration data. If
        ``True`` use the seven predefined lock mass values. Else use the
        specified lock mass m/z values.
    :param scanRange: a list of positive and negative int values
        When using MSn spectrum identifications for aquiring calibration data
        the m/z value of the peptide precursor is used to find the corresponding
        ion m/z value in the preceeding MS1 scan. By adding positive and
        negative integers to the ``scanRange`` parameter the subsequent and
        preceeding MS1 scans are also used to aquire calibration data from the
        same peptide precursor. By default "[-1, 0, 1]".
    :param massTolerance: float or int
        The maximal allowed deviation for matching two masses. By default
        "10 * 1e-6" which corresponds to 10ppm.
    :param toleranceMode: "relative" or "absolute"
        Specifies how the ``massTolerance`` value is applied, by default
        "relative
    :param topIntMatches: bool, by default False
    :param useTopIntIons: bool,  by default False

    :returns: {'siId': numpy.array([str, ...]),
               'rt': numpy.array([float, ...]),
               'obsMz': numpy.array([float, ...]),
               'excMz': numpy.array([float, ...]),
               'relDev': numpy.array([float, ...]),
               'absDev': numpy.array([float, ...]),
               'int': numpy.array([float, ...])
               }
    """
    scanRange = kwargs.get('scanRange', [-1, 0, 1])
    massTolerance = kwargs.get('massTolerance', 10*1e-6)
    toleranceMode = kwargs.get('toleranceMode', 'relative')
    topIntMatches = kwargs.get('topIntMatches', False)
    useTopIntIons = kwargs.get('useTopIntIons', False)

    ms1Arrays = msrunContainer.getArrays(['rt'], specfiles=specfile, sort='rt',
                                         selector=lambda si: si.msLevel==1
                                         )
    ms1Arrays['calibrationMz'] = [list() for x in range(len(ms1Arrays['id']))]

    if lockMass is not None:
        if lockMass == True: # Set default lock mass values
            lockMass = [445.12002, 519.13882, 593.15761, 667.1764, 536.16536,
            		    610.18416, 684.20295]
        lockMassTuples = [(lockMassMz, lockMassMz) for lockMassMz in lockMass]
        for ms1ArrayPos in range(len(ms1Arrays['id'])):
            ms1Arrays['calibrationMz'][ms1ArrayPos].extend(lockMassTuples)

    if siiArrays is not None:
        precursorArrayLookup = {__: _ for _, __ in enumerate(ms1Arrays['id'])}
        for obsMz, excMz, precursorId in zip(siiArrays['obsMz'],
                                             siiArrays['excMz'],
                                             siiArrays['precursorId']
                                             ):
            #find array position of the precursor scan
            precursorArrayPos = precursorArrayLookup[precursorId]
            #also look for the peptide precursor ions in earlier and later MS1
            #scans by modifying the precursorArrayPos according to the numbers
            #in "scanRange"
            for posModifier in scanRange:
                _arrayPos = precursorArrayPos + posModifier
                try:
                    ms1Arrays['calibrationMz'][_arrayPos].append((obsMz, excMz))
                except IndexError:
                    #An IndexError can occur because of the "scanRange" extension
                    #at the end and the beginning of the MS1 scans
                    pass

    calibrationDataMs1 = {_: [] for _ in ['siId', 'rt', 'obsMz', 'excMz',
                                          'relDev', 'absDev', 'int', 'iit']
                          }
    for siId, rtMs1, calibrationMz in zip(ms1Arrays['id'],
                                          ms1Arrays['rt'],
                                          ms1Arrays['calibrationMz']
                                          ):
        if len(calibrationMz) == 0:
            continue

        ionMzListMs1 = msrunContainer.saic[specfile][siId].arrays['mz']
        ionIntListMs1 = msrunContainer.saic[specfile][siId].arrays['i']

        if useTopIntIons:
            ionIntMask = ionIntListMs1.argsort()[::-1][:useTopIntIons]
            ionMzListMs1 = ionMzListMs1[ionIntMask]
            ionIntListMs1 = ionIntListMs1[ionIntMask]

        currCalibrationData = ddict(list)
        for obsMz, excMz in calibrationMz:
            #todo: would be faster by using bisect
            limHigh = obsMz * (1+massTolerance)
            limLow = obsMz * (1-massTolerance)
            pL = bisect.bisect_left(ionMzListMs1, limLow)
            pH = bisect.bisect_right(ionMzListMs1, limHigh)

            if pH - pL <= 0:
                continue
            ionMatchMask = abs(ionMzListMs1[pL:pH] - obsMz).argmin()
            ms1Mz = ionMzListMs1[pL:pH][ionMatchMask]
            ms1Int = ionIntListMs1[pL:pH][ionMatchMask]

            #note: rel dev calculation changed, test dataFit functions!!!
            relDevObs = (1 - ms1Mz / obsMz)
            absDevObs = obsMz - ms1Mz
            relDevExc = (1 - ms1Mz / excMz)
            absDevExc = excMz - ms1Mz
            if abs(relDevObs) <= massTolerance:
                currCalibrationData['siId'].append(siId)
                currCalibrationData['rt'].append(rtMs1)
                currCalibrationData['obsMz'].append(ms1Mz)
                currCalibrationData['excMz'].append(excMz)
                currCalibrationData['relDev'].append(relDevExc)
                currCalibrationData['absDev'].append(absDevExc)
                currCalibrationData['int'].append(ms1Int)
                currCalibrationData['iit'].append(msrunContainer.sic[specfile][siId].iit)
        if len(currCalibrationData['siId']) == 0:
            continue
        for key in currCalibrationData.keys():
            calibrationDataMs1[key].extend(currCalibrationData[key])

    # Convert calibration data into numpy arrays
    for key in calibrationDataMs1.keys():
        calibrationDataMs1[key] = numpy.array(calibrationDataMs1[key])

    return calibrationDataMs1


def timecalMs1DataMedian(msrunContainer, specfile, calibrationData,
                         minDataPoints=50, deviationKey='relDev'):
    """Generates a calibration value for each MS1 scan by calculating the median
    deviation

    :param msrunContainer: intance of :class:`maspy.core.MsrunContainer`
    :param specfile: filename of an ms-run file, used to generate an calibration
        value for each MS1 spectrum item.
    :param calibrationData: a dictionary of ``numpy.arrays`` containing
        calibration data, as returned by :func:`aquireMs1CalibrationData()`
    :param minDataPoints: The minimal number of data points necessary to
        calculate the calibration value, default value is "50". The calibration
        value for each scan is calulated as the median of all calibration data
        points present for this scan. However, if the number of data points is
        less then specified by ``minDataPoints` the data points of the
        preceeding and subsequent scans are added until the minimal number of
        data points is reached.
    :param deviationKey: the ``calibrationData`` key which contains the
        calibration data that should be used.

    :returns: a dictionary containing the calibration values for each MS1
        ``Si``. ``{si.id: {'calibValue': float, 'n': int, 'data': list}``
    """
    corrData = dict()
    _posDict = dict()
    pos = 0
    for si in msrunContainer.getItems(specfiles=specfile, sort='rt',
                                      selector=lambda si: si.msLevel==1
                                      ):
        corrData[si.id] =  {'calibValue': float(), 'n': int(), 'data': list()}
        _posDict[pos] = si.id
        pos += 1

    for siId, deviation in zip(calibrationData['siId'],
                               calibrationData[deviationKey]):
        corrData[siId]['data'].append(deviation)
        corrData[siId]['n'] += 1

    for pos in range(len(corrData)):
        entry = corrData[_posDict[pos]]

        _data = [entry['data']]
        _n = entry['n']
        expansion = 0
        while _n < minDataPoints:
            expansion += 1
            try:
                expData = corrData[_posDict[pos+expansion]]['data']
                _data.append(expData)
                _n += corrData[_posDict[pos+expansion]]['n']
            except KeyError:
                pass
            try:
                expData = corrData[_posDict[pos-expansion]]['data']
                _data.append(expData)
                _n += corrData[_posDict[pos-expansion]]['n']
            except KeyError:
                pass

        if len(entry['data']) > 0:
            median = numpy.median(entry['data'])
            factor = 1
        else:
            median = float()
            factor = 0

        for expData in _data[1:]:
            if len(expData) > 0:
                median += numpy.median(expData) * 0.5
                factor += 0.5
        median = median / factor
        entry['calibValue'] = median
    return corrData


def applyTimeCalMs1(msrunContainer, specfile, correctionData, **kwargs):
    """Applies correction values to the MS1 ion m/z arrays in order to correct
    for a time dependent m/z error.

    :param msrunContainer: intance of :class:`maspy.core.MsrunContainer`,
        containing the :class:`maspy.core.Sai` items of the "specfile".
    :param specfile: filename of an ms-run file to which the m/z calibration
        should be applied
    :param correctionData: a dictionary containing the calibration values for
        each MS1 ``Si``, as returned by :func:`timecalMs1DataMedian()`.
        ``{si.id: {'calibValue': float}``
    :param toleranceMode: "relative" or "absolute"
        Specifies how the ``massTolerance`` value is applied, by default
        "relative".
    """
    toleranceMode = kwargs.get('toleranceMode', 'relative')

    if toleranceMode == 'relative':
        for siId in correctionData:
            calibValue = correctionData[siId]['calibValue']
            msrunContainer.saic[specfile][siId].arrays['mz'] *= (1 + calibValue)
    elif toleranceMode == 'absolute':
        for siId in correctionData:
            calibValue = correctionData[siId]['calibValue']
            msrunContainer.saic[specfile][siId].arrays['mz'] += calibValue
    else:
        raise Exception('#TODO: a proper exception text')


def applyMassCalMs1(msrunContainer, specfile, dataFit, **kwargs):
    """Applies a correction function to the MS1 ion m/z arrays in order to
    correct for a m/z dependent m/z error.

    :param msrunContainer: intance of :class:`maspy.core.MsrunContainer`,
        containing the :class:`maspy.core.Sai` items of the "specfile".
    :param specfile: filename of an ms-run file to which the m/z calibration
        should be applied
    :param dataFit: a :class:`maspy.auxiliary.DataFit` object, containing
        processed calibration data.
    :param toleranceMode: "relative" or "absolute"
        Specifies how the ``massTolerance`` value is applied, by default
        "relative".
    """
    toleranceMode = kwargs.get('toleranceMode', 'relative')

    if toleranceMode == 'relative':
        for si in msrunContainer.getItems(specfile,
                                          selector=lambda si: si.msLevel==1):
            mzArr = msrunContainer.saic[specfile][si.id].arrays['mz']
            corrArr = dataFit.corrArray(mzArr)
            mzArr *= (1 + corrArr)
    elif toleranceMode == 'absolute':
        for si in msrunContainer.getItems(specfile,
                                          selector=lambda si: si.msLevel==1):
            mzArr = msrunContainer.saic[specfile][si.id].arrays['mz']
            corrArr = dataFit.corrArray(mzArr)
            mzArr += corrArr
    else:
        raise Exception('#TODO: a proper exception text')


