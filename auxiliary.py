import os
import copy
from collections import defaultdict as ddict
import functools
import math
import numpy
import operator
import random

from scipy.interpolate import LSQUnivariateSpline

from pyms.mit_stats import runningMode as runningMode
from pyms.mit_stats import runningMean as runningMean
from pyms.mit_stats import runningMedian as runningMedian

# Define constants #
atomicMassH = 1.00782504
atomicMassProton = 1.00727646677

unimodToMassDict = dict()
unimodToMassDict['36'] = 28.031300 # Dimethyl light label
unimodToMassDict['199'] = 32.056407 # Dimethyl medium label
unimodToMassDict['4'] = 57.021464 # Carbamidomethylation
unimodToMassDict['374'] = -1.007825 # Half of a disulfide bridge
unimodToMassDict['7'] = 0.984016 # Deamidated
unimodToMassDict['188'] = 6.020129 # Label:13C(6)
unimodToMassDict['35'] = 15.994915 # Oxidation
unimodToMassDict['21'] = 79.966331 # Phospho
unimodToMassDict['1'] = 42.010565 # Acetyl
unimodToMassDict['27'] = -18.010565 # Glu->pyro-Glu
unimodToMassDict['28'] = -17.026549 # Gln->pyro-Glu
unimodToMassDict['121'] = 114.042927 # GG, ubiquitinlyation residue
unimodToMassDict['DSS'] = 138.068 # Xlink:DSS / BS3
unimodToMassDict['1020'] = 156.078644 # Xlink:DSS, Water-quenched monolink of of DSS/BS3 crosslinker
unimodToMassDict['1356'] = 212.008590 # phosphate-ribosylation: R, (D, E)
unimodToMassDict['213'] = 541.061110 # ADP-Ribosyl, R
unimodToMassDict['5'] = 43.005814 # Carbamyl, pep-n / K / R
unimodToMassDict['3'] = 226.077598 # Biotin K
unimodToMassDict['*'] = 0.0 #Place holder for the second position of a dipeptide modification like a crosslink.

xTandemMassToUniModDict = copy.deepcopy(unimodToMassDict)
xTandemMassToUniModDict[4] = 57.02147
xTandemMassToUniModDict[374] = -1.00783
xTandemMassToUniModDict[1] = 42.01057
xTandemMassToUniModDict = dict([(round(mass, 5), unimod) for unimod, mass in xTandemMassToUniModDict.items()])

"""
This dict contains regular expressions for cleavage rules of the most popular proteolytic enzymes.
The rules were copied from pyteomics http://pythonhosted.org/pyteomics/ and originally taken from
`PeptideCutter tool <http://ca.expasy.org/tools/peptidecutter/peptidecutter_enzymes.html>`_ at Expasy.
"""
expasy_rules = {'arg-c': 'R',
                'asp-n': '\\w(?=D)',
                'bnps-skatole': 'W',
                'caspase 1': '(?<=[FWYL]\\w[HAT])D(?=[^PEDQKR])',
                'caspase 10': '(?<=IEA)D',
                'caspase 2': '(?<=DVA)D(?=[^PEDQKR])',
                'caspase 3': '(?<=DMQ)D(?=[^PEDQKR])',
                'caspase 4': '(?<=LEV)D(?=[^PEDQKR])',
                'caspase 5': '(?<=[LW]EH)D',
                'caspase 6': '(?<=VE[HI])D(?=[^PEDQKR])',
                'caspase 7': '(?<=DEV)D(?=[^PEDQKR])',
                'caspase 8': '(?<=[IL]ET)D(?=[^PEDQKR])',
                'caspase 9': '(?<=LEH)D',
                'chymotrypsin high specificity': '([FY](?=[^P]))|(W(?=[^MP]))',
                'chymotrypsin low specificity': '([FLY](?=[^P]))|(W(?=[^MP]))|(M(?=[^PY]))|(H(?=[^DMPW]))',
                'clostripain': 'R',
                'cnbr': 'M',
                'enterokinase': '(?<=[DE]{3})K',
                'factor xa': '(?<=[AFGILTVM][DE]G)R',
                'formic acid': 'D',
                'glutamyl endopeptidase': 'E',
                'granzyme b': '(?<=IEP)D',
                'hydroxylamine': 'N(?=G)',
                'iodosobenzoic acid': 'W',
                'lysc': 'K',
                'ntcb': '\\w(?=C)',
                'pepsin ph1.3': '((?<=[^HKR][^P])[^R](?=[FLWY][^P]))|((?<=[^HKR][^P])[FLWY](?=\\w[^P]))',
                'pepsin ph2.0': '((?<=[^HKR][^P])[^R](?=[FL][^P]))|((?<=[^HKR][^P])[FL](?=\\w[^P]))',
                'proline endopeptidase': '(?<=[HKR])P(?=[^P])',
                'proteinase k': '[AEFILTVWY]',
                'staphylococcal peptidase i': '(?<=[^E])E',
                'thermolysin': '[^DE](?=[AFILMV])',
                'thrombin': '((?<=G)R(?=G))|((?<=[AFGILTVM][AFGILTVWA]P)R(?=[^DE][^DE]))',
                'trypsin': '([KR](?=[^P]))|((?<=W)K(?=P))|((?<=M)R(?=P))',
                'trypsin simple': '[KR]'
                }

class Factorial():
    def __getitem__(self, n):
        try:
            return getattr(self, str(n))
        except AttributeError:
            setattr(self, str(n), math.factorial(int(n)))
            return getattr(self, str(n))


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
        self.dependentVar = numpy.array(self.dependentVarInput, dtype='float64')
        self.independentVar = numpy.array(self.independentVarInput, dtype='float64')

        sortMask = self.independentVar.argsort()
        self.dependentVar = self.dependentVar[sortMask]
        self.independentVar = self.independentVar[sortMask]

        if dataAveraging:
            if dataAveraging == 'median':
                averagedData = averagingData(self.dependentVar, windowSize=windowSize, averagingType=dataAveraging)
            elif dataAveraging =='mean':
                averagedData = averagingData(self.dependentVar, windowSize=windowSize, averagingType=dataAveraging)
            averagedData = numpy.array(averagedData, dtype='float64')

            missingNumHigh = (self.independentVar.size - averagedData.size) / 2
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
    if averagingType == 'median':
        if windowSize is None:
            windowSize = int(float(len(array)) / (50)) if int(float(len(array) / (50))) > 100 else 100
        averagedData = runningMedian(array, windowSize)
    elif averagingType == 'mean':
        if windowSize is None:
            windowSize = int(float(len(array)) / (50)) if int(float(len(array) / (50))) > 100 else 100
        averagedData = runningMean(array, None, windowSize)
    return averagedData


def returnSplineList(dependentVar, independentVar, subsetPercentage=0.4, cycles=10,
                     minKnotPoints=10, initialKnots=200, splineOrder=2, terminalExpansion=0.1
                     ):
    """Expects sorted arrays.

    :ivar terminalExpansion: expand subsets on both sides
    """
    expansions = ddict(list)
    expansionArea = (independentVar[-1] - independentVar[0]) * terminalExpansion
    for i in range(100):
        expansions['indUp'].append(independentVar[-1] + expansionArea/100*i)
        expansions['indDown'].append(independentVar[0] - expansionArea/100*(100-i+1))
        expansions['depUp'].append(dependentVar[-1])
        expansions['depDown'].append(dependentVar[0])

    dependentVar = numpy.array(expansions['depDown'] + list(dependentVar) + expansions['depUp'], dtype='float64')
    independentVar = numpy.array(expansions['indDown'] + list(independentVar) + expansions['indUp'], dtype='float64')

    splineList = list()
    for cycle in range(cycles):
        subset = sorted(random.sample(xrange(len(dependentVar)), int(len(dependentVar)*subsetPercentage )) )
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
            numHits = len(independentSubset[(independentSubset >= lastKnot) & (independentSubset <= nextKnot)] )
            if numHits >= minKnotPoints:
                newKnotList.append( nextKnot )
                lastKnot = nextKnot
        knots = newKnotList

        spline = LSQUnivariateSpline(independentSubset, dependentSubset, knots, k=splineOrder)
        splineList.append(spline)
    return splineList


def calcMhFromMz(mz, charge):
    """Calculate the MH+ value from mz and charge.

    :type mz: float
    :type charge: int
    """
    mh = (mz * charge) - (atomicMassProton * (charge-1) )
    return mh


def calcMzFromMh(mh,charge):
    """Calculate the mz value from MH+ and charge.

    :type mz: float
    :type charge: int
    """
    mz = ( mh + (atomicMassProton * (charge-1) ) ) / charge
    return mz


def calcMzFromMass(mass, charge):
    """Calculate the mz value of a peptide from its mass and charge.

    :type mass: float
    :type charge: int
    """
    mz = (mass + (atomicMassProton * charge) ) / charge
    return mz


def calcMassFromMz(mz, charge):
    """Calculate the mass of a peptide from its mz and charge.

    :type mz: float
    :type charge: int
    """
    mass = (mz - atomicMassProton) * charge
    return mass


def searchFileLocation(targetFileName, targetFileExtension, rootDirectory):
    """Search for file with specified file extension in all subfolders of specified rootDirectory, returns first matching instance.

    :type targetFileName: str
    :type rootDirectory: str
    :type targetFileExtension: str
    """
    expectedFileName = targetFileName.split('.')[0] + '.' + targetFileExtension
    targetFilePath = None

    for dirpath, dirnames, filenames in os.walk( rootDirectory ):
        for filename in filenames:
            if filename == expectedFileName:
                targetFilePath = os.path.join(dirpath, filename).replace('\\','/')
                break
        if targetFilePath is not None:
            break
    return targetFilePath


def matchingFilePaths(targetFileName, rootDirectory, targetFileExtension=None, selector=None):
    """Search for files in all subfolders of specified rootDirectory, return filepaths of all matching instances.
    :param targetFileName: filename to search for, only the string before the last '.' is used for filename matching
    :param rootDirectory: search directory, including all subdirectories
    :param targetFileExtension: string after the last '.' in the filename, has to be identical if specified

    :param selector: a function which is called with the value of targetFileName
    and has to return True (include value) or False (discard value).
    If no selector is specified, equality to targetFileName is used.
    """
    if targetFileName.find('.') != -1:
        dotPosition = [x for x in findAllSubstrings(targetFileName, '.')][-1]
        targetFileName = targetFileName[:dotPosition]
    matchExtensions = False if targetFileExtension is None else True
    targetFilePaths = list()

    if selector is None:
        selector = functools.partial(operator.eq, targetFileName)

    for dirpath, dirnames, filenames in os.walk( rootDirectory ):
        for filename in filenames:
            if filename.find('.') != -1:
                dotPosition = [x for x in findAllSubstrings(filename, '.')][-1]
                filenameWithoutExtension = filename[:dotPosition]
            else:
                filenameWithoutExtension = filename

            if selector(filenameWithoutExtension):
                if matchExtensions:
                    if not filename.endswith(''.join(('.', targetFileExtension))):
                        continue
                targetFilePaths.append(os.path.join(dirpath, filename).replace('\\','/'))
    return targetFilePaths


def findAllSubstrings(a_str, sub):
    start = 0
    while True:
        start = a_str.find(sub, start)
        if start == -1: return
        yield start
        start +=  1# use start += 1 to find overlapping matches, otherwise len(sub)


def toList(potentialNoList, types=(str, int, float)):
    """ Converts a string / int to a list

    :type potentialString: (str, list)
    """
    if isinstance(potentialNoList, types):
        return [potentialNoList]
    else:
        return potentialNoList

