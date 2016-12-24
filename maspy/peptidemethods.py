"""
provides functions to work with peptide
  sequences, mass to charge ratios and modifications and calvulation
  of masses.
"""
######################## Python 2 and 3 compatibility #########################
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
import itertools
import re

import pyteomics.mass

import maspy.constants


def digestInSilico(proteinSequence, cleavageRule='[KR]', missedCleavage=0,
                   removeNtermM=True, minLength=5, maxLength=55):
    """Returns a list of peptide sequences and cleavage information derived
    from an in silico digestion of a polypeptide.

    :param proteinSequence: amino acid sequence of the poly peptide to be
        digested
    :param cleavageRule: cleavage rule expressed in a regular expression, see
        :attr:`maspy.constants.expasy_rules`
    :param missedCleavage: number of allowed missed cleavage sites
    :param removeNtermM: booo, True to consider also peptides with the
        N-terminal methionine of the protein removed
    :param minLength: int, only yield peptides with length >= minLength
    :param maxLength: int, only yield peptides with length <= maxLength

    :returns: a list of resulting peptide enries. Protein positions start with
        ``1`` and end with ``len(proteinSequence``. ::

            [(peptide amino acid sequence,
              {'startPos': int, 'endPos': int, 'missedCleavage': int}
              ), ...
             ]

    .. note::
        This is a regex example for specifying N-terminal cleavage at lysine
        sites ``\\w(?=[K])``
    """
    passFilter = lambda startPos, endPos: (endPos - startPos >= minLength and
                                           endPos - startPos <= maxLength
                                           )
    _regexCleave = re.finditer(cleavageRule, proteinSequence)

    cleavagePosList = set(itertools.chain(map(lambda x: x.end(), _regexCleave)))
    cleavagePosList.add(len(proteinSequence))
    cleavagePosList = sorted(list(cleavagePosList))
    #Add end of protein as cleavage site if protein doesn't end with specififed
    #cleavage positions
    numCleavageSites = len(cleavagePosList)

    if missedCleavage >= numCleavageSites:
        missedCleavage = numCleavageSites -1

    digestionresults = list()
    #Generate protein n-terminal peptides after methionine removal
    if removeNtermM and proteinSequence[0] == 'M':
        for cleavagePos in range(0, missedCleavage+1):
            startPos = 1
            endPos = cleavagePosList[cleavagePos]
            if passFilter(startPos, endPos):
                sequence = proteinSequence[startPos:endPos]
                info = dict()
                info['startPos'] = startPos+1
                info['endPos'] = endPos
                info['missedCleavage'] = cleavagePos
                digestionresults.append((sequence, info))

    #Generate protein n-terminal peptides
    if cleavagePosList[0] != 0:
        for cleavagePos in range(0, missedCleavage+1):
            startPos = 0
            endPos = cleavagePosList[cleavagePos]
            if passFilter(startPos, endPos):
                sequence = proteinSequence[startPos:endPos]
                info = dict()
                info['startPos'] = startPos+1
                info['endPos'] = endPos
                info['missedCleavage'] = cleavagePos
                digestionresults.append((sequence, info))

    #Generate all remaining peptides, including the c-terminal peptides
    lastCleavagePos = 0
    while lastCleavagePos < numCleavageSites:
        for missedCleavage in range(0, missedCleavage+1):
            nextCleavagePos = lastCleavagePos + missedCleavage + 1
            if nextCleavagePos < numCleavageSites:
                startPos = cleavagePosList[lastCleavagePos]
                endPos = cleavagePosList[nextCleavagePos]
                if passFilter(startPos, endPos):
                    sequence = proteinSequence[startPos:endPos]
                    info = dict()
                    info['startPos'] = startPos+1
                    info['endPos'] = endPos
                    info['missedCleavage'] = missedCleavage
                    digestionresults.append((sequence, info))
        lastCleavagePos += 1

    return digestionresults


# --- Functions to work with peptide sequences --- #
def calcPeptideMass(peptide, **kwargs):
    """Calculate the mass of a peptide.

    :param aaMass: A dictionary with the monoisotopic masses of amino acid
        residues, by default :attr:`maspy.constants.aaMass`
    :param aaModMass: A dictionary with the monoisotopic mass changes of
        modications, by default :attr:`maspy.constants.aaModMass`
    :param elementMass: A dictionary with the masses of chemical elements, by
        default ``pyteomics.mass.nist_mass``
    :param peptide: peptide sequence, modifications have to be written in the
        format "[modificationId]" and "modificationId" has to be present in
        :attr:`maspy.constants.aaModMass`

    #TODO: change to a more efficient way of calculating the modified mass, by
    first extracting all present modifications and then looking up their masses.
    """
    aaMass = kwargs.get('aaMass', maspy.constants.aaMass)
    aaModMass = kwargs.get('aaModMass', maspy.constants.aaModMass)
    elementMass = kwargs.get('elementMass', pyteomics.mass.nist_mass)

    addModMass = float()
    unmodPeptide = peptide
    for modId, modMass in viewitems(aaModMass):
        modSymbol = '[' + modId + ']'
        numMod = peptide.count(modSymbol)
        if numMod > 0:
            unmodPeptide = unmodPeptide.replace(modSymbol, '')
            addModMass += modMass * numMod

    if unmodPeptide.find('[') != -1:
        print(unmodPeptide)
        raise Exception('The peptide contains modification, ' +
                        'not present in maspy.constants.aaModMass'
                        )

    unmodPeptideMass = sum(aaMass[i] for i in unmodPeptide)
    unmodPeptideMass += elementMass['H'][0][0]*2 + elementMass['O'][0][0]
    modPeptideMass = unmodPeptideMass + addModMass
    return modPeptideMass


def removeModifications(peptide):
    """Removes all modifications from a peptide string and return the plain
    amino acid sequence.

    :param peptide: peptide sequence, modifications have to be written in the
        format "[modificationName]"
    :param peptide: str

    :returns: amino acid sequence of ``peptide`` without any modifications
    """
    while peptide.find('[') != -1:
        peptide = peptide.split('[', 1)[0] + peptide.split(']', 1)[1]
    return peptide


def returnModPositions(peptide, indexStart=1, removeModString='UNIMOD:'):
    """Determines the amino acid positions of all present modifications.

    :param peptide: peptide sequence, modifications have to be written in the
        format "[modificationName]"
    :param indexStart: returned amino acids positions of the peptide start with
        this number (first amino acid position = indexStart)
    :param removeModString: string to remove from the returned modification name

    :return: {modificationName:[position1, position2, ...], ...}

    #TODO: adapt removeModString to the new unimod ids in
    #maspy.constants.aaModComp ("UNIMOD:X" -> "u:X") -> also change unit tests.
    """
    unidmodPositionDict = dict()
    while peptide.find('[') != -1:
        currModification = peptide.split('[')[1].split(']')[0]
        currPosition = peptide.find('[') - 1
        if currPosition == -1: # move n-terminal modifications to first position
            currPosition = 0
        currPosition += indexStart

        peptide = peptide.replace('['+currModification+']', '', 1)

        if removeModString:
            currModification = currModification.replace(removeModString, '')
        unidmodPositionDict.setdefault(currModification,list())
        unidmodPositionDict[currModification].append(currPosition)
    return unidmodPositionDict


# --- Functions to transform mass to mz values --- #
def calcMhFromMz(mz, charge):
    """Calculate the MH+ value from mz and charge.

    :param mz: float, mass to charge ratio (Dalton / charge)
    :param charge: int, charge state

    :returns: mass to charge ratio of the mono protonated ion (charge = 1)
    """
    mh = (mz * charge) - (maspy.constants.atomicMassProton * (charge-1) )
    return mh


def calcMzFromMh(mh, charge):
    """Calculate the mz value from MH+ and charge.

    :param mh: float, mass to charge ratio (Dalton / charge) of the mono
        protonated ion
    :param charge: int, charge state

    :returns: mass to charge ratio of the specified charge state
    """
    mz = (mh + (maspy.constants.atomicMassProton * (charge-1))) / charge
    return mz


def calcMzFromMass(mass, charge):
    """Calculate the mz value of a peptide from its mass and charge.

    :param mass: float, exact non protonated mass
    :param charge: int, charge state

    :returns: mass to charge ratio of the specified charge state
    """
    mz = (mass + (maspy.constants.atomicMassProton * charge)) / charge
    return mz


def calcMassFromMz(mz, charge):
    """Calculate the mass of a peptide from its mz and charge.

    :param mz: float, mass to charge ratio (Dalton / charge)
    :param charge: int, charge state

    :returns: non protonated mass (charge = 0)
    """
    mass = (mz - maspy.constants.atomicMassProton) * charge
    return mass
