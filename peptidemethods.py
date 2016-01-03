from __future__ import print_function

import itertools
import re

import pyteomics.mass

import maspy.constants

def digestInSilico(proteinSequence, cleavageRule='[KR]', missedCleavages=0, removeNtermM=True, minLength=5, maxLength=40):
    """Yields peptides derived from an in silico digest of a polypeptide.

    :param proteinSequence: amino acid sequence of the protein to be digested
    :param cleavageRule: TODO add description
    :param missedCleavages: number of allowed missed cleavage sites
    :param removeNtermM: If True, consider peptides with the n-terminal methionine of the protein removed
    :param minLength: only yield peptides with length >= minLength
    :param maxLength: only yield peptides with length <= maxLength

    Note: An example for specifying N-terminal cleavage at Lysine sites: \\w(?=[K])
    """
    # Return in silico digested peptides, peptide start position, peptide end position
    # Peptide position start at 1 and end at len(proteinSequence)
    passFilter = lambda startPos, endPos: (endPos - startPos >= minLength and endPos - startPos <= maxLength)

    cleavagePosList = set(itertools.chain(map(lambda x: x.end(), re.finditer(cleavageRule, proteinSequence))))
    cleavagePosList.add(len(proteinSequence))
    cleavagePosList = sorted(list(cleavagePosList))
    # Add end of protein as cleavage site if protein doesn't end with specififed cleavage positions
    numCleavageSites = len(cleavagePosList)

    if missedCleavages >= numCleavageSites:
        missedCleavages = numCleavageSites -1

    #Yield protein n-terminal peptides after methionine removal
    if removeNtermM and proteinSequence[0] == 'M':
        for cleavagePos in range(0,missedCleavages+1):
            startPos = 1
            endPos = cleavagePosList[cleavagePos]
            if passFilter(startPos, endPos):
                sequence = proteinSequence[startPos:endPos]
                info = dict()
                info['startPos'] = startPos+1
                info['endPos'] = endPos
                info['missedCleavage'] = cleavagePos
                yield sequence, info

    #Yield protein n-terminal peptides
    if cleavagePosList[0] != 0:
        for cleavagePos in range(0,missedCleavages+1):
            startPos = 0
            endPos = cleavagePosList[cleavagePos]
            if passFilter(startPos, endPos):
                sequence = proteinSequence[startPos:endPos]
                info = dict()
                info['startPos'] = startPos+1
                info['endPos'] = endPos
                info['missedCleavage'] = cleavagePos
                yield sequence, info

    #Yield all remaining peptides, including the c-terminal peptides
    lastCleavagePos = 0
    while lastCleavagePos < numCleavageSites:
        for missedCleavage in range(0, missedCleavages+1):
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
                    yield sequence, info
        lastCleavagePos += 1


# --- Functions to work with peptide sequences --- #
def calcPeptideMass(peptide):
    """Calculate the mass of a peptide. (Should be changed to allow for modifications not present in unimod.org)

    :ivar peptide: peptide sequence, modifications have to be written in the format "[modificationName]"
    and 'modificationName' has to be present in maspy.constants.unimodToMassDict
    """
    unimodMassDict = maspy.constants.unimodToMassDict

    additionalModMass = float()
    unmodPeptide = peptide
    for unimodNumber, unimodMass in unimodMassDict.items():
        try:
            int(unimodNumber)
        except ValueError:
            unimodSymbol = '[' + unimodNumber + ']'
        else:
            unimodSymbol = '[UNIMOD:' + unimodNumber + ']'
        numMod = peptide.count(unimodSymbol)
        unmodPeptide = unmodPeptide.replace(unimodSymbol, '')
        additionalModMass += unimodMass * numMod

    if unmodPeptide.find('[') != -1:
        print(unmodPeptide)
        raise Exception()

    unmodPeptideMass = pyteomics.mass.calculate_mass(unmodPeptide, charge=0)
    modPeptideMass = unmodPeptideMass + additionalModMass
    return modPeptideMass


def removeModifications(peptide):
    """Removes all modifications from a peptide string

    :ivar peptide: peptide sequence, modifications have to be written in the format "[modificationName]"
    :type peptide: str
    """
    while peptide.find('[') != -1:
        peptide = peptide.split('[', 1)[0] + peptide.split(']', 1)[1]
    return peptide


def returnModPositions(peptide, indexStart=1, removeModString='UNIMOD:'):
    """Determines the amino acid positions of all present modifications.

    :ivar peptide: peptide sequence, modifications have to be written in the format "[modificationName]"
    :ivar indexStart: returned amino acids positions of the peptide start with this number (1st amino acid position = indexStart)
    :ivar removeModString: string to remove from the returned modification name

    :return: {modificationName:[position1, position2, ...], ...}
    """
    unidmodPositionDict = dict()
    while peptide.find('[') != -1:
        currModification = peptide.split('[')[1].split(']')[0]
        currPosition = peptide.find('[') - 1
        if currPosition == -1: # move n-terminal modifications to first position
            currPosition = 0
        currPosition += indexStart

        peptide = peptide.replace('['+currModification+']', '', 1)

        if isinstance(removeModString, str):
            currModification = currModification.replace(removeModString, '')
        unidmodPositionDict.setdefault(currModification,list())
        unidmodPositionDict[currModification].append(currPosition)
    return unidmodPositionDict


def calcMhFromMz(mz, charge):
    """Calculate the MH+ value from mz and charge.

    :type mz: float
    :type charge: int
    """
    mh = (mz * charge) - (maspy.constants.atomicMassProton * (charge-1) )
    return mh


def calcMzFromMh(mh,charge):
    """Calculate the mz value from MH+ and charge.

    :type mz: float
    :type charge: int
    """
    mz = ( mh + (maspy.constants.atomicMassProton * (charge-1) ) ) / charge
    return mz


def calcMzFromMass(mass, charge):
    """Calculate the mz value of a peptide from its mass and charge.

    :type mass: float
    :type charge: int
    """
    mz = (mass + (maspy.constants.atomicMassProton * charge) ) / charge
    return mz


def calcMassFromMz(mz, charge):
    """Calculate the mass of a peptide from its mz and charge.

    :type mz: float
    :type charge: int
    """
    mass = (mz - maspy.constants.atomicMassProton) * charge
    return mass
