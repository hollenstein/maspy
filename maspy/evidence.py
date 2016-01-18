from __future__ import print_function, division
from future.utils import viewkeys, viewvalues, viewitems, listvalues, listitems

from maspy.auxiliary import lazyAttribute
import maspy.proteindb
import maspy.peptidemethods

class PeptideEvidence(maspy.proteindb.PeptideSequence):
    """Summarizes all the evidence (:class:`SpectrumIdentificationItem`) for a certain peptide.
    for other parameters see :class:`PeptideSequence`

    :ivar peptide: amino acid sequence of the peptide including modifications (written in brackets 'AADIT[modification]SLYK')
    :ivar sequence: amino acid sequence of the peptide, corresponds to peptideRef of mzidentml files
    :ivar bestId: containerId of best scoring Sii item
    :ivar siiIds: containerIds of all added Sii items
    :ivar score: best score of all added Sii items
    :ivar scores: scores of all added Sii items
    """
    def __init__(self, peptide, sequence=None):
        sequence = maspy.peptidemethods.removeModifications(peptide) if sequence is None else sequence
        super(PeptideEvidence, self).__init__(sequence)

        self.peptide = peptide
        self.sequence = sequence
        self.bestId = tuple()
        self.score = float()
        self.siiIds = list()
        self.scores = list()

    @lazyAttribute
    def mass(self):
        return maspy.peptidemethods.calcPeptideMass(peptide)

    @classmethod
    def fromPeptideSequence(cls, peptide, peptideSequence):
        newInstance = cls(peptide, peptideSequence.sequence)
        newInstance.isUnique = peptideSequence.isUnique
        newInstance.missedCleavage = peptideSequence.missedCleavage
        newInstance.proteins = peptideSequence.proteins
        newInstance.proteinPositions = peptideSequence.proteinPositions
        return newInstance


class ProteinEvidence(maspy.proteindb.ProteinSequence):
    """ Summarizes all the PeptideEvidence information for a certain protein
    for other paremeters see :class:`ProteinSequence`

    :ivar uniquePsmCount: the sum of PSMs of all unique peptides
    :ivar sharedPsmCount: the sum of PSMs of all shared peptides
    :ivar isValid: should evaluate to True or False, None if unspecified - used to filter data.
    """
    def __init__(self, identifier, sequence, name=None):
        super(ProteinEvidence, self).__init__(identifier, sequence, name)
        self.uniquePsmCount = int()
        self.sharedPsmCount = int()
        self.isValid = None


class EvidenceContainer(object):
    """Container to collect peptide evidence from :class:`SiiContainer` and summarize to protein evidence.

    :ivar db: :class:`ProteinDatabase`, representation of an in silico digest of a fasta file.
    :ivar proteinEvidence: {proteinId: :class:`ProteinEvidence`, ...}
    :ivar peptideEvidence: {peptide: :class:`PeptideEvidence`, ...}

    :ivar uniqueProteins: list of protein ids which have at least one unique peptideEvidence entry.
    :ivar scoreKey: specify what attribute of :class:`SiiItem` should be used as a peptide identification score
    :ivar largerBetter: boolean, True if a larger peptide identification score is better
    """
    def __init__(self, proteinDatabase):
        assert isinstance(proteinDatabase, maspy.proteindb.ProteinDatabase)
        self.db = proteinDatabase
        self.proteinEvidence = dict()
        self.peptideEvidence = dict()

        self.uniqueProteins = list()
        self.scoreKey = None
        self.largerBetter = None

    def calculateCoverage(self):
        """Calcualte the sequence coverage masks for all ProteinEvidence() elements.

        For detailed description see :func:`maspy.proteindb.ProteinDatabase._calculateCoverageMasks`
        """
        maspy.proteindb.ProteinDatabase._calculateCoverageMasks(self.proteinEvidence, self.peptideEvidence)


def generateEvidence(evContainer, siiContainer, scoreKey, largerBetter, peptideKey='peptide'):
    """Summarize all experimental evidence from :class:`SiiContainer` and generate an :class:`EvidenceContainer`.

    :ivar evContainer: :class:`EvidenceContainer` to store the evidence information
    :ivar siiContainer: :class`SiiContainer` instance that contains the experimental evidence
    :ivar peptideKey: 'sequence' -> ignore modification, uses sequence as unique key in
    peptide evidence; 'peptide' -> consider modified peptides as unique entries
    :ivar scoreKey: score attribtue name of :class:`SiiItem`
    :ivar largerBetter: boolean, True if a larger score is better

    NOTE: Peptides which sequence is not present in evContainer.db.peptides are ignored and generate an error message.
    """
    evContainer.scoreKey = scoreKey
    evContainer.largerBetter = largerBetter
    evContainer.peptideKey = peptideKey

    #Summarize psm information into unique peptides
    for sii in siiContainer.getItems(sort=scoreKey, reverse=largerBetter):
        peptide = getattr(sii, peptideKey)
        siiScore = getattr(sii, scoreKey)
        if peptide in evContainer.peptideEvidence:
            evContainer.peptideEvidence[peptide].siiIds.append(sii.containerId)
            evContainer.peptideEvidence[peptide].scores.append(siiScore)
        else:
            try:
                pepEv = PeptideEvidence.fromPeptideSequence(peptide, evContainer.db.peptides[sii.sequence])
                pepEv.bestId = sii.containerId
                pepEv.siiIds.append(sii.containerId)
                pepEv.score = siiScore
                pepEv.scores.append(siiScore)
                evContainer.peptideEvidence[peptide] = pepEv
            except KeyError:
                #TODO: sequences without a protein entry should be accessible via a own list
                #-> one print stateing the number of peptides without a protein in the end
                print('Sequence not found in evContainer.db.peptides: ', sii.sequence, sii.containerId)
                pass

    #Assemble peptide evidence into protein evidence
    for peptide, pepEv in viewitems(evContainer.peptideEvidence):
        for protein in pepEv.proteins:
            if protein in evContainer.proteinEvidence:
                proteinEv = evContainer.proteinEvidence[protein]
            else:
                proteinEv = ProteinEvidence(protein, evContainer.db.proteins[protein].sequence,
                                            evContainer.db.proteins[protein].name
                                            )
                evContainer.proteinEvidence[protein] = proteinEv

            if pepEv.isUnique:
                proteinEv.uniquePeptides.add(peptide)
                proteinEv.uniquePsmCount += len(pepEv.siiIds)
            else:
                proteinEv.sharedPeptides.add(peptide)
                proteinEv.sharedPsmCount += len(pepEv.siiIds)

    #Define proteins which have unique evidence
    evContainer.uniqueProteins = list()
    for proteinEv in viewvalues(evContainer.proteinEvidence):
        if len(proteinEv.uniquePeptides) > 0:
            proteinEv.isUnique = True
            evContainer.uniqueProteins.append(proteinEv.id)
        else:
            proteinEv.isUnique = False
