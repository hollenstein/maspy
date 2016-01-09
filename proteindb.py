from __future__ import print_function

import itertools
import re

import numpy
import pyteomics.mass

from maspy.auxiliary import lazyAttribute
import maspy.peptidemethods

# --- Protein and peptide related classes --- #
class PeptideSequence(object):
    """Describes a peptide as derived by digestion of one or multiple proteins, can't contain any modified amino acids.
    see also :class:`PeptideEvidence`

    :ivar sequence: amino acid sequence of the peptide
    :ivar missedCleavage: number of missed cleavages, dependens on enzyme specificity
    :ivar proteinList: protein ids that generate this peptide under certain digest condition
    :ivar proteinPositions: {proteinId:(startPosition, endPositions) ...} startposition and endposition of peptide in protein
    eg. 'AADITSLYK' IN 'TAKAADITSLYKEETR':(4,12)
    :ivar mass: peptide mass in Daltons
    :ivar length: number of amino acids
    """
    def __init__(self, sequence, mc=None):
        self.sequence = sequence

        self.missedCleavage = mc
        self.isUnique = None
        self.proteins = set()
        self.proteinPositions = dict()

    @lazyAttribute
    def length(self):
        return len(self.sequence)
    @lazyAttribute
    def mass(self):
        return pyteomics.mass.calculate_mass(self.sequence, charge=0)


class ProteinSequence(object):
    """Describes a protein.

    :ivar id: identifier of the protein eg. UniprotId
    :ivar name: name of the protein
    :ivar sequence: amino acid sequence of the protein
    :ivar isUnique: boolean, True if at least one unique peptide can be assigned to the protein
    :ivar uniquePeptides: a set of peptides which can be unambiguously assigned to this protein
    :ivar sharedPeptides: a set of peptides which are shared between different proteins
    :ivar mass: protein mass in Daltons
    :ivar length: number of amino acids
    :ivar coverageUnique: the number of amino acids in the protein sequence that are coverd by unique peptides
    :ivar coverageShared: the number of amino acids in the protein sequence that are coverd by unique or shared peptides
    """
    def __init__(self, identifier, sequence, name=str()):
        self.id = identifier
        self.name = name
        self.sequence = sequence

        self.isUnique = None
        self.uniquePeptides = set()
        self.sharedPeptides = set()

    @lazyAttribute
    def mass(self):
        return pyteomics.mass.calculate_mass(self.sequence, charge=0)

    @lazyAttribute
    def length(self):
        return len(self.sequence)


class ProteinDatabase(object):
    """Describes proteins and peptides generated by an in silico digestion of proteins.

    :ivar peptides: {sequence:PeptideSequence(), ...} contains elements of :class:`PeptideSequence` derived by an in silico digest of the proteins
    :ivar proteins: {proteinId:Protein(), proteinId:Protein()}, used to access :class:`ProteinSequence` elements by their id
    :ivar proteinNames: {proteinName:Protein(), proteinName:Protein()}, alternative way to access :class:`ProteinSequence` elements by their names
    """
    # A container for Protein or ProteinEvidence objects
    def __init__(self):
        self.peptides = dict()
        self.proteins = dict()
        self.proteinNames = dict()

    def __getitem__(self, key):
        """Uses key to return protein entries :class:`Protein`.

        :ivar key: either a protein id or a protein name
        """
        if key in self.proteins:
            return self.proteins[key]
        elif key in self.proteinNames:
            return self.proteinNames[key]
        else:
            raise KeyError(key)

    def calculateCoverage(self):
        """Calcualte the sequence coverage masks for all ProteinEvidence() elements.

        For detailed description see :func:`_calculateCoverageMasks`
        """
        self._calculateCoverageMasks(self.proteins, self.peptides)

    @staticmethod
    def _calculateCoverageMasks(proteindb, peptidedb):
        """Calcualte the sequence coverage masks for all proteindb elements.
        Private method used by :class:`ProteinDatabase` and :class:`EvidenceContainer`

        A coverage mask is a numpy boolean array with the length of the protein sequence.
        Each protein position that has been covered in at least one peptide is set to True.
        Coverage masks are calculated for unique and for shared peptides. Peptides are
        matched to proteins according to positions derived by the digestion of the FASTA file.

        Alternatively peptides could also be matched to proteins just by sequence as it is
        done in :func:`pyteomics.parser.coverage`, but this is not the case here.

        :ivar :attr:`PeptideSequence.coverageMaskUnique`: coverage mask of unique peptides
        :ivar :attr:`PeptideSequence.coverageMaskShared`: coverage mask of shared peptides
        """
        for proteinId, proteinEntry in proteindb.items():
            coverageMaskUnique = numpy.zeros(proteinEntry.length, dtype='bool')
            for peptide in proteinEntry.uniquePeptides:
                startPos, endPos = peptidedb[peptide].proteinPositions[proteinId]
                coverageMaskUnique[startPos-1:endPos] = True
            coverageMaskShared = numpy.zeros(proteinEntry.length, dtype='bool')
            for peptide in proteinEntry.sharedPeptides:
                startPos, endPos = peptidedb[peptide].proteinPositions[proteinId]
                coverageMaskShared[startPos-1:endPos] = True
            setattr(proteinEntry, 'coverageMaskUnique', coverageMaskUnique)
            setattr(proteinEntry, 'coverageMaskShared', coverageMaskShared)


# --- import for ProteinDatabase class --- #
def importProteinDatabase(filePath, proteindb=None, decoyTag='[decoy]', contaminationTag='[cont]', headerParser=None, forceId=False,
                          cleavageRule='[KR]', minLength=5, maxLength=40, missedCleavage=2, ignoreIsoleucine=False, removeNtermM=True,
                          ):
    """Generates a :class:`ProteinContainer` and :class:`PeptideContainer` by in silico digestion of proteins from a fasta file.

    :param filePath: File path
    :param proteindb: optional an existing :class:`ProteinDatabase` can be specified, otherwise a new instance is generated and returned
    :param ignoreIsoleucine: If True, treat I and L in peptide sequence as indistinguishable
    :param decoyTag: If a fasta file contains decoy protein entries, they should be specified with a sequence tag
    :param contaminationTag: If a fasta file contains contamination protein entries, they should be specified with a sequence tag
    :param headerParser: optional a headerParser can be specified TODO: describe how a parser looks like
    :param forceId: If set True if no id can be extracted from the fasta header, the whole header sequence is used as a protein id
    :param cleavageRule: TODO add description
    :param minLength: only yield peptides with length >= minLength
    :param maxLength: only yield peptides with length <= maxLength
    :param ignoreIsoleucine: If True treat I and L in peptide sequence as indistinguishable
    :param missedCleavages: number of allowed missed digestion sites
    :param removeNtermM: If True, consider peptides with the n-terminal methionine of the protein removed

    See also :func:`maspy.peptidemethods.digestInSilico`
    """
    proteindb = ProteinDatabase() if proteindb is None else proteindb
    fastaRead = _readFastaFile(filePath)

    for header, sequence, comments in fastaRead:
        proteinTags = list()
        if header.startswith(decoyTag):
            isDecoy = True
            header = header.replace(decoyTag, '')
            proteinTags.append(decoyTag)
        else:
            isDecoy = False
        if header.startswith(contaminationTag):
            isCont = True
            header = header.replace(contaminationTag, '')
            proteinTags.append(contaminationTag)
        else:
            isCont = False

        headerInfo = _extractFastaHeader(header, headerParser, forceId)
        proteinId = ''.join(itertools.chain(proteinTags, [headerInfo['id']]))
        if 'name' in headerInfo:
            proteinName = ''.join(itertools.chain(proteinTags, [headerInfo['name']]))
        else:
            proteinName = proteinId

        if proteinId not in proteindb.proteins:
            protein = ProteinSequence(proteinId, sequence)
            protein.name = proteinName
            protein.fastaHeader = header
            protein.fastaInfo = headerInfo
            proteindb.proteins[protein.id] = protein

        for unmodPeptide, info in maspy.peptidemethods.digestInSilico(sequence, cleavageRule, missedCleavage,
                                                                      removeNtermM, minLength, maxLength
                                                                      ):
            if ignoreIsoleucine:
                unmodPeptideNoIsoleucine = unmodPeptide.replace('I', 'L')
                if unmodPeptideNoIsoleucine in proteindb.peptides:
                    currPeptide = proteindb.peptides[unmodPeptideNoIsoleucine]
                else:
                    currPeptide = PeptideSequence(unmodPeptideNoIsoleucine, mc=info['missedCleavage'])
                    proteindb.peptides[unmodPeptideNoIsoleucine] = currPeptide

                if unmodPeptide not in proteindb.peptides:
                    proteindb.peptides[unmodPeptide] = currPeptide
            else:
                if unmodPeptide in proteindb.peptides:
                    currPeptide = proteindb.peptides[unmodPeptide]
                else:
                    currPeptide = PeptideSequence(unmodPeptide, mc=info['missedCleavage'])
                    proteindb.peptides[unmodPeptide] = currPeptide

            if proteinId not in currPeptide.proteins:
                currPeptide.proteins.add(proteinId)
                currPeptide.proteinPositions[proteinId] = (info['startPos'], info['endPos'])

    for peptide, peptideEntry in proteindb.peptides.items():
        numProteinMatches = len(peptideEntry.proteins)
        if numProteinMatches == 1:
            peptideEntry.isUnique = True
        elif numProteinMatches > 1:
            peptideEntry.isUnique = False
        else:
            print('No protein matches in proteindb for peptide sequence: ', peptide)

        for proteinId in peptideEntry.proteins:
            if peptideEntry.isUnique:
                proteindb.proteins[proteinId].uniquePeptides.add(peptide)
            else:
                proteindb.proteins[proteinId].sharedPeptides.add(peptide)

    for proteinEntry in proteindb.proteins.values():
        if len(proteinEntry.uniquePeptides) > 0:
            proteinEntry.isUnique = True
        else:
            proteinEntry.isUnique = False
    #Note: TODO, altough isoleucin is ignored, the protein entry should only show the actually present ILE / LEU occurence
    return proteindb


def _readFastaFile(fastaFileLocation):
    """Reads a fasta file and seperates entries into 'header' and 'sequence'.)

    :param fastaFileLocation: File path of the fasta file
    :returns : A list of protein entry touples: [(fasta header, sequence, comments), ...]
    Comments are optional entries of fasta files between the fasta header and the sequence,
    starting with either ";" or ">".

    See also :func:`returnDigestedFasta` and :func:`maspy.peptidemethods.digestInSilico`
    """
    fastaPattern = '(?P<header>([>;].+\n)+)(?P<sequence>[A-Z\*\n]+)' #[\*\n)
    with open(fastaFileLocation, 'r') as openfile:
        readfile = openfile.read()
        entries = list()

        rePattern = re.compile(fastaPattern, re.VERBOSE)
        reResult = rePattern.finditer(str(readfile))

        comments = str()
        isHeader = True
        for entry in reResult:
            if isHeader:
                header = entry.group('header').replace('>', '').strip()
                if entry.group('sequence').strip():
                    sequence = entry.group('sequence').replace('\n', '').replace('\r', '').strip('*')
                    entries.append([header, sequence, comments])
                else:
                    isHeader = False
            else:
                comments += entry.group('header')
                if entry.group('sequence').strip():
                    sequence = entry.group('sequence').replace('\n', '').replace('\r', '').strip('*')
                    entries.append((header, sequence, comments))
                    isHeader = True
                    comments = str()
    return entries


def _extractFastaHeader(fastaHeader, parser=None, forceId=False):
    """
    :param parser: is a function that takes a fastaHeader string and returns a dictionary, containing at least the key 'id'
    if None is specified uses the parser function from pyteomics :func:`pyteomics.fasta.parse()`
    """
    if parser is None:
        try:
            headerInfo = pyteomics.fasta.parse(fastaHeader)
        except pyteomics.auxiliary.PyteomicsError as pyteomicsError:
            #If exceptError is set True, it forces to return the whole header as id
            if forceId:
                headerInfo = {'id':fastaHeader}
            else:
                raise pyteomicsError
    else:
        headerInfo = parser(fastaHeader)
    return headerInfo


def fastaParseSgd(header):
    ID, name, description = re.match('([\S]+)\s([\S]+).+(\".+\")', header).groups()
    info = {'id':ID, 'name':name, 'description':description}
    return info
