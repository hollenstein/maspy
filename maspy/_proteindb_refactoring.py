"""
The protein database module provides an interface to protein entries from
fasta files. Using the function importProteinDatabase() imports protein
sequences; which are digested in silico by a specified cleavage rule. Proteins
and peptides can be accessed via the returned ProteinDatabase instance.
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

import io
import re

import numpy
import pyteomics.fasta

import maspy.errors
import maspy.peptidemethods


class ProteinDatabase(object):
    """A database of protein entries from a fasta file.

    In addition to protein entries, it stores peptides generate by an in silico
    digestion of the protein sequences. Also contains the protein and peptide
    relationship.

    :ivar proteins: dict, links proteinIds to ProteinEntry instances
    :ivar peptides: dict, links peptide sequences to PeptideEntry instances
    :ivar ignoreIsoleucine: bool,
    """
    def __init__(self, ignoreIsoleucine=False):
        self.proteins = dict()
        self.peptides = dict()
        self.ignoreIsoleucine = ignoreIsoleucine

        #FUTURE: remove or change... #
        self.info = {'name': '', 'mc': 0, 'cleavageRule': '', 'minLength': 0,
                     'maxLength': 0, 'ignoreIsoleucine': False,
                     'removeNtermM': False}
        self.proteinNames = dict()

        #FUTURE: implement
        self._proteinNameLookUp = dict()

    def getProtein(self, proteinId):
        #Return one protein entry
        raise NotImplementedError

    def _getProteinByName(self, proteinName):
        #After importing a protein database, fill up the self._proteinNameLookup
        #with weak reference. This function can then be used to retrieve protein
        #entries from the weak refs.
        raise NotImplementedError
        #return self._proteinNameLookup[proteinName]()

    def getPeptide(self, peptide):
        #Return one peptide entry
        raise NotImplementedError

    def getProteins(self, proteinId):
        #Iterator for all proteins
        raise NotImplementedError

    def getPeptides(self, peptide):
        #Iterator for all peptides
        raise NotImplementedError

    def getStdSequence(self, sequence):
        """Transform a peptide sequence into a standard sequence, i.e. replace
        all isoleucines with leucines if ignore isoleucine is set.

        :param sequence: an amino acid sequence
        :returns: standard sequence
        """
        if self.ignoreIsoleucine:
            sequence = sequence.replace('I', 'L')
        return sequence

    def _addProtein(self, proteinId, proteinName, sequence, fastaHeader,
                    headerInfo, isDecoy=False, isContaminant=False):
        """#TODO"""
        proteinEntry = ProteinEntry(
            proteinId, proteinName, sequence, fastaHeader, headerInfo,
            isDecoy=isDecoy, isContaminant=isContaminant
        )
        self.proteins[proteinEntry.id] = proteinEntry

    def _addPeptide(self, sequence, proteinId, digestInfo):
        """Add a peptide to the protein database.

        :param sequence: str, amino acid sequence
        :param proteinId: str, proteinId
        :param digestInfo: dict, contains information about the in silico digest
            must contain the keys 'missedCleavage', 'startPos' and 'endPos'
        """
        stdSequence = self.getStdSequence(sequence)

        if stdSequence not in self.peptides:
            self.peptides[stdSequence] = PeptideEntry(
                stdSequence, mc=digestInfo['missedCleavage']
            )
        if sequence not in self.peptides:
            self.peptides[sequence] = self.peptides[stdSequence]

        if proteinId not in self.peptides[stdSequence].proteins:
            #FUTURE: peptide can appear at multiple positions per protein.
            #peptideEntry.addSource(proteinId, startPos, endPos)
            self.peptides[stdSequence].proteins.add(proteinId)
            self.peptides[stdSequence].proteinPositions[proteinId] = (
                                    digestInfo['startPos'], digestInfo['endPos']
                                    )
            self.proteins[proteinId].peptides.add(sequence)

    # --- The methods below are not covered by unit tests --- #

    #FUTURE: remove the direct getitem from proteindb
    def __getitem__(self, key):
        """Uses key to return protein entries.

        :param key: either a protein id or a protein name

        :returns: :class:`Protein`
        """
        if key in self.proteins:
            return self.proteins[key]
        elif key in self.proteinNames:
            return self.proteinNames[key]
        else:
            raise KeyError(key)

    #FUTURE: the generated file is too large and takes too long to import...
    #Maybe only save protein entries, and digestion positions of peptides
    def save(self, path, compress=True):
        raise NotImplementedError
    def load(cls, path, name):
        raise NotImplementedError

    #FUTURE: the way coverage mask calculation is executed may be changed...
    def calculateCoverage(self):
        """Calcualte the sequence coverage masks for all protein entries.

        For a detailed description see :func:`_calculateCoverageMasks()`
        """
        self._calculateCoverageMasks(self.proteins, self.peptides)

    @staticmethod
    def _calculateCoverageMasks(proteindb, peptidedb):
        """Calcualte the sequence coverage masks for all proteindb elements.
        Private method used by :class:`ProteinDatabase`.

        A coverage mask is a numpy boolean array with the length of the protein
        sequence. Each protein position that has been covered in at least one
        peptide is set to True. Coverage masks are calculated for unique and for
        shared peptides. Peptides are matched to proteins according to positions
        derived by the digestion of the FASTA file.

        Alternatively peptides could also be matched to proteins just by
        sequence as it is done in :func:`pyteomics.parser.coverage`, but this is
        not the case here.

        :param proteindb: a dictionary containing :class:`ProteinSequence`
            entries, for example ``ProteinDatabase.proteins``
        :param proteindb: a dictionary containing :class:`PeptideSequence`
            entries, for example ``ProteinDatabase.peptides``

        Sets two attributes for each ``ProteinSequence`` entry:
            ``.coverageMaskUnique`` = coverage mask of unique peptides
            ``.coverageMaskShared`` = coverage mask of shared peptides
        """
        for proteinId, proteinEntry in viewitems(proteindb):
            coverageMaskUnique = numpy.zeros(proteinEntry.length(), dtype='bool')
            for peptide in proteinEntry.uniquePeptides:
                startPos, endPos = peptidedb[peptide].proteinPositions[proteinId]
                coverageMaskUnique[startPos-1:endPos] = True
            coverageMaskShared = numpy.zeros(proteinEntry.length(), dtype='bool')
            for peptide in proteinEntry.sharedPeptides:
                startPos, endPos = peptidedb[peptide].proteinPositions[proteinId]
                coverageMaskShared[startPos-1:endPos] = True
            setattr(proteinEntry, 'coverageMaskUnique', coverageMaskUnique)
            setattr(proteinEntry, 'coverageMaskShared', coverageMaskShared)


class ProteinEntry(object):
    """Representation of a protein.

    :ivar id: str, identifier of the protein, for example a uniprot id.
    :ivar name: str, name of the protein
    :ivar sequence: str, amino acid sequence of the protein
    :ivar fastaHeader: str, the proteins faster header line
    :ivar fastaInfo: dict, the interpreted fasta header as generated when
        using a faster header parsing function, see :func:`fastaParseSgd()`.
    :ivar isDecoy: bool, True if the protein was recognised as a decoy entry.
    :ivar isContaminant: bool, True of the protein was recognised as a
        contaminant entry.
    :ivar peptides: set, peptides derived by the in silico digestion
    :ivar isUnique: bool, True if at least one unique peptide can be assigned to
        the protein
    :ivar uniquePeptides: a set of peptides which can be unambiguously assigned
        to this protein
    :ivar sharedPeptides: a set of peptides which are shared between different
        proteins
    """
    def __init__(self, identifier, name, sequence, fastaHeader, headerInfo,
                 isDecoy=False, isContaminant=False):
        self.id = identifier
        self.name = name
        self.sequence = sequence
        self.isDecoy = isDecoy
        self.isContaminant = isContaminant
        self.fastaHeader = fastaHeader
        self.headerInfo = headerInfo #Formerly fastaInfo

        self.peptides = set()
        self.isUnique = None

        #FUTURE: remove these attributes
        self.uniquePeptides = set()
        self.sharedPeptides = set()

    def mass(self):
        """Returns the number of amino acids of the amino acid sequence."""
        return maspy.peptidemethods.calcPeptideMass(self.sequence)

    def length(self):
        """Returns the mass of the amino acid sequence in dalton."""
        return len(self.sequence)


class PeptideEntry(object):
    """Describes a peptide sequence, derived by the in silico digestion of one
    or several proteins.

    :param sequence: str(), amino acid sequence of the peptide
    :param missedCleavage: int(), number of missed cleavage events, dependens on
        the applied enzyme cleavage rule.
    :param proteins: set(), parent protein IDs which have generated this peptide
        during the in silico digestion.
   :param proteinPositions: dict(), protein IDs point to the start position and
        end position of a peptide in the protein amino acid sequence. One based
        index, i.e. the first amino acid position in a protein is "1".
        ``{proteinId: (startPosition, endPositions) ...}``
    """
    def __init__(self, sequence, mc=None):
        self.sequence = sequence
        self.missedCleavage = mc
        self.isUnique = None
        self.proteins = set()
        self.proteinPositions = dict()

    def length(self):
        """Returns the number of amino acids in the sequence."""
        return len(self.sequence)

    def mass(self):
        """Returns the mass of the amino acid sequence in dalton."""
        return maspy.peptidemethods.calcPeptideMass(self.sequence)


def importProteinDatabase(filePath, proteindb=None, headerParser=None,
        forceId=False, decoyTag='[decoy]', contaminationTag='[cont]',
        ignoreIsoleucine=False,  cleavageRule='[KR]', minLength=5, maxLength=40,
        missedCleavage=0, removeNtermM=False):
    """Generates a :class:`ProteinDatabase` by in silico digestion of proteins
    from a fasta file.

    :param filePath: File path
    :param proteindb: optional an existing :class:`ProteinDatabase` can be
        specified, otherwise a new instance is generated and returned
    :param decoyTag: If a fasta file contains decoy protein entries, they should
        be specified with a sequence tag
    :param contaminationTag: If a fasta file contains contamination protein
        entries, they should be specified with a sequence tag
    :param headerParser: optional, allows specifying an individual headerParser
    :param forceId: bool, if True and no id can be extracted from the fasta
        header the whole header sequence is used as a protein id instead of
        raising an exception.
    :param cleavageRule: cleavage rule expressed in a regular expression, see
        :attr:`maspy.constants.expasy_rules`
    :param missedCleavage: number of allowed missed cleavage sites
    :param removeNtermM: bool, True to consider also peptides with the
        N-terminal Methionine of the protein removed
    :param minLength: int, only yield peptides with length >= minLength
    :param maxLength: int, only yield peptides with length <= maxLength
    :param ignoreIsoleucine: bool, if True treat Isoleucine and Leucine in
        peptide sequences as indistinguishable

    See also :func:`maspy.peptidemethods.digestInSilico`
    """
    if proteindb is None:
        proteindb = ProteinDatabase(ignoreIsoleucine=ignoreIsoleucine)

    # - Add protein entries to the protein database - #
    for fastaHeader, sequence in _readFastaFile(filePath):
        #TODO: function, make protein entry or something like this.
        header, isDecoy = _removeHeaderTag(fastaHeader, decoyTag)
        header, isContaminant = _removeHeaderTag(header, contaminationTag)

        headerInfo = _parseFastaHeader(header, headerParser, forceId)
        proteinId = _idFromHeaderInfo(headerInfo, isDecoy, decoyTag)
        proteinName = _nameFromHeaderInfo(headerInfo, isDecoy, decoyTag)

        if not isDecoy:
            isDecoy = _proteinTagPresent(header, decoyTag)
        if not isContaminant:
            isContaminant = _proteinTagPresent(header, contaminationTag)

        proteindb._addProtein(
            proteinId, proteinName, sequence, fastaHeader, headerInfo,
            isDecoy=isDecoy, isContaminant=isContaminant
        )

    # - Perform an insilico digestion and add peptides to the proteindb - #
    for proteinId in proteindb.proteins:
        sequence = proteindb.proteins[proteinId].sequence
        digestInfo = maspy.peptidemethods.digestInSilico(
            sequence, cleavageRule, missedCleavage, removeNtermM,
            minLength, maxLength
        )
        for sequence, info in digestInfo:
            proteindb._addPeptide(sequence, proteinId, info)

    #Define wheter a peptide is unique to one protein entry
    for peptide, peptideEntry in viewitems(proteindb.peptides):
        if len(peptideEntry.proteins) == 1:
            peptideEntry.isUnique = True
        else:
            peptideEntry.isUnique = False

    #Add peptide as shared or unique to its protein entries
    for peptide, peptideEntry in viewitems(proteindb.peptides):
        for proteinId in peptideEntry.proteins:
            if peptideEntry.isUnique:
                proteindb.proteins[proteinId].uniquePeptides.add(peptide)
            else:
                proteindb.proteins[proteinId].sharedPeptides.add(peptide)

    #Define unique proteins, i.e. have at least one unique peptide
    for proteinEntry in viewvalues(proteindb.proteins):
        if len(proteinEntry.uniquePeptides) > 0:
            proteinEntry.isUnique = True
        else:
            proteinEntry.isUnique = False

    return proteindb


def fastaParseSgd(header):
    """Custom parser for fasta headers in the Saccharomyces Genome Database
    (SGD) format, see www.yeastgenome.org

    :param header: str, protein entry header from a fasta file

    :returns: dict, parsed header
    """
    rePattern = '([\S]+)\s([\S]+).+(\".+\")'
    ID, name, description = re.match(rePattern, header).groups()
    info = {'id': ID, 'name': name, 'description': description}
    return info


def fastaParserSpectraClusterPy(header):
    """Custom parser for fasta headers adapted from
    https://github.com/spectra-cluster/spectra-cluster-py

    :param header: str, protein entry header from a fasta file

    :returns: dict, parsed header
    """

    isUniprot = lambda h: h[0:3] in ['sp|', 'tr|', 'up|']

    if isUniprot(header):
        start = 3
        end = header.find('|', start)
    else:
        start = 0
        breakPositions = [header.find(' '), header.find('|')]
        breakPositions = [i if i > 0 else len(header) for i in breakPositions]
        end = min(breakPositions)

    return {'id': header[start:end]}


def _readFastaFile(filepath):
    """Read a FASTA file and yields tuples of 'header' and 'sequence' entries.

    :param filepath: file path of the FASTA file

    :yields: FASTA entries in the format ('header', 'sequence').
        The 'header' string does not contain the '>' and trailing white spaces.
        The 'sequence' string does not contain trailing white spaces, a '*' at
            the end of the sequence is removed.

    See also :func:`importProteinDatabase` and
    :func:`maspy.peptidemethods.digestInSilico`.
    """
    processSequences = lambda i: ''.join([s.rstrip() for s in i]).rstrip('*')
    processHeaderLine = lambda line: line[1:].rstrip()
    with io.open(filepath, 'r') as openfile:
        #Iterate through lines until the first header is encountered
        try:
            line = openfile.next()
            while line[0] != '>':
                line = openfile.next()
            header = processHeaderLine(line)
            sequences = list()
        except StopIteration:
            errorText = 'File does not contain fasta entries.'
            raise maspy.errors.FileFormatError(errorText)
        #Read remaining lines
        for line in openfile:
            if line[0] == '>':
                yield header, processSequences(sequences)
                header = processHeaderLine(line)
                sequences = list()
            else:
                sequences.append(line)

        #Yield last entry
        if sequences:
            yield header, processSequences(sequences)


def _removeHeaderTag(header, tag):
    """Removes a tag from the beginning of a header string.

    :param header: str
    :param tag: str
    :returns: (str, bool), header without the tag and a bool that indicates
        wheter the tag was present.
    """
    if header.startswith(tag):
        tagPresent = True
        header = header[len(tag):]
    else:
        tagPresent = False
    return header, tagPresent


def _parseFastaHeader(fastaHeader, parser=None, forceId=False):
    """Parses a fasta header and returns extracted information in a dictionary.

    Unless a custom parser is specified, a ``Pyteomics`` function is used, which
    provides parsers for the formats of UniProtKB, UniRef, UniParc and  UniMES
    (UniProt Metagenomic and Environmental Sequences), described at
    `www.uniprot.org <http://www.uniprot.org/help/fasta-headers>_`.

    :param fastaHeader: str, protein entry header from a fasta file
    :param parser: is a function that takes a fastaHeader string and returns a
        dictionary, containing at least the key "id". If None the parser
        function from pyteomics ``pyteomics.fasta.parse()`` is used.
    :param forceId: bool, if True and no id can be extracted from the fasta
        header the whole header sequence is used as a protein id instead of
        raising an exception.

    :returns: dict, describing a fasta header. Minimally contains an 'id' key.
    """
    if parser is None:
        try:
            headerInfo = pyteomics.fasta.parse(fastaHeader)
        except pyteomics.auxiliary.PyteomicsError as raisedPyteomicsError:
            #If forceId is set True, the whole header is used as id
            if forceId:
                headerInfo = {'id': fastaHeader}
            else:
                raise raisedPyteomicsError
    else:
        headerInfo = parser(fastaHeader)
    return headerInfo


def _idFromHeaderInfo(headerInfo, isDecoy, decoyTag):
    """Generates a protein id from headerInfo. If "isDecoy" is True, the
    "decoyTag" is added to beginning of the generated protein id.

    :param headerInfo: dict, must contain a key "id"
    :param isDecoy: bool, determines if the "decoyTag" is added or not.
    :param decoyTag: str, a tag that identifies decoy / reverse protein entries.
    :returns: str, protein id
    """
    proteinId = headerInfo['id']
    if isDecoy:
        proteinId = ''.join((decoyTag, proteinId))

    return proteinId


def _nameFromHeaderInfo(headerInfo, isDecoy, decoyTag):
    """Generates a protein name from headerInfo. If "isDecoy" is True, the
    "decoyTag" is added to beginning of the generated protein name.

    :param headerInfo: dict, must contain a key "name" or "id"
    :param isDecoy: bool, determines if the "decoyTag" is added or not.
    :param decoyTag: str, a tag that identifies decoy / reverse protein entries.
    :returns: str, protein name
    """
    if 'name' in headerInfo:
        proteinName = headerInfo['name']
    else:
        proteinName = headerInfo['id']
    if isDecoy:
        proteinName = ''.join((decoyTag, proteinName))

    return proteinName


def _proteinTagPresent(fastaHeader, tag):
    """Checks wheter a tag string is present in the fastaHeader.

    :param fastaHeader: str, protein entry header from a fasta file
    :returns: bool, True if tag is present in fastaHeader
    """
    return (tag in fastaHeader)
