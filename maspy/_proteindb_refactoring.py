"""
Refactoring of the protein data base module (maspy.proteindb)
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

import pyteomics.fasta

import maspy.errors
import maspy.peptidemethods


class ProteinDatabase(object):
    def __init__(self):
        self.proteins = dict()
        self.peptides = dict()

    def addProtein(self, proteinEntry):
        self.proteins[proteinEntry.id] = proteinEntry

    def addPeptides(self, peptide):
        pass
        """
        pepSequence, info
        if ignoreIsoleucine:
            pepSequenceRepr = pepSequence.replace('I', 'L')
        else:
            pepSequenceRepr = pepSequence

        if pepSequenceRepr in proteindb.peptides:
            peptideEntry = proteindb.peptides[pepSequenceRepr]
        else:
            peptideEntry = PeptideEntry(
                pepSequenceRepr, mc=info['missedCleavage']
            )
            proteindb.peptides[pepSequenceRepr] = peptideEntry

        if pepSequence not in proteindb.peptides:
            proteindb.peptides[pepSequence] = peptideEntry

        if proteinId not in peptideEntry.proteins:
            #FUTURE:
            #peptideEntry.addSource(proteinId, startPos, endPos)

            #TODO: peptide can appear at multiple positions per protein.
            peptideEntry.proteins.add(proteinId)
            peptideEntry.proteinPositions[proteinId] = (
                                        info['startPos'], info['endPos']
                                        )
        """


class ProteinEntry(object):
    def __init__(self, identifier, name, sequence, headerInfo, fastaHeader,
                 isDecoy=False, isContaminant=False):
        self.id = identifier
        self.name = name
        self.sequence = sequence
        self.isDecoy = isDecoy
        self.isCont = isContaminant
        self.headerInfo = headerInfo
        self.fastaHeader = fastaHeader

        #FUTURE: remove these attributes
        self.uniquePeptides = set()
        self.sharedPeptides = set()


class PeptideEntry(object):
    def __init__(self, sequence, mc=None):
        self.sequence = sequence
        self.missedCleavage = mc
        self.isUnique = None
        self.proteins = set()
        self.proteinPositions = dict()


def nan():
    """
    header='sp|ID001|geneId1_taxon Protein description 1 OS=organism GN=gene_name1 PE=1 SV=1'
    decoyTag='[rev]'; contTag='[cont]'; headerParser=None; forceId=False; proteindb=None
    filePath='D:/Dropbox/python/msfunctions_0_3/fasta/uniprot_ecoli_K12_refproteome.fasta'
    ignoreIsoleucine = True

    cleavageRule='[KR]'; missedCleavage=0; removeNtermM=True;
    minLength=5; maxLength=55
    """
    proteindb = ProteinDatabase() if proteindb is None else proteindb
    for header, sequence in _readFastaFile(filePath):

        header, isDecoy = _removeHeaderTag(header, decoyTag)
        header, isContaminant = _removeHeaderTag(header, contTag)

        headerInfo = _parseFastaHeader(header, headerParser, forceId)
        proteinId = _idFromHeaderInfo(headerInfo, isDecoy, decoyTag)
        proteinName = _nameFromHeaderInfo(headerInfo, isDecoy, decoyTag)

        proteinEntry = ProteinEntry(
            proteinId, proteinName, sequence, headerInfo, header,
            isDecoy=isDecoy, isContaminant=isContaminant
        )

        proteindb.addProtein(proteinEntry)

        #Perform the insilico digestion
        #Add peptides to the protein database
        digestInfo = maspy.peptidemethods.digestInSilico(
            sequence, cleavageRule, missedCleavage, removeNtermM,
            minLength, maxLength
        )
        #CONTINUE REFACTORING FROM HERE

        #Add peptides to the protein database
        for pepSequence, info in digestInfo:
            if ignoreIsoleucine:
                pepSequenceRepr = pepSequence.replace('I', 'L')
            else:
                pepSequenceRepr = pepSequence

            if pepSequenceRepr in proteindb.peptides:
                peptideEntry = proteindb.peptides[pepSequenceRepr]
            else:
                peptideEntry = PeptideEntry(
                    pepSequenceRepr, mc=info['missedCleavage']
                )
                proteindb.peptides[pepSequenceRepr] = peptideEntry

            if pepSequence not in proteindb.peptides:
                proteindb.peptides[pepSequence] = peptideEntry

            if proteinId not in peptideEntry.proteins:
                #FUTURE:
                #peptideEntry.addSource(proteinId, startPos, endPos)

                #TODO: peptide can appear at multiple positions per protein.
                peptideEntry.proteins.add(proteinId)
                peptideEntry.proteinPositions[proteinId] = (
                                            info['startPos'], info['endPos']
                                            )

    #Add peptide entries to the protein entries, define wheter a peptide can be
    #uniquely assigend to a single protein (.isUnique = True).
    for peptide, peptideEntry in viewitems(proteindb.peptides):
        numProteinMatches = len(peptideEntry.proteins)
        if numProteinMatches == 1:
            peptideEntry.isUnique = True
        elif numProteinMatches > 1:
            peptideEntry.isUnique = False
        else:
            raise Exception('No protein matches in proteindb for peptide' +
                            'sequence: ' + peptide)

        for proteinId in peptideEntry.proteins:
            if peptideEntry.isUnique:
                proteindb.proteins[proteinId].uniquePeptides.add(peptide)
            else:
                proteindb.proteins[proteinId].sharedPeptides.add(peptide)

    #Check protein entries if the digestions generated at least one peptide that
    #is uniquely assigned to the protein (.isUnique = True)
    for proteinEntry in viewvalues(proteindb.proteins):
        if len(proteinEntry.uniquePeptides) > 0:
            proteinEntry.isUnique = True
        else:
            proteinEntry.isUnique = False
    #Note: TODO, altough isoleucin is ignored, the protein entry should only
    #show the actually present ILE / LEU occurence, not any possibilities
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

