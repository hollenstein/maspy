"""
The protein database module allows the import of protein sequences from fasta
files, parsing of fasta entry headers and performing in silico digestion by
specified cleavage rules to generate peptides.
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
import itertools
import json
import re
import zipfile

import numpy
import pyteomics
import pyteomics.fasta

import maspy.auxiliary as aux
import maspy.peptidemethods


# --- Protein and peptide related classes --- #
class PeptideSequence(object):
    """Describes a peptide as derived by digestion of one or multiple proteins,
    can't contain any modified amino acids.

    :param sequence: amino acid sequence of the peptide
    :param missedCleavage: number of missed cleavages, dependens on enzyme
        specificity
    :param proteins: protein ids that generate this peptide under certain
        digest condition
    :param proteinPositions: start position and end position of a peptide in a
        protein sequence. One based index, ie the first protein position is "1".
        ``{proteinId:(startPosition, endPositions) ...}``
    """
    __slots__ = ['sequence', 'missedCleavage', 'isUnique', 'proteins',
                 'proteinPositions']

    def __init__(self, sequence, mc=None):
        self.sequence = sequence
        self.missedCleavage = mc
        self.isUnique = None
        self.proteins = set()
        self.proteinPositions = dict()

    def length(self):
        """Returns the number of amino acids of the polypeptide sequence."""
        return len(self.sequence)

    def mass(self):
        """Returns the mass of the polypeptide sequence in Dalton."""
        return maspy.peptidemethods.calcPeptideMass(self.sequence)

    def _reprJSON(self):
        """Returns a JSON serializable represenation of a ``PeptideSequence``
        class instance. Use :func:`maspy.proteindb.PeptideSequence._fromJSON()`
        to generate a new ``PeptideSequence`` instance from the return value.

        :returns: a JSON serializable python object
        """
        return {'__PepSeq__': [self.sequence, self.missedCleavage,
                               self.isUnique, list(self.proteins),
                               self.proteinPositions]}

    @classmethod
    def _fromJSON(cls, jsonobject):
        """Generates a new instance of :class:`maspy.proteindb.PeptideSequence`
        from a decoded JSON object (as generated by
        :func:`maspy.proteindb.PeptideSequence._reprJSON()`).

        :param jsonobject: decoded JSON object

        :returns: a new instance of :class:`PeptideSequence`
        """
        newInstance = cls(jsonobject[0], jsonobject[1])
        newInstance.isUnique = jsonobject[2]
        newInstance.proteins = set(jsonobject[3])
        newInstance.proteinPositions = jsonobject[4]
        return newInstance

    @staticmethod
    def jsonHook(encoded):
        """Custom JSON decoder that allows construction of a new
        ``PeptideSequence`` instance from a decoded JSON object.

        :param encoded: a JSON decoded object literal (a dict)

        :returns: "encoded" or :class:`PeptideSequence`
        """
        if '__PepSeq__' in encoded:
            return PeptideSequence._fromJSON(encoded['__PepSeq__'])
        else:
            return encoded


class ProteinSequence(object):
    """Describes a protein.

    :ivar id: identifier of the protein, for example a uniprot id.
    :ivar name: name of the protein
    :ivar sequence: amino acid sequence of the protein
    :ivar fastaHeader: str(), the proteins faster header line
    :ivar fastaInfo: dict(), the interpreted fasta header as generated when
        using a faster header parsing function, see :func:`fastaParseSgd()`.
    :ivar isUnique: bool, True if at least one unique peptide can be assigned to
        the protein
    :ivar uniquePeptides: a set of peptides which can be unambiguously assigned
        to this protein
    :ivar sharedPeptides: a set of peptides which are shared between different
        proteins
    :ivar coverageUnique: the number of amino acids in the protein sequence that
        are coverd by unique peptides
    :ivar coverageShared: the number of amino acids in the protein sequence that
        are coverd by unique or shared peptides
    """
    def __init__(self, identifier, sequence, name=str()):
        self.id = identifier
        self.name = name
        self.sequence = sequence

        self.isUnique = None
        self.uniquePeptides = set()
        self.sharedPeptides = set()

    def mass(self):
        """Returns the number of amino acids of the polypeptide sequence."""
        return maspy.peptidemethods.calcPeptideMass(self.sequence)

    def length(self):
        """Returns the mass of the polypeptide sequence in dalton."""
        return len(self.sequence)

    def _reprJSON(self):
        """Returns a JSON serializable represenation of a ``ProteinSequence``
        class instance. Use :func:`maspy.proteindb.ProteinSequence._fromJSON()`
        to generate a new ``ProteinSequence`` instance from the return value.

        :returns: a JSON serializable python object
        """

        jsonDict = self.__dict__
        jsonDict['uniquePeptides'] = list(jsonDict['uniquePeptides'])
        jsonDict['sharedPeptides'] = list(jsonDict['sharedPeptides'])
        return {'__ProtSeq__': jsonDict}

    @classmethod
    def _fromJSON(cls, jsonobject):
        """Generates a new instance of :class:`maspy.proteindb.ProteinSequence`
        from a decoded JSON object (as generated by
        :func:`maspy.proteindb.ProteinSequence._reprJSON()`).

        :param jsonobject: decoded JSON object

        :returns: a new instance of :class:`ProteinSequence`
        """
        newInstance = cls(None, None)
        newInstance.__dict__.update(jsonobject)
        newInstance.uniquePeptides = set(newInstance.uniquePeptides)
        newInstance.sharedPeptides = set(newInstance.sharedPeptides)
        return newInstance

    @staticmethod
    def jsonHook(encoded):
        """Custom JSON decoder that allows construction of a new
        ``ProteinSequence`` instance from a decoded JSON object.

        :param encoded: a JSON decoded object literal (a dict)

        :returns: "encoded" or :class:`ProteinSequence`
        """
        if '__ProtSeq__' in encoded:
            return ProteinSequence._fromJSON(encoded['__ProtSeq__'])
        else:
            return encoded


class ProteinDatabase(object):
    """Describes proteins and peptides generated by an in silico digestion of
    proteins.

    :ivar peptides: {sequence:PeptideSequence(), ...} contains elements of
        :class:`PeptideSequence` derived by an in silico digest of the proteins
    :ivar proteins: {proteinId:Protein(), proteinId:Protein()}, used to access
        :class:`ProteinSequence` elements by their id
    :ivar proteinNames: {proteinName:Protein(), proteinName:Protein()},
        alternative way to access :class:`ProteinSequence` elements by their
        names. Must be populated manually
    :ivar info: a dictionary containing information about the protein database
        and parameters specified for the in silico digestion of the protein
        entries. ::

            {'name': str, 'mc': str, 'cleavageRule': str, 'minLength': int
             'maxLength': int, 'ignoreIsoleucine': bool, 'removeNtermM': bool
             }

        **name**: a descriptive name of the protein database, used as the file
            name when saving the protein database to the hard disk
        **mc**: number of allowed missed cleavage sites
        **cleavageRule**: cleavage rule expressed in a regular expression
        **minLength**: minimal peptide length
        **maxLength**: maximal peptide length
        **ignoreIsoleucine**: if True Isoleucine and Leucinge in peptide
            sequences are treated as indistinguishable.
        **removeNtermM**: if True also peptides with the N-terminal Methionine
            of the protein removed are considered.
    """
    def __init__(self):
        self.peptides = dict()
        self.proteins = dict()
        self.proteinNames = dict()
        self.info = {'name': '', 'mc': 0, 'cleavageRule': '', 'minLength': 0,
                     'maxLength': 0, 'ignoreIsoleucine': False,
                     'removeNtermM': False}

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

    def save(self, path, compress=True):
        """Writes the ``.proteins`` and ``.peptides`` entries to the hard disk
        as a ``proteindb`` file.

        .. note::
            If ``.save()`` is called and no ``proteindb`` file is present in the
            specified path a new files is generated, otherwise the old file is
            replaced.

        :param path: filedirectory to which the ``proteindb`` file is written.
            The output file name is specified by ``self.info['name']``
        :param compress: bool, True to use zip file compression
        """
        with aux.PartiallySafeReplace() as msr:
            filename = self.info['name'] + '.proteindb'
            filepath = aux.joinpath(path, filename)
            with msr.open(filepath, mode='w+b') as openfile:
                self._writeContainer(openfile, compress=compress)

    def _writeContainer(self, filelike, compress=True):
        """Writes the ``.proteins`` and ``.peptides`` entries to the
        ``proteindb`` format. In addition it also dumps the ``self.info`` entry
        to the zipfile with the filename ``info``. For details see
        :func:`maspy.auxiliary.writeJsonZipfile()`

        :param filelike: path to a file (str) or a file-like object
        :param compress: bool, True to use zip file compression
        """
        aux.writeJsonZipfile(filelike, self.proteins, compress, 'w', 'proteins')
        aux.writeJsonZipfile(filelike, self.peptides, compress, 'a', 'peptides')
        zipcomp = zipfile.ZIP_DEFLATED if compress else zipfile.ZIP_STORED
        with zipfile.ZipFile(filelike, 'a', allowZip64=True) as containerFile:
            infodata = {key: value for key, value in
                        viewitems(self.info) if key != 'path'
                        }
            containerFile.writestr('info', json.dumps(infodata, zipcomp))

    @classmethod
    def load(cls, path, name):
        """Imports the specified ``proteindb`` file from the hard disk.

        :param path: filedirectory of the ``proteindb`` file
        :param name: filename without the file extension ".proteindb"

        .. note:: this generates rather large files, which actually take longer
            to import than to newly generate. Maybe saving / loading should be
            limited to the protein database whitout in silico digestion
            information.
        """

        filepath = aux.joinpath(path, name + '.proteindb')
        with zipfile.ZipFile(filepath, 'r', allowZip64=True) as containerZip:
            #Convert the zipfile data into a str object, necessary since
            #containerZip.read() returns a bytes object.
            proteinsString = io.TextIOWrapper(containerZip.open('proteins'),
                                              encoding='utf-8'
                                              ).read()
            peptidesString = io.TextIOWrapper(containerZip.open('peptides'),
                                              encoding='utf-8'
                                              ).read()
            infoString = io.TextIOWrapper(containerZip.open('info'),
                                          encoding='utf-8'
                                          ).read()
        newInstance = cls()
        newInstance.proteins = json.loads(proteinsString,
                                          object_hook=ProteinSequence.jsonHook)
        newInstance.peptides = json.loads(peptidesString,
                                          object_hook=PeptideSequence.jsonHook)
        newInstance.info.update(json.loads(infoString))
        return newInstance

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


# --- import for ProteinDatabase class --- #
def importProteinDatabase(filePath, proteindb=None, decoyTag='[decoy]',
        contaminationTag='[cont]', headerParser=None, forceId=False,
        cleavageRule='[KR]', minLength=5, maxLength=40, missedCleavage=2,
        ignoreIsoleucine=False, removeNtermM=True):
    """Generates a :class:`ProteinDatabase` by in silico digestion of proteins
    from a fasta file.

    :param filePath: File path
    :param proteindb: optional an existing :class:`ProteinDatabase` can be
        specified, otherwise a new instance is generated and returned
    :param decoyTag: If a fasta file contains decoy protein entries, they should
        be specified with a sequence tag
    :param contaminationTag: If a fasta file contains contamination protein
        entries, they should be specified with a sequence tag
    :param headerParser: optional a headerParser can be specified
        #TODO: describe how a parser looks like
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
    proteindb = ProteinDatabase() if proteindb is None else proteindb
    fastaRead = _readFastaFile(filePath)

    for header, sequence in fastaRead:
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
            proteinName = ''.join(itertools.chain(proteinTags,
                                                  [headerInfo['name']]
                                                  )
                                  )
        else:
            proteinName = proteinId

        if proteinId not in proteindb.proteins:
            protein = ProteinSequence(proteinId, sequence)
            protein.name = proteinName
            protein.fastaHeader = header
            protein.fastaInfo = headerInfo
            proteindb.proteins[protein.id] = protein

        #Perform the insilico digestion
        _digestion = maspy.peptidemethods.digestInSilico(sequence, cleavageRule,
                                                         missedCleavage,
                                                         removeNtermM,
                                                         minLength, maxLength
                                                         )

        #Add peptides to the protein database
        for unmodPeptide, info in _digestion:
            if ignoreIsoleucine:
                unmodPeptideNoIsoleucine = unmodPeptide.replace('I', 'L')
                if unmodPeptideNoIsoleucine in proteindb.peptides:
                    currPeptide = proteindb.peptides[unmodPeptideNoIsoleucine]
                else:
                    currPeptide = PeptideSequence(unmodPeptideNoIsoleucine,
                                                  mc=info['missedCleavage']
                                                  )
                    proteindb.peptides[unmodPeptideNoIsoleucine] = currPeptide

                if unmodPeptide not in proteindb.peptides:
                    proteindb.peptides[unmodPeptide] = currPeptide
            else:
                if unmodPeptide in proteindb.peptides:
                    currPeptide = proteindb.peptides[unmodPeptide]
                else:
                    currPeptide = PeptideSequence(unmodPeptide,
                                                  mc=info['missedCleavage']
                                                  )
                    proteindb.peptides[unmodPeptide] = currPeptide

            if proteinId not in currPeptide.proteins:
                currPeptide.proteins.add(proteinId)
                #TODO: change that a peptide can appear multiple times in a
                #  protein sequence.
                currPeptide.proteinPositions[proteinId] = (info['startPos'],
                                                           info['endPos']
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
    with io.open(filepath) as openfile:
        #Iterate through lines until the first header is encountered
        try:
            line = next(openfile)
            while line[0] != '>':
                line = next(openfile)
            header = processHeaderLine(line)
            sequences = list()
        except StopIteration:
            errorText = 'File does not contain fasta entries.'
            raise maspy.errors.FileFormatError(errorText)

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


def _extractFastaHeader(fastaHeader, parser=None, forceId=False):
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

    :returns: dict, describing a fasta header
    """
    if parser is None:
        try:
            headerInfo = pyteomics.fasta.parse(fastaHeader)
        except pyteomics.auxiliary.PyteomicsError as pyteomicsError:
            #If forceId is set True, it uses the whole header as an id
            if forceId:
                headerInfo = {'id': fastaHeader}
            else:
                raise pyteomicsError
    else:
        headerInfo = parser(fastaHeader)
    return headerInfo


def fastaParseSgd(header):
    """Custom parser for fasta headers in the SGD format, see
    www.yeastgenome.org.

    :param header: str, protein entry header from a fasta file

    :returns: dict, parsed header
    """
    rePattern = '([\S]+)\s([\S]+).+(\".+\")'
    ID, name, description = re.match(rePattern, header).groups()
    info = {'id':ID, 'name':name, 'description':description}
    return info
