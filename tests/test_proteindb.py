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

import os
import sys
import tempfile

import pyteomics.auxiliary
import unittest

sys.path.append(os.path.abspath('..'))
import maspy.errors
import maspy._proteindb_refactoring as MODULE


class SetupFastaTestFiles(unittest.TestCase):
    def setUp(self):
        uniprotString = '\n'.join([
            '>sp|ID001|geneId1_taxon Protein description 1 OS=organism GN=gene_name1 PE=1 SV=1',
            'MERCGWV ',
            'HQFDPVK ',
            'NHQPQVT',
            'CYPGNKP*',
            '>sp|ID002|geneId2_taxon Protein description 2 OS=organism GN=gene_name2 PE=1 SV=1',
            'MYTLNWQPPYDWSWMLGFLAARAVSSVETVADSYYARSLAVGEYRGVVTAIPDIARHTLH',
            '>sp|ID003|geneId3_taxon Protein description 3 OS=organism GN=gene_name2 PE=1 SV=1',
        ])
        uniprotHeaders = [
            'sp|ID001|geneId1_taxon Protein description 1 OS=organism GN=gene_name1 PE=1 SV=1',
            'sp|ID002|geneId2_taxon Protein description 2 OS=organism GN=gene_name2 PE=1 SV=1'
        ]
        uniprotFastaSequences = [
            'MERCGWVHQFDPVKNHQPQVTCYPGNKP',
            'MYTLNWQPPYDWSWMLGFLAARAVSSVETVADSYYARSLAVGEYRGVVTAIPDIARHTLH'
        ]
        uniprotHeaderInfo = [
            {'GN': 'gene_name1', 'OS': 'organism', 'PE': 1, 'SV': 1, 'db': 'sp',
             'entry': 'geneId1_taxon', 'gene_id': 'geneId1', 'id': 'ID001',
             'name': 'Protein description 1', 'taxon': 'taxon'},
            {'GN': 'gene_name2', 'OS': 'organism', 'PE': 1, 'SV': 1, 'db': 'sp',
             'entry': 'geneId2_taxon', 'gene_id': 'geneId2', 'id': 'ID002',
             'name': 'Protein description 2', 'taxon': 'taxon'}
        ]
        sgdHeaders = [
            'ID001 NAME001 SGDID:S000000001, Chr I from 10-11, Verified ORF, "Protein description 1"',
            'ID002-A NAME002 SGDID:S000000002, Chr XVI from 10-20,30-50, "Protein description 2"',
        ]
        sgdHeaderInfo = [
            {'id': 'ID001', 'name': 'NAME001', 'description': '"Protein description 1"'},
            {'id': 'ID002-A', 'name': 'NAME002', 'description': '"Protein description 2"'}
        ]

        proteinEntry = MODULE.ProteinEntry(
            uniprotHeaderInfo[0]['id'], uniprotHeaderInfo[0]['name'],
            uniprotFastaSequences[0], uniprotHeaders[0], uniprotHeaderInfo[0],
            isDecoy=False, isContaminant=False
        )

        with tempfile.NamedTemporaryFile(delete=False, mode='w') as uniprotTempfile:
            uniprotTempfile.file.write(uniprotString)
        with tempfile.NamedTemporaryFile(delete=False, mode='w') as noFastaTempfile:
            noFastaTempfile.file.write('\n'.join(uniprotFastaSequences))

        self.uniprotHeaders = uniprotHeaders
        self.uniprotHeaderInfo = uniprotHeaderInfo
        self.uniprotFastaSequences = uniprotFastaSequences

        self.sgdHeaders = sgdHeaders
        self.sgdHeaderInfo = sgdHeaderInfo

        self.uniprotTempfile = uniprotTempfile
        self.noFastaTempfile = noFastaTempfile

        self.proteinEntry = proteinEntry

    def tearDown(self):
        # Remove the directory after the test
        os.remove(self.uniprotTempfile.name)
        os.remove(self.noFastaTempfile.name)
        self.assertFalse(os.path.exists(self.uniprotTempfile.name))
        self.assertFalse(os.path.exists(self.noFastaTempfile.name))


class TestBasicProteindbFunctions(SetupFastaTestFiles):
    #Tests for _readFastaFile()
    def test_readUniprotFastaFile(self):
        fastaEntryIter = MODULE._readFastaFile(self.uniprotTempfile.name)
        for i, (header, sequence) in enumerate(fastaEntryIter):
            self.assertEqual(header, self.uniprotHeaders[i])
            self.assertEqual(sequence, self.uniprotFastaSequences[i])

    def test_readEmptyFastaFile(self):
        fastaEntryIter = MODULE._readFastaFile(self.noFastaTempfile.name)
        with self.assertRaises(maspy.errors.FileFormatError):
            next(fastaEntryIter)

    #Test fastaParseSgd()
    def test_fastaParseSgd(self):
        for headerString, expectedHeaderInfo in zip(self.sgdHeaders,
                                                    self.sgdHeaderInfo):
            headerInfo = MODULE.fastaParseSgd(headerString)
            self.assertDictEqual(headerInfo, expectedHeaderInfo)

    #Test fastaParserSpectraClusterPy()
    def test_fastaParserSpectraClusterPy(self):
        headers = [
            'ID_001|geneId1_taxon Protein description 1',
            'ID_001 geneId1_taxon Protein description 1',
            'ID_001 geneId1_taxon|Protein|description 1',
            'ID_001',
        ]
        for headerString in headers:
            headerInfo = MODULE.fastaParserSpectraClusterPy(headerString)
            self.assertEqual(headerInfo['id'], 'ID_001')
        #Entry contains a '|' but not a ' ' symbol


    #Test _parseFastaHeader()
    def test_parseFastaHeader_uniprotHeader(self):
        for headerString, expectedHeaderInfo in zip(self.uniprotHeaders,
                                                    self.uniprotHeaderInfo):
            headerInfo = MODULE._parseFastaHeader(headerString)
            self.assertDictEqual(headerInfo, expectedHeaderInfo)

    def test_parseFastaHeader_unparsableRaisesException(self):
        with self.assertRaises(pyteomics.auxiliary.PyteomicsError):
            MODULE._parseFastaHeader(self.sgdHeaders[0])

    def test_parseFastaHeader_forceId(self):
        for headerString in self.sgdHeaders:
            headerInfo = MODULE._parseFastaHeader(headerString, forceId=True)
            self.assertDictEqual(headerInfo, {'id': headerString})

    def test_parseFastaHeader_sgdHeader(self):
        for headerString, expectedHeaderInfo in zip(self.sgdHeaders,
                                                    self.sgdHeaderInfo):
            headerInfo = MODULE._parseFastaHeader(headerString,
                                                  parser=MODULE.fastaParseSgd)
            self.assertDictEqual(headerInfo, expectedHeaderInfo)

    #Test processing functions
    def test_removeHeaderTag(self):
        expectedTag = '[rev]'
        expectedHeader = self.uniprotHeaders[0]

        combinedHeader = ''.join([expectedTag, expectedHeader])
        header, tagIsPresent = MODULE._removeHeaderTag(combinedHeader, expectedTag)
        self.assertEqual(header, expectedHeader)
        self.assertTrue(tagIsPresent)

        header, tagIsPresent = MODULE._removeHeaderTag(expectedHeader, expectedTag)
        self.assertEqual(header, expectedHeader)
        self.assertFalse(tagIsPresent)

    def test_idFromHeaderInfo(self):
        headerInfo = {'id': 'ID001'}
        decoyTag = '[decoy]'

        proteinId = MODULE._idFromHeaderInfo(headerInfo, True, decoyTag)
        self.assertEqual(proteinId, '[decoy]ID001')

        proteinId = MODULE._idFromHeaderInfo(headerInfo, False, decoyTag)
        self.assertEqual(proteinId, 'ID001')

    def test_nameFromHeaderInfo(self):
        headerInfo = {'id': 'ID001', 'name': 'NAME001'}
        decoyTag = '[decoy]'

        proteinId = MODULE._nameFromHeaderInfo(headerInfo, True, decoyTag)
        self.assertEqual(proteinId, '[decoy]NAME001')

        proteinId = MODULE._nameFromHeaderInfo(headerInfo, False, decoyTag)
        self.assertEqual(proteinId, 'NAME001')

        headerInfo = {'id': 'ID001'}
        decoyTag = '[decoy]'

        proteinId = MODULE._nameFromHeaderInfo(headerInfo, True, decoyTag)
        self.assertEqual(proteinId, '[decoy]ID001')

        proteinId = MODULE._nameFromHeaderInfo(headerInfo, False, decoyTag)
        self.assertEqual(proteinId, 'ID001')

    def test_proteinTagPresent(self):
        tag = '[decoy]'
        for header in self.uniprotHeaders:
            self.assertFalse(MODULE._proteinTagPresent(header, tag))
        decoyHeader = self.uniprotHeaders[0].replace('ID001', tag+'ID001')
        self.assertTrue(MODULE._proteinTagPresent(decoyHeader, tag))


class TestProteinEntryClass(SetupFastaTestFiles):
    def test_calculateLength(self):
        protein = self.proteinEntry
        self.assertEqual(protein.length(), len(protein.sequence))

    def test_calculateMass(self):
        protein = self.proteinEntry
        self.assertAlmostEqual(protein.mass(), 3294.52736188385)


class TestPeptideEntryClass(SetupFastaTestFiles):
    def test_calculateLength(self):
        peptide = MODULE.PeptideEntry(self.uniprotFastaSequences[0])
        self.assertEqual(peptide.length(), len(peptide.sequence))

    def test_calculateMass(self):
        peptide = MODULE.PeptideEntry(self.uniprotFastaSequences[0])
        self.assertAlmostEqual(peptide.mass(), 3294.52736188385)


class TestProteindbClass(SetupFastaTestFiles):
    def test_addProtein(self):
        proteindb = MODULE.ProteinDatabase()
        protein = self.proteinEntry
        proteindb._addProtein(
            protein.id, protein.name, protein.sequence, protein.fastaHeader,
            protein.headerInfo, isDecoy=protein.isDecoy,
            isContaminant=protein.isContaminant
        )
        self.assertIn(protein.id, proteindb.proteins)
        self.assertEqual(proteindb.proteins[protein.id].sequence, protein.sequence)

    def test_getStdSequence(self):
        proteindb_ignoreTrue = MODULE.ProteinDatabase(ignoreIsoleucine=True)
        proteindb_ignoreFalse = MODULE.ProteinDatabase(ignoreIsoleucine=False)
        sequence = 'PEPTIDEK'
        stdSequence = 'PEPTLDEK'
        self.assertEqual(proteindb_ignoreTrue.getStdSequence(sequence), stdSequence)
        self.assertEqual(proteindb_ignoreFalse.getStdSequence(sequence), sequence)

    def test_addPeptide(self):
        sequence = 'PEPTIDEK'
        stdSequence = 'PEPTLDEK'
        proteinId = self.proteinEntry.id
        info = {'missedCleavage': 0, 'startPos': 1, 'endPos': 10}

        # Without ignoreIsoleucine
        proteindb = MODULE.ProteinDatabase(ignoreIsoleucine=False)
        proteindb.proteins[self.proteinEntry.id] = self.proteinEntry
        proteindb._addPeptide(sequence, proteinId, info)

        self.assertIn(sequence, proteindb.peptides)
        self.assertNotIn(stdSequence, proteindb.peptides)

        self.assertIn(sequence, proteindb.proteins[proteinId].peptides)
        self.assertNotIn(stdSequence, proteindb.proteins[proteinId].peptides)

        self.assertEqual(proteindb.peptides[sequence].sequence, sequence)
        self.assertEqual(proteindb.peptides[sequence].missedCleavage,info['missedCleavage'])
        self.assertIn(proteinId, proteindb.peptides[sequence].proteins)
        self.assertDictEqual({proteinId: (1, 10)},
                             proteindb.peptides[sequence].proteinPositions
                             )

        # With ignoreIsoleucine
        proteindb = MODULE.ProteinDatabase(ignoreIsoleucine=True)
        proteindb.proteins[self.proteinEntry.id] = self.proteinEntry
        proteindb._addPeptide(sequence, proteinId, info)

        self.assertIn(sequence, proteindb.peptides)
        self.assertIn(stdSequence, proteindb.peptides)

        self.assertIn(sequence, proteindb.proteins[proteinId].peptides)
        self.assertNotIn(stdSequence, proteindb.proteins[proteinId].peptides)

        self.assertEqual(proteindb.peptides[stdSequence].sequence, stdSequence)
        self.assertEqual(proteindb.peptides[sequence].sequence, stdSequence)
        self.assertEqual(proteindb.peptides[stdSequence].missedCleavage,info['missedCleavage'])
        self.assertIn(proteinId, proteindb.peptides[stdSequence].proteins)
        self.assertDictEqual({proteinId: (1, 10)},
                             proteindb.peptides[stdSequence].proteinPositions
                             )







