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

from collections import defaultdict as ddict
import numpy
import unittest
import maspy.inference as MODULE


class TestSimpleProteinInferenceFunctions(unittest.TestCase):
    def setUp(self):
        proteinToPeptides = {
            '01_P1': {'01_p1', '01_p2'},
            '02_P1': {'02_p1', '02_p2'},
            '02_P2': {'02_p2'},
            '03_P1': {'03_p1', '03_p2'},
            '03_P2': {'03_p1', '03_p2'},
            '04_P1': {'04_p1', '04_p2'},
            '04_P2': {'04_p1', '04_p2'},
            '04_P3': {'04_p2'},
            }
        peptideToProteins = {
            '01_p1': {'01_P1'},
            '01_p2': {'01_P1'},
            '02_p1': {'02_P1'},
            '02_p2': {'02_P1', '02_P2'},
            '03_p1': {'03_P1', '03_P2'},
            '03_p2': {'03_P1', '03_P2'},
            '04_p1': {'04_P1', '04_P2'},
            '04_p2': {'04_P1', '04_P2', '04_P3'}
            }
        mergedProteinToPeptides = {
            '01_P1': {'01_p1', '01_p2'},
            '02_P1': {'02_p1', '02_p2'},
            '02_P2': {'02_p2'},
            tuple(['03_P1', '03_P2']): {'03_p1', '03_p2'},
            tuple(['04_P1', '04_P2']): {'04_p1', '04_p2'},
            '04_P3': {'04_p2'},
            }
        mergedPeptideToProteins = {
            '01_p1': {'01_P1'},
            '01_p2': {'01_P1'},
            '02_p1': {'02_P1'},
            '02_p2': {'02_P1', '02_P2'},
            '03_p1': {tuple(['03_P1', '03_P2'])},
            '03_p2': {tuple(['03_P1', '03_P2'])},
            '04_p1': {tuple(['04_P1', '04_P2'])},
            '04_p2': {tuple(['04_P1', '04_P2']), '04_P3'}
            }

        observedPeptides = set(peptideToProteins)
        expectedGroups = {p.split('_')[0] for p in proteinToPeptides}
        equalProteins = [{'03_P1', '03_P2'}, {'04_P1', '04_P2'}]

        uniqueProteins = {'01_P1', '02_P1'}
        uniqueMergedProteins = {'01_P1', '02_P1', tuple(['03_P1', '03_P2']),
                                tuple(['04_P1', '04_P2'])
                                }

        subsetProteins = [
            ('02_P2', {'02_P1'}), ('03_P1', {'03_P2'}), ('03_P2', {'03_P1'}),
            ('04_P1', {'04_P2'}), ('04_P2', {'04_P1'}),
            ('04_P3', {'04_P1', '04_P2'})
            ]
        subsetMergedProteins = [
            ('02_P2', {'02_P1'}), ('04_P3', {tuple(['04_P1', '04_P2'])})
            ]

        self.proteinToPeptides = proteinToPeptides
        self.mergedProteinToPeptides = mergedProteinToPeptides
        self.peptideToProteins = peptideToProteins
        self.mergedPeptideToProteins = mergedPeptideToProteins
        self.observedPeptides = observedPeptides

        self.expectedGroups = expectedGroups
        self.equalProteins = equalProteins
        self.equalProteinTuples = sorted([tuple(p) for p in equalProteins])

        self.uniqueProteins = uniqueProteins
        self.uniqueMergedProteins = uniqueMergedProteins

        self.subsetProteins = sorted(subsetProteins)
        self.subsetMergedProteins = sorted(subsetMergedProteins)

    def test_groupConnectedPeptides(self):
        peptideGroups = MODULE._groupConnectedPeptides(self.peptideToProteins,
                                                       self.proteinToPeptides)
        self.assertEqual(len(peptideGroups), len(self.expectedGroups))

        allGroupedPeptides = reduce(lambda s1, s2: s1.union(s2), peptideGroups)
        self.assertSetEqual(allGroupedPeptides, self.observedPeptides)

        #Assert that all grouped peptides are actually from the same group
        for group in peptideGroups:
            self.assertIsInstance(group, set)
            self.assertEqual(len({p.split('_')[0] for p in group}), 1)

    def test_findEqualEvidenceProteins(self):
        equalProteins = MODULE._findEqualEvidenceProteins(self.proteinToPeptides,
                                                          self.proteinToPeptides)
        equalProteinTuples = sorted([tuple(p) for p in equalProteins])
        self.assertListEqual(equalProteinTuples, self.equalProteinTuples)

    def test_mergeProteinEntries(self):
        mergedProteinToPeptides = MODULE._mergeProteinEntries(self.equalProteinTuples,
                                                              self.proteinToPeptides)
        self.assertDictEqual(mergedProteinToPeptides, self.mergedProteinToPeptides)

    def test_invertMapping(self):
        peptideToProteins = MODULE._invertMapping(self.proteinToPeptides)
        self.assertDictEqual(peptideToProteins, self.peptideToProteins)

        proteinToPeptides = MODULE._invertMapping(self.peptideToProteins)
        self.assertDictEqual(proteinToPeptides, self.proteinToPeptides)

    def test_reducedProtToPeps(self):
        proteins = ['02_P1', '03_P1']
        reduceProteinToPeptides = MODULE._reducedProtToPeps(self.proteinToPeptides, proteins)
        for protein in proteins:
            self.assertNotIn(protein, reduceProteinToPeptides)
            self.assertIn(protein, self.proteinToPeptides)

    def test_findUniqueProteins(self):
        uniqueProteins = MODULE._findUniqueProteins(self.peptideToProteins)
        self.assertSetEqual(uniqueProteins, self.uniqueProteins)
        uniqueMergedProteins = MODULE._findUniqueProteins(self.mergedPeptideToProteins)
        self.assertSetEqual(uniqueMergedProteins, self.uniqueMergedProteins)

    def test_findSubsetProteins(self):
        subsetProteins = MODULE._findSubsetProteins(self.proteinToPeptides,
                                                    self.proteinToPeptides,
                                                    self.peptideToProteins)
        self.assertListEqual(sorted(subsetProteins), self.subsetProteins)
        subsetMergedProteins = MODULE._findSubsetProteins(self.mergedProteinToPeptides,
                                                          self.mergedProteinToPeptides,
                                                          self.mergedPeptideToProteins)
        self.assertListEqual(sorted(subsetMergedProteins), self.subsetMergedProteins)

    def test_getValueCounts(self):
        peptideFrequency = MODULE._getValueCounts(self.peptideToProteins)
        for peptide, counts in viewitems(peptideFrequency):
            self.assertEqual(counts, len(self.peptideToProteins[peptide]))

        protPepCounts = MODULE._getValueCounts(self.proteinToPeptides)
        for protein, counts in viewitems(protPepCounts):
            self.assertEqual(counts, len(self.proteinToPeptides[protein]))


class TestFindRedundantProteins(unittest.TestCase):
    def setUp(self):
        proteinToPeptides = {
            '05_P1': {'05_p1', '05_p2', '05_p3'},
            '05_P2': {'05_p1', '05_p4'},
            '05_P3': {'05_p2', '05_p3', '05_p4'},
            '06_P1': {'06_p1', '06_p2', '06_p3'},
            '06_P2': {'06_p2', '06_p3'},
            '06_P3': {'06_p2', '06_p4'},
            '06_P4': {'06_p2', '06_p3', '06_p4'},
            '06_P5': {'06_p2', '06_p4'},
            '07_P1': {'07_p1', '07_p2'},
            '07_P2': {'07_p1', '07_p3', '07_p4'},
            '07_P3': {'07_p2', '07_p3'},
            '07_P4': {'07_p3', '07_p5'}
            }
        peptideToProteins = {
            '05_p1': {'05_P1', '05_P2'},
            '05_p2': {'05_P1', '05_P3'},
            '05_p3': {'05_P1', '05_P3'},
            '05_p4': {'05_P2', '05_P3'},
            '06_p1': {'06_P1'},
            '06_p2': {'06_P1', '06_P2' , '06_P3', '06_P4', '06_P5'},
            '06_p3': {'06_P1', '06_P2' , '06_P4'},
            '06_p4': {'06_P3', '06_P4', '06_P5'},
            '07_p1': {'07_P1', '07_P2'},
            '07_p2': {'07_P1', '07_P3'},
            '07_p3': {'07_P2', '07_P3', '07_P4'},
            '07_p4': {'07_P2'},
            '07_p5': {'07_P4'},
            }
        redundantProteins = {'05_P2', '06_P2', '06_P3', '06_P5', '07_P3'}

        self.proteinToPeptides = proteinToPeptides
        self.peptideToProteins = peptideToProteins
        self.redundantProteins = redundantProteins

    def test_findRedundantProteins(self):
        redundantProteins = MODULE._findRedundantProteins(self.proteinToPeptides,
                                                          self.peptideToProteins)
        self.assertSetEqual(redundantProteins, self.redundantProteins)
        #The group "07" includes two non unique proteins with an equal number of
        #peptides: "07_P1", "07_P3". The protein with the higher peptide
        #frequencies should be removed; "07_P3".

        #Remove the redundant proteins from protToPeps and pepToProts, assert
        #that all proteins are unique and each peptide is still explained by at
        #least one peptide.
        protToPeps = MODULE._reducedProtToPeps(self.proteinToPeptides, redundantProteins)
        pepToProts = MODULE._invertMapping(protToPeps)
        uniqueProteins = MODULE._findUniqueProteins(pepToProts)
        self.assertSetEqual(uniqueProteins, set(protToPeps))

        for peptide, proteins in viewitems(pepToProts):
            self.assertGreater(len(proteins), 0)
