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

import sys

from collections import defaultdict as ddict
import numpy
import unittest
import maspy.inference as MODULE


class SetupMappingExamples(unittest.TestCase):
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
            '07_P4': {'07_p3', '07_p5'},
            '08_P1': {'08_p1', '08_p2'},
            '08_P2': {'08_p3', '08_p4'},
            '08_P3': {'08_p2', '08_p3'},
            }
        peptideToProteins = {
            '01_p1': {'01_P1'},
            '01_p2': {'01_P1'},
            '02_p1': {'02_P1'},
            '02_p2': {'02_P1', '02_P2'},
            '03_p1': {'03_P1', '03_P2'},
            '03_p2': {'03_P1', '03_P2'},
            '04_p1': {'04_P1', '04_P2'},
            '04_p2': {'04_P1', '04_P2', '04_P3'},
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
            '08_p1': {'08_P1'},
            '08_p2': {'08_P1', '08_P3'},
            '08_p3': {'08_P2', '08_P3'},
            '08_p4': {'08_P2'},
            }
        #Expected protein to peptide mapping after merging "samesetProteins"
        mergedProteinToPeptides = {
            '01_P1': {'01_p1', '01_p2'},
            '02_P1': {'02_p1', '02_p2'},
            '02_P2': {'02_p2'},
            tuple(['03_P1', '03_P2']): {'03_p1', '03_p2'},
            tuple(['04_P1', '04_P2']): {'04_p1', '04_p2'},
            '04_P3': {'04_p2'},
            '05_P1': {'05_p1', '05_p2', '05_p3'},
            '05_P2': {'05_p1', '05_p4'},
            '05_P3': {'05_p2', '05_p3', '05_p4'},
            '06_P1': {'06_p1', '06_p2', '06_p3'},
            '06_P2': {'06_p2', '06_p3'},
            tuple(['06_P3', '06_P5']): {'06_p2', '06_p4'},
            '06_P4': {'06_p2', '06_p3', '06_p4'},
            '07_P1': {'07_p1', '07_p2'},
            '07_P2': {'07_p1', '07_p3', '07_p4'},
            '07_P3': {'07_p2', '07_p3'},
            '07_P4': {'07_p3', '07_p5'},
            '08_P1': {'08_p1', '08_p2'},
            '08_P2': {'08_p3', '08_p4'},
            '08_P3': {'08_p2', '08_p3'},
        }
        #Expected peptide to protein mapping after merging "samesetProteins"
        mergedPeptideToProteins = {
            '01_p1': {'01_P1'},
            '01_p2': {'01_P1'},
            '02_p1': {'02_P1'},
            '02_p2': {'02_P1', '02_P2'},
            '03_p1': {tuple(['03_P1', '03_P2'])},
            '03_p2': {tuple(['03_P1', '03_P2'])},
            '04_p1': {tuple(['04_P1', '04_P2'])},
            '04_p2': {tuple(['04_P1', '04_P2']), '04_P3'},
            '05_p1': {'05_P1', '05_P2'},
            '05_p2': {'05_P1', '05_P3'},
            '05_p3': {'05_P1', '05_P3'},
            '05_p4': {'05_P2', '05_P3'},
            '06_p1': {'06_P1'},
            '06_p2': {'06_P1', '06_P2' , '06_P4', tuple(['06_P3', '06_P5'])},
            '06_p3': {'06_P1', '06_P2' , '06_P4'},
            '06_p4': {'06_P4', tuple(['06_P3', '06_P5'])},
            '07_p1': {'07_P1', '07_P2'},
            '07_p2': {'07_P1', '07_P3'},
            '07_p3': {'07_P2', '07_P3', '07_P4'},
            '07_p4': {'07_P2'},
            '07_p5': {'07_P4'},
            '08_p1': {'08_P1'},
            '08_p2': {'08_P1', '08_P3'},
            '08_p3': {'08_P2', '08_P3'},
            '08_p4': {'08_P2'},
            }
        #All peptides that are present in the protein and peptide mappings
        observedPeptides = set(peptideToProteins)
        #All proteins that are present in the protein and peptide mappings
        potentialProteins = set(proteinToPeptides)
        #Expected number of "clusters" (ie. groups of connected proteins)
        expectedClusters = {p.split('_')[0] for p in proteinToPeptides}
        #All "samesetProteins"
        equalProteins = [('03_P1', '03_P2'), ('04_P1', '04_P2'),
                         ('06_P3', '06_P5')]
        #All "uniqueProtein" entries, ie. proteins with at least one peptide that
        #is only connected to one single protein entry.
        uniqueProteins = {'01_P1', '02_P1', '06_P1', '07_P2', '07_P4',
                          '08_P1', '08_P2'}
        #All unique protein entries after merging "samesetProteins"
        uniqueMergedProteins = {'01_P1', '02_P1', tuple(['03_P1', '03_P2']),
                                tuple(['04_P1', '04_P2']), '06_P1', '07_P2',
                                '07_P4', '08_P1', '08_P2'
                                }
        #All expected "subsetProteins"
        subsetProteins = [
            ('02_P2', {'02_P1'}), ('04_P3', {'04_P1', '04_P2'}),
            ('06_P2', {'06_P1', '06_P4'}), ('06_P3', {'06_P4'}),
            ('06_P5', {'06_P4'})
        ]
        #All expected "subsetProteins" after merging "samesetProteins"
        subsetMergedProteins = [
            ('02_P2', {'02_P1'}), ('04_P3', {tuple(['04_P1', '04_P2'])}),
            ('06_P2', {'06_P1', '06_P4'}),
            (tuple(['06_P3', '06_P5']), {'06_P4'})
        ]
        #All proteins expected to be defined as "redundantProtein" by the
        #function inference._findRedundantProteins(), for details see docstring
        redundantProteins = {'02_P2', '03_P2', '04_P2', '04_P3',
                             '05_P2', '06_P2', '06_P3', '06_P5',
                             '07_P3', '08_P3'}

        ########################################################################
        # This section describes the expected results of the protein grouping  #
        # procedure performed by the function inference.mappingBasedGrouping() #
        ########################################################################

        #Expected grouping results on the proteinGroup level
        groupingGroupResults = {
            '01_P1': {
                'representative': '01_P1',
                'leading': {'01_P1'},
                'subset': set(),
                'subsumable': set(),
            },
            '02_P1': {
                'representative': '02_P1',
                'leading': {'02_P1'},
                'subset': {'02_P2'},
                'subsumable': set(),
            },
            '03_P1': {
                'representative': '03_P1',
                'leading': {'03_P1', '03_P2'},
                'subset': set(),
                'subsumable': set(),
            },
            '04_P1': {
                'representative': '04_P1',
                'leading': {'04_P1', '04_P2'},
                'subset': {'04_P3'},
                'subsumable': set(),
            },
            '05_P1': {
                'representative': '05_P1',
                'leading': {'05_P1'},
                'subset': set(),
                'subsumable': {'05_P2'},
            },
            '05_P3': {
                'representative': '05_P3',
                'leading': {'05_P3'},
                'subset': set(),
                'subsumable': {'05_P2'},
            },
            '06_P1': {
                'representative': '06_P1',
                'leading': {'06_P1'},
                'subset': {'06_P2'},
                'subsumable': set(),
            },
            '06_P4': {
                'representative': '06_P4',
                'leading': {'06_P4'},
                'subset': {'06_P2', '06_P3', '06_P5'},
                'subsumable': set(),
            },
            '07_P1': {
                'representative': '07_P1',
                'leading': {'07_P1'},
                'subset': set(),
                'subsumable': {'07_P3'},
            },
            '07_P2': {
                'representative': '07_P2',
                'leading': {'07_P2'},
                'subset': set(),
                'subsumable': {'07_P3'},
            },
            '07_P4': {
                'representative': '07_P4',
                'leading': {'07_P4'},
                'subset': set(),
                'subsumable': {'07_P3'},
            },
            '08_P1': {
                'representative': '08_P1',
                'leading': {'08_P1'},
                'subset': set(),
                'subsumable': {'08_P3'},
            },
            '08_P2': {
                'representative': '08_P2',
                'leading': {'08_P2'},
                'subset': set(),
                'subsumable': {'08_P3'},
            },
        }
        for groupingResults in listvalues(groupingGroupResults):
            groupingResults['proteins'] = set(groupingResults['representative'])
            for key in ['leading', 'subset']:
                groupingResults['proteins'].update(groupingResults[key])

        #Expected grouping results on the protein level
        groupingProteinResults = {
            '01_P1': {'subsetOf': set(), 'sameset': set(), 'isSubsumable': False},
            '02_P1': {'subsetOf': set(), 'sameset': set(), 'isSubsumable': False},
            '02_P2': {'subsetOf': {'02_P1'}, 'sameset': set(), 'isSubsumable': False},
            '03_P1': {'subsetOf': set(), 'sameset': {'03_P1', '03_P2'}, 'isSubsumable': False},
            '03_P2': {'subsetOf': set(), 'sameset': {'03_P1', '03_P2'}, 'isSubsumable': False},
            '04_P1': {'subsetOf': set(), 'sameset': {'04_P1', '04_P2'}, 'isSubsumable': False},
            '04_P2': {'subsetOf': set(), 'sameset': {'04_P1', '04_P2'}, 'isSubsumable': False},
            '04_P3': {'subsetOf': {'04_P1', '04_P2'}, 'sameset': set(), 'isSubsumable': False},
            '05_P1': {'subsetOf': set(), 'sameset': set(), 'isSubsumable': False},
            '05_P2': {'subsetOf': set(), 'sameset': set(), 'isSubsumable': True},
            '05_P3': {'subsetOf': set(), 'sameset': set(), 'isSubsumable': False},
            '06_P1': {'subsetOf': set(), 'sameset': set(), 'isSubsumable': False},
            '06_P2': {'subsetOf': {'06_P1', '06_P4'}, 'sameset': set(), 'isSubsumable': False},
            '06_P3': {'subsetOf': {'06_P4'}, 'sameset': {'06_P3', '06_P5'}, 'isSubsumable': False},
            '06_P4': {'subsetOf': set(), 'sameset': set(), 'isSubsumable': False},
            '06_P5': {'subsetOf': {'06_P4'}, 'sameset': {'06_P3', '06_P5'}, 'isSubsumable': False},
            '07_P1': {'subsetOf': set(), 'sameset': set(), 'isSubsumable': False},
            '07_P2': {'subsetOf': set(), 'sameset': set(), 'isSubsumable': False},
            '07_P3': {'subsetOf': set(), 'sameset': set(), 'isSubsumable': True},
            '07_P4': {'subsetOf': set(), 'sameset': set(), 'isSubsumable': False},
            '08_P1': {'subsetOf': set(), 'sameset': set(), 'isSubsumable': False},
            '08_P2': {'subsetOf': set(), 'sameset': set(), 'isSubsumable': False},
            '08_P3': {'subsetOf': set(), 'sameset': set(), 'isSubsumable': True},
        }

        #Expected grouping results on the peptide level
        groupingPeptideResults = {
            'uniquePeptides': {
                '01_p1', '01_p2', '02_p1', '03_p1', '03_p2', '04_p1', '06_p1',
                '07_p4', '07_p5', '08_p1', '08_p4'
            },
            'groupUniquePeptides': {
                '02_p2', '04_p2', '06_p4'
            },
            'groupSubsumablePeptides': {
                '05_p1', '05_p4', '07_p2', '08_p2', '08_p3'
            },
            'sharedPeptides': {
                '05_p2', '05_p3', '06_p2', '06_p3', '07_p1', '07_p3'
            },
        }

        #Store all test data
        self.proteinToPeptides = proteinToPeptides
        self.mergedProteinToPeptides = mergedProteinToPeptides
        self.peptideToProteins = peptideToProteins
        self.mergedPeptideToProteins = mergedPeptideToProteins
        self.observedPeptides = observedPeptides
        self.potentialProteins = potentialProteins

        self.expectedClusters = expectedClusters
        self.equalProteins = equalProteins
        self.equalProteinTuples = sorted([tuple(p) for p in equalProteins])

        self.uniqueProteins = uniqueProteins
        self.uniqueMergedProteins = uniqueMergedProteins

        self.subsetProteins = sorted(subsetProteins)
        self.subsetMergedProteins = sorted(subsetMergedProteins)
        self.redundantProteins = redundantProteins

        self.groupingGroupResults = groupingGroupResults
        self.groupingProteinResults = groupingProteinResults
        self.groupingPeptideResults = groupingPeptideResults


class TestProteinGroupClass(unittest.TestCase):
    def TestAddProteins(self):
        proteinGroup = MODULE.ProteinGroup('groupId', 'leader1')
        proteinGroup.addLeadingProteins(['leader1', 'leader2'])
        proteinGroup.addLeadingProteins(['leader3'])
        proteinGroup.addSubsetProteins({'subset1'})
        proteinGroup.addSubsetProteins({'subset3', 'subset2'})
        proteinGroup.addSubsumableProteins('subsumable1')
        proteinGroup.addSubsumableProteins('subsumable2')
        proteinGroup.addSubsumableProteins('subsumable3')

        self.assertSetEqual(proteinGroup.leading, {'leader1', 'leader2', 'leader3'})
        self.assertSetEqual(proteinGroup.subset, {'subset1', 'subset2', 'subset3'})
        self.assertSetEqual(proteinGroup.subsumableProteins,
                            {'subsumable1', 'subsumable2', 'subsumable3'}
                            )
        self.assertSetEqual(proteinGroup.proteins, {'leader1', 'leader2', 'leader3',
                                                    'subset1', 'subset2', 'subset3',
                                                    }
                            )

    #class TestProteinInferenceClass(unittest.TestCase):
    def TestAddProteinGroup(self):
        inference = MODULE.ProteinInference(dict())
        groupId = inference.addProteinGroup('leader1')
        self.assertEqual(inference.groups[groupId].representative, 'leader1')
        self.assertEqual(inference._nextGroupId, 2)
        groupId = inference.addProteinGroup('leader1')
        self.assertEqual(inference._nextGroupId, 3)

    def TestAddProteinsToGroups(self):
        inference = MODULE.ProteinInference(dict())
        id_1 = inference.addProteinGroup('leader1')
        id_2 = inference.addProteinGroup('leader1')

        inference.addLeadingToGroups('leader3', id_1)
        inference.addSubsetToGroups(['subset1', 'subset2'], id_1)
        inference.addLeadingToGroups({'leader1', 'leader2'}, id_2)
        inference.addSubsumableToGroups('subsumable1', {id_1, id_2})

        group1 = inference.groups[id_1]
        self.assertSetEqual(group1.leading, {'leader1', 'leader3'})
        self.assertSetEqual(group1.subset, {'subset1', 'subset2'})
        self.assertSetEqual(group1.subsumableProteins, {'subsumable1'})

        group2 = inference.groups[id_2]
        self.assertSetEqual(group2.leading, {'leader1', 'leader2'})
        self.assertSetEqual(group2.subset, set())
        self.assertSetEqual(group2.subsumableProteins, {'subsumable1'})

        self.assertSetEqual({g.id for g in inference.getGroups('leader1')}, {id_1, id_2})
        self.assertSetEqual({g.id for g in inference.getGroups('leader2')}, {id_2})
        self.assertSetEqual({g.id for g in inference.getGroups('leader3')}, {id_1})
        self.assertSetEqual({g.id for g in inference.getGroups('subset1')}, {id_1})
        self.assertSetEqual({g.id for g in inference.getGroups('subset2')}, {id_1})
        self.assertSetEqual({g.id for g in inference.getGroups('subsumable1')}, {id_1, id_2})



class TestProteinInferenceFunctions(SetupMappingExamples):
    def test_findProteinClusters(self):
        proteinClusters = MODULE._findProteinClusters(self.proteinToPeptides,
                                                     self.peptideToProteins)
        self.assertEqual(len(proteinClusters), len(self.expectedClusters))
        allClusteredProteins = reduce(lambda c1, c2: c1.union(c2), proteinClusters)
        self.assertSetEqual(allClusteredProteins, self.potentialProteins)

        #Assert that all proteins withing a cluster belong the correct cluster,
        #which is indicated by the number infront of the first "_" of the protein names.
        for cluster in proteinClusters:
            self.assertIsInstance(cluster, set)
            self.assertEqual(len({p.split('_')[0] for p in cluster}), 1)

    def test_findSamesetProteins(self):
        equalProteins = MODULE._findSamesetProteins(self.proteinToPeptides,
                                                    proteins=self.proteinToPeptides)
        equalProteinTuples = sorted(equalProteins)
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

    def test_flattenMergedProteins(self):
        flatProteins = MODULE._flattenMergedProteins(self.mergedProteinToPeptides)
        self.assertSetEqual(flatProteins, self.potentialProteins)

    def test_mappingGetValueSet(self):
        valueSet1 = MODULE._mappingGetValueSet(self.proteinToPeptides, ['07_P1'])
        self.assertSetEqual(valueSet1, {'07_p1', '07_p2'})
        valueSet2 = MODULE._mappingGetValueSet(self.proteinToPeptides, ['07_P1', '07_P2'])
        self.assertSetEqual(valueSet2, {'07_p1', '07_p2', '07_p3', '07_p4'})
        valueSet3 = MODULE._mappingGetValueSet(self.proteinToPeptides, ['07_P2', '07_P4'])
        self.assertSetEqual(valueSet3, {'07_p1', '07_p3', '07_p4', '07_p5'})
        valueSet4 = MODULE._mappingGetValueSet(self.proteinToPeptides,
                                               ['07_P1', '07_P2', '07_P3', '07_P4'])
        self.assertSetEqual(valueSet4, {'07_p1', '07_p2', '07_p3', '07_p4', '07_p5'})

    def test_reducedProtToPeps(self):
        proteins = ['02_P1', '03_P1']
        reduceProteinToPeptides = MODULE._reducedProtToPeps(self.proteinToPeptides, proteins)
        for protein in proteins:
            self.assertNotIn(protein, reduceProteinToPeptides)
            self.assertIn(protein, self.proteinToPeptides)

    def test_findUniqueMappingValues(self):
        uniqueProteins = MODULE._findUniqueMappingValues(self.peptideToProteins)
        self.assertSetEqual(uniqueProteins, self.uniqueProteins)
        uniqueMergedProteins = MODULE._findUniqueMappingValues(self.mergedPeptideToProteins)
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
        uniqueProteins = MODULE._findUniqueMappingValues(pepToProts)
        self.assertSetEqual(uniqueProteins, set(protToPeps))

        for peptide, proteins in viewitems(pepToProts):
            self.assertGreater(len(proteins), 0)

    def test_getValueCounts(self):
        peptideFrequency = MODULE._getValueCounts(self.peptideToProteins)
        for peptide, counts in viewitems(peptideFrequency):
            self.assertEqual(counts, len(self.peptideToProteins[peptide]))

        protPepCounts = MODULE._getValueCounts(self.proteinToPeptides)
        for protein, counts in viewitems(protPepCounts):
            self.assertEqual(counts, len(self.proteinToPeptides[protein]))


class TestMappingBasedGrouping(SetupMappingExamples):
    def test_mappingBasedGrouping_proteinGroups(self):
        proteinInference = MODULE.mappingBasedGrouping(self.proteinToPeptides)

        # - Tests - #
        #number of protein clusters must be equal to the calculated number
        numClusters = len({p.split('_')[0] for p in self.proteinToPeptides})
        #self.assertEqual(, numClusters)

        allProteins = set()
        for group in viewvalues(proteinInference.groups):
            #Each protein group must contain exactly one "representative" protein
            self.assertIn(group.representative, group.leading)
            #Each protein group must at least contain one "leading" protein
            self.assertGreaterEqual(len(group.leading), 1)
            #Collect all group proteins and assosicated proteins (=subsumable)
            allProteins.update(group.proteins)
            allProteins.update(group.subsumableProteins)

        #The union of all proteinGroup.proteins must be equal to potentialProteins
        self.assertSetEqual(allProteins, self.potentialProteins)

        #Evaluate protein groups
        for reprProtein, groupingResults in viewitems(self.groupingGroupResults):
            groups = proteinInference.getGroups(reprProtein)
            self.assertEqual(len(groups), 1)
            group = groups[0]
            self.assertEqual(group.leading, groupingResults['leading'])
            self.assertSetEqual(group.subset, groupingResults['subset'])
            self.assertSetEqual(group.subsumableProteins, groupingResults['subsumable'])

    def test_mappingBasedGrouping_proteins(self):
        proteinInference = MODULE.mappingBasedGrouping(self.proteinToPeptides)
        #Compore protein entries to expected results
        for proteinId, expectedResults in viewitems(self.groupingProteinResults):
            proteinEntry = proteinInference.proteins[proteinId]
            self.assertSetEqual(proteinEntry.isSubset, expectedResults['subsetOf'])
            self.assertSetEqual(proteinEntry.isSameset, expectedResults['sameset'])
            self.assertEqual(proteinEntry.isSubsumable, expectedResults['isSubsumable'])
        #Assert agreement of protein entry relations and protein groups.
        for proteinId, proteinEntry in viewitems(proteinInference.proteins):
            if proteinEntry.isLeading:
                proteinGroup = proteinInference.getGroups(proteinId)[0]
                self.assertIn(proteinId, proteinGroup.leading)
            if proteinEntry.isSubset:
                proteinGroup = proteinInference.getGroups(proteinId)[0]
                self.assertIn(proteinId, proteinGroup.subset)
            if proteinEntry.isSubsumable:
                for proteinGroup in proteinInference.getGroups(proteinId):
                    self.assertIn(proteinId, proteinGroup.subsumableProteins)

    def test_mappingBasedGrouping_proteinPeptideMapping(self):
        proteinInference = MODULE.mappingBasedGrouping(self.proteinToPeptides)
        #Todo, check peptide assignment

        for peptide, proteinIds in viewitems(proteinInference.pepToProts):
            for proteinId in proteinIds:
                proteinEntry = proteinInference.proteins[proteinId]

                #Assert that all peptides have been characterized
                combinedProteinPeps = set()
                combinedProteinPeps.update(proteinEntry.uniquePeptides)
                combinedProteinPeps.update(proteinEntry.groupUniquePeptides)
                combinedProteinPeps.update(proteinEntry.groupSubsumablePeptides)
                combinedProteinPeps.update(proteinEntry.sharedPeptides)
                self.assertSetEqual(combinedProteinPeps, proteinEntry.peptides)

                #Assert that characterization of peptides is as expected
                for peptide in proteinEntry.uniquePeptides:
                    self.assertIn(peptide, self.groupingPeptideResults['uniquePeptides'])
                for peptide in proteinEntry.groupUniquePeptides:
                    self.assertIn(peptide, self.groupingPeptideResults['groupUniquePeptides'])
                for peptide in proteinEntry.groupSubsumablePeptides:
                    self.assertIn(peptide, self.groupingPeptideResults['groupSubsumablePeptides'])
                for peptide in proteinEntry.sharedPeptides:
                    self.assertIn(peptide, self.groupingPeptideResults['sharedPeptides'])























