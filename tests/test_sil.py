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

import unittest

sys.path.append(os.path.abspath('..'))
import maspy.sil as MODULE


class TestIsotopicLabelMethods(unittest.TestCase):
    def setUp(self):
        self.labelDescriptor = MODULE.LabelDescriptor()
        self.labelDescriptor.addLabel({'K':'', 'R':''})
        self.labelDescriptor.addLabel({'K':'u:188', 'R':'u:188'})
        self.labelDescriptor.addLabel({'nTerm':'u:36', 'K':['u:188','u:36'], 'R':['u:188','u:36']}, {'u:1':'u:36'})
        self.labelDescriptor.addLabel({'nTerm':'u:199', 'K':['u:188','u:199'], 'R':['u:188','u:199']}, {'u:1':'u:199'})

    def test_complex_returnLabelState(self):
        peptide0 = 'KASFDJK'
        peptide1 = 'K[u:188]ASFDJK[u:188]'
        peptide2 = '[u:36]K[u:188][u:36]ASFDJK[u:188][u:36]'
        peptide3 = '[u:199]K[u:188][u:199]ASFDJK[u:188][u:199]'
        peptide4 = '[u:1]ASFDJ'
        peptide5 = '[u:1]K[u:188]ASFDJK[u:188]'
        peptide6 = 'KASFDJ[u:188]'
        peptide7 = '[u:1]KASFDJK[u:188]'
        peptide8 = '[u:1]ASF[u:188]DJ'

        self.assertEqual(MODULE.returnLabelState(peptide0, self.labelDescriptor), 0)
        self.assertEqual(MODULE.returnLabelState(peptide1, self.labelDescriptor), 1)
        self.assertEqual(MODULE.returnLabelState(peptide2, self.labelDescriptor), 2)
        self.assertEqual(MODULE.returnLabelState(peptide3, self.labelDescriptor), 3)
        self.assertEqual(MODULE.returnLabelState(peptide4, self.labelDescriptor), -3)
        self.assertEqual(MODULE.returnLabelState(peptide5, self.labelDescriptor), 1)
        self.assertEqual(MODULE.returnLabelState(peptide6, self.labelDescriptor), -2)
        self.assertEqual(MODULE.returnLabelState(peptide7, self.labelDescriptor), -2)
        self.assertEqual(MODULE.returnLabelState(peptide8, self.labelDescriptor), -2)

    def test_simple_returnLabelState(self):
        labelDescriptor = MODULE.LabelDescriptor()
        labelDescriptor.addLabel({'K':'', 'R':''})
        labelDescriptor.addLabel({'K':'u:188', 'R':'u:188'})

        peptide0 = 'ASFDJK'
        peptide1 = 'ASFDJK[u:188]'
        peptide2 = 'ASFDJ'

        self.assertEqual(MODULE.returnLabelState(peptide0, labelDescriptor), 0)
        self.assertEqual(MODULE.returnLabelState(peptide1, labelDescriptor), 1)
        self.assertEqual(MODULE.returnLabelState(peptide2, labelDescriptor), -1)

    def test_returnLabelStateMassDifferences(self):
        peptide0 = 'KASFDJK'
        peptide1 = '[u:1]KASFDJK'
        peptide2 = 'KASFDJ[u:188]'

        result = dict([(acc, round(mass, 6)) for acc, mass in
                       MODULE.returnLabelStateMassDifferences(peptide0, self.labelDescriptor).items()])
        self.assertEqual(result, {1: 12.040258, 2: 96.134158, 3: 108.209479})

        result = dict([(acc, round(mass, 6)) for acc, mass in
                       MODULE.returnLabelStateMassDifferences(peptide1, self.labelDescriptor).items()])
        self.assertEqual(result, {1: 12.040258, 2: 68.102858, 3: 76.153072})

        self.assertEqual(MODULE.returnLabelStateMassDifferences(peptide2, self.labelDescriptor), {})

    def test_modSymbolsFromLabelInfo(self):
        self.assertEqual(MODULE.modSymbolsFromLabelInfo(self.labelDescriptor), {'u:188', 'u:199', 'u:36'})

    def test_modAminoacidsFromLabelInfo(self):
        self.assertEqual(MODULE.modAminoacidsFromLabelInfo(self.labelDescriptor), {'K', 'R', 'nTerm'})

    def test_expectedLabelPosition(self):
        peptide0 = 'KASFDJK'
        peptide1 = '[u:1]KASFDJK'

        self.assertEqual(MODULE.expectedLabelPosition(peptide0, self.labelDescriptor.labels[0]), {})
        self.assertEqual(MODULE.expectedLabelPosition(peptide0, self.labelDescriptor.labels[1]),
                         {0: ['u:188'], 6: ['u:188']})
        self.assertEqual(MODULE.expectedLabelPosition(peptide0, self.labelDescriptor.labels[2]),
                         {0: ['u:188', 'u:36', 'u:36'], 6: ['u:188', 'u:36']})
        self.assertEqual(MODULE.expectedLabelPosition(peptide0, self.labelDescriptor.labels[3]),
                         {0: ['u:188', 'u:199', 'u:199'], 6: ['u:188', 'u:199']})

        self.assertEqual(MODULE.expectedLabelPosition(peptide1, self.labelDescriptor.labels[0]), {})
        self.assertEqual(MODULE.expectedLabelPosition(peptide1, self.labelDescriptor.labels[1]),
                         {0: ['u:188'], 6: ['u:188']})
        self.assertEqual(MODULE.expectedLabelPosition(peptide1, self.labelDescriptor.labels[2]),
                         {0: ['u:188', 'u:36'], 6: ['u:188', 'u:36']})
        self.assertEqual(MODULE.expectedLabelPosition(peptide1, self.labelDescriptor.labels[3]),
                         {0: ['u:188', 'u:199'], 6: ['u:188', 'u:199']})
