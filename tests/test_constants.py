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
import unittest
import maspy.constants as module

class TestMasses(unittest.TestCase):
    def test_aaComp(self):
        self.assertEqual(module.aaComp['A'], module.COMPOSITION({'H': 5, 'C': 3, 'O': 1, 'N': 1}))

    def test_aaMass(self):
        self.assertEqual(round(module.aaMass['A'], 6), 71.037114)

    def test_aaModComp(self):
        self.assertEqual(module.aaModComp['u:1'], module.COMPOSITION({'C': 2, 'H': 2, 'O': 1}))

    def test_aaModMass(self):
        self.assertEqual(round(module.aaModMass['u:1'], 6), 42.010565)

