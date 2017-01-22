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

