import sys
sys.path.append('D:/Dropbox/python')
import unittest
import maspy.peptidemethods as module

class TestTransformMassToMzMethods(unittest.TestCase):
    def test_calcMhFromMz(self):
        self.assertEqual(module.calcMhFromMz(327.8465636835834, 3), 981.5251381172102)
        self.assertEqual(module.calcMhFromMz(491.26620729199004, 2), 981.5251381172101)

    def test_calcMzFromMh(self):
        self.assertEqual(module.calcMzFromMh(981.5251381172101, 3), 327.8465636835834)
        self.assertEqual(module.calcMzFromMh(981.5251381172101, 2), 491.26620729199004)

    def test_calcMzFromMass(self):
        self.assertEqual(module.calcMzFromMass(980.5178616504401, 3), 327.8465636835834)
        self.assertEqual(module.calcMzFromMass(980.5178616504401, 2), 491.26620729199004)

    def test_calcMassFromMz(self):
        self.assertEqual(module.calcMassFromMz(491.26620729199004, 2), 980.5178616504401)
        self.assertEqual(module.calcMassFromMz(327.8465636835834, 3), 980.5178616504402)


class TestPeptideSequenceMethods(unittest.TestCase):
    def test_calcPeptideMass(self):
        self.assertEqual(module.calcPeptideMass('AADITSLYK'), 980.5178616504401)
        self.assertEqual(module.calcPeptideMass('AADITS[UNIMOD:21]LYK'), 1060.48419265044)

    def test_removeModifications(self):
        self.assertEqual(module.removeModifications('[UNIMOD:1]AADI[DSS]T[UNIMOD:21]SLYK'), 'AADITSLYK')

    def test_returnModPositions(self):
        self.assertDictEqual(module.returnModPositions('[UNIMOD:1]AADI[DSS]T[UNIMOD:21]SLYK'), {'1': [1], '21': [5], 'DSS': [4]})
        self.assertDictEqual(module.returnModPositions('[UNIMOD:1]AADI[DSS]T[UNIMOD:21]SLYK', removeModString=False), {'DSS': [4], 'UNIMOD:1': [1], 'UNIMOD:21': [5]})
        self.assertDictEqual(module.returnModPositions('[UNIMOD:1]AADI[DSS]T[UNIMOD:21]SLYK', indexStart=51), {'1': [51], '21': [55], 'DSS': [54]})


class TestInsilicoDigestionMethods(unittest.TestCase):
    def test_digestInSilico(self):
        pass
        #module.digestInSilico(proteinSequence, cleavageRule='[KR]', missedCleavages=0, removeNtermM=True, minLength=5, maxLength=40)
