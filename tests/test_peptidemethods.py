import sys
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
        self.assertEqual(round(module.calcPeptideMass('AADITSLYK'), 6), 980.517862)
        self.assertEqual(round(module.calcPeptideMass('AADITS[UNIMOD:21]LYK'), 6), 1060.484192)

    def test_removeModifications(self):
        self.assertEqual(module.removeModifications('[UNIMOD:1]AADI[DSS]T[UNIMOD:21]SLYK'), 'AADITSLYK')

    def test_returnModPositions(self):
        self.assertDictEqual(module.returnModPositions('[UNIMOD:1]AADI[DSS]T[UNIMOD:21]SLYK'), {'1': [1], '21': [5], 'DSS': [4]})
        self.assertDictEqual(module.returnModPositions('[UNIMOD:1]AADI[DSS]T[UNIMOD:21]SLYK', removeModString=False), {'DSS': [4], 'UNIMOD:1': [1], 'UNIMOD:21': [5]})
        self.assertDictEqual(module.returnModPositions('[UNIMOD:1]AADI[DSS]T[UNIMOD:21]SLYK', indexStart=51), {'1': [51], '21': [55], 'DSS': [54]})


class TestInsilicoDigestionMethods(unittest.TestCase):
    def test_digestInSilico(self):
        polypeptide = 'MSFLNASCTLCDEPISNRRKGEKIIELACGHLSHQECLIISGEIIELACGHLSHQECLIISECLIISGEIIELACGHLSHQECLIISECLIISGEIIELACGHLSHQECLIISRADAI'

        #cleavage rules
        params = {'cleavageRule':'[R]', 'missedCleavage':0, 'removeNtermM':False, 'minLength':1, 'maxLength':1000}
        self.assertTrue(len(module.digestInSilico(polypeptide, **params)) == polypeptide.count('R') + 1)

        params = {'cleavageRule':'[K]', 'missedCleavage':0, 'removeNtermM':False, 'minLength':1, 'maxLength':1000}
        self.assertTrue(len(module.digestInSilico(polypeptide, **params)) == polypeptide.count('K') + 1)

        params = {'cleavageRule':'[R]', 'missedCleavage':0, 'removeNtermM':False, 'minLength':1, 'maxLength':1000}
        peptides = [pep for pep, info in module.digestInSilico(polypeptide, **params)]
        self.assertIn('MSFLNASCTLCDEPISNR', peptides)
        self.assertNotIn('MSFLNASCTLCDEPISNRRK', peptides)

        params = {'cleavageRule':'[K]', 'missedCleavage':0, 'removeNtermM':False, 'minLength':1, 'maxLength':1000}
        peptides = [pep for pep, info in module.digestInSilico(polypeptide, **params)]
        self.assertIn('MSFLNASCTLCDEPISNRRK', peptides)
        self.assertNotIn('MSFLNASCTLCDEPISNR', peptides)
        self.assertNotIn('MSFLNASCTLCDEPISNRRKGEK', peptides)
        self.assertNotIn('SFLNASCTLCDEPISNRRK', peptides)

        params = {'cleavageRule':'\w(?=[R])', 'missedCleavage':0, 'removeNtermM':False, 'minLength':1, 'maxLength':1000}
        peptides = [pep for pep, info in module.digestInSilico(polypeptide, **params)]
        self.assertIn('MSFLNASCTLCDEPISN', peptides)
        self.assertNotIn('MSFLNASCTLCDEPISNR', peptides)

        params = {'cleavageRule':'[K]', 'missedCleavage':1, 'removeNtermM':False, 'minLength':23, 'maxLength':23}
        digestionresults = module.digestInSilico(polypeptide, **params)
        self.assertEqual(digestionresults, [('MSFLNASCTLCDEPISNRRKGEK', {'endPos': 23, 'missedCleavage': 1, 'startPos': 1})])

        #min and max length
        params = {'cleavageRule':'[K]', 'missedCleavage':0, 'removeNtermM':False, 'minLength':5, 'maxLength':1000}
        peptides = [pep for pep, info in module.digestInSilico(polypeptide, **params)]
        self.assertIn('MSFLNASCTLCDEPISNRRK', peptides)
        self.assertNotIn('GEK', peptides)

        params = {'cleavageRule':'[K]', 'missedCleavage':0, 'removeNtermM':False, 'minLength':1, 'maxLength':95}
        peptides = [pep for pep, info in module.digestInSilico(polypeptide, **params)]
        self.assertEqual(len(peptides), 3)

        params = {'cleavageRule':'[K]', 'missedCleavage':0, 'removeNtermM':False, 'minLength':1, 'maxLength':94}
        peptides = [pep for pep, info in module.digestInSilico(polypeptide, **params)]
        self.assertEqual(len(peptides), 2)

        #missed cleavage
        params = {'cleavageRule':'[K]', 'missedCleavage':1, 'removeNtermM':False, 'minLength':5, 'maxLength':1000}
        peptides = [pep for pep, info in module.digestInSilico(polypeptide, **params)]
        self.assertEqual(len(peptides), 4)
        self.assertIn('MSFLNASCTLCDEPISNRRKGEK', peptides)

        #remove nterminal methionine
        params = {'cleavageRule':'[K]', 'missedCleavage':1, 'removeNtermM':True, 'minLength':5, 'maxLength':40}
        peptides = [pep for pep, info in module.digestInSilico(polypeptide, **params)]
        self.assertEqual(len(peptides), 4)
        self.assertIn('SFLNASCTLCDEPISNRRK', peptides)
        self.assertIn('SFLNASCTLCDEPISNRRKGEK', peptides)
