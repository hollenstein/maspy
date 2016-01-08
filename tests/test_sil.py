import sys
sys.path.append('D:/Dropbox/python')
import unittest
import pyms.sil as module


class TestIsotopicLabelMethods(unittest.TestCase):
    def setUp(self):
        self.labelDescriptor = module.LabelDescriptor()
        self.labelDescriptor.addLabel({'K':'', 'R':''})
        self.labelDescriptor.addLabel({'K':'188', 'R':'188'})
        self.labelDescriptor.addLabel({'nTerm':'36', 'K':['188','36'], 'R':['188','36']}, {'1':'36'})
        self.labelDescriptor.addLabel({'nTerm':'199', 'K':['188','199'], 'R':['188','199']}, {'1':'199'})

    def test_complex_returnLabelState(self):
        peptide0 = 'KASFDJK'
        peptide1 = 'K[UNIMOD:188]ASFDJK[UNIMOD:188]'
        peptide2 = '[UNIMOD:36]K[UNIMOD:188][UNIMOD:36]ASFDJK[UNIMOD:188][UNIMOD:36]'
        peptide3 = '[UNIMOD:199]K[UNIMOD:188][UNIMOD:199]ASFDJK[UNIMOD:188][UNIMOD:199]'
        peptide4 = '[UNIMOD:1]ASFDJ'
        peptide5 = '[UNIMOD:1]K[UNIMOD:188]ASFDJK[UNIMOD:188]'
        peptide6 = 'KASFDJ[UNIMOD:188]'
        peptide7 = '[UNIMOD:1]KASFDJK[UNIMOD:188]'
        peptide8 = '[UNIMOD:1]ASF[UNIMOD:188]DJ'

        self.assertEqual(module.returnLabelState(peptide0, self.labelDescriptor), 0)
        self.assertEqual(module.returnLabelState(peptide1, self.labelDescriptor), 1)
        self.assertEqual(module.returnLabelState(peptide2, self.labelDescriptor), 2)
        self.assertEqual(module.returnLabelState(peptide3, self.labelDescriptor), 3)
        self.assertEqual(module.returnLabelState(peptide4, self.labelDescriptor), -3)
        self.assertEqual(module.returnLabelState(peptide5, self.labelDescriptor), 1)
        self.assertEqual(module.returnLabelState(peptide6, self.labelDescriptor), -2)
        self.assertEqual(module.returnLabelState(peptide7, self.labelDescriptor), -2)
        self.assertEqual(module.returnLabelState(peptide8, self.labelDescriptor), -2)

    def test_simple_returnLabelState(self):
        labelDescriptor = module.LabelDescriptor()
        labelDescriptor.addLabel({'K':'', 'R':''})
        labelDescriptor.addLabel({'K':'188', 'R':'188'})

        peptide0 = 'ASFDJK'
        peptide1 = 'ASFDJK[UNIMOD:188]'
        peptide2 = 'ASFDJ'

        self.assertEqual(module.returnLabelState(peptide0, labelDescriptor), 0)
        self.assertEqual(module.returnLabelState(peptide1, labelDescriptor), 1)
        self.assertEqual(module.returnLabelState(peptide2, labelDescriptor), -1)

    def test_returnLabelStateMassDifferences(self):
        peptide0 = 'KASFDJK'
        peptide1 = '[UNIMOD:1]KASFDJK'
        peptide2 = 'KASFDJ[UNIMOD:188]'

        self.assertEqual(module.returnLabelStateMassDifferences(peptide0, self.labelDescriptor), {1: 12.040258, 2: 96.134158, 3: 108.20947899999999})
        self.assertEqual(module.returnLabelStateMassDifferences(peptide1, self.labelDescriptor), {1: 12.040258, 2: 68.102858, 3: 76.153072})
        self.assertEqual(module.returnLabelStateMassDifferences(peptide2, self.labelDescriptor), {})

    def test_modSymbolsFromLabelInfo(self):
        self.assertEqual(module.modSymbolsFromLabelInfo(self.labelDescriptor), {'188', '199', '36'})

    def test_modAminoacidsFromLabelInfo(self):
        self.assertEqual(module.modAminoacidsFromLabelInfo(self.labelDescriptor), {'K', 'R', 'nTerm'})

    def test_expectedLabelPosition(self):
        peptide0 = 'KASFDJK'
        peptide1 = '[UNIMOD:1]KASFDJK'

        self.assertEqual(module.expectedLabelPosition(peptide0, self.labelDescriptor.labels[0]), {})
        self.assertEqual(module.expectedLabelPosition(peptide0, self.labelDescriptor.labels[1]), {0: ['188'], 6: ['188']})
        self.assertEqual(module.expectedLabelPosition(peptide0, self.labelDescriptor.labels[2]), {0: ['188', '36', '36'], 6: ['188', '36']})
        self.assertEqual(module.expectedLabelPosition(peptide0, self.labelDescriptor.labels[3]), {0: ['188', '199', '199'], 6: ['188', '199']})

        self.assertEqual(module.expectedLabelPosition(peptide1, self.labelDescriptor.labels[0]), {})
        self.assertEqual(module.expectedLabelPosition(peptide1, self.labelDescriptor.labels[1]), {0: ['188'], 6: ['188']})
        self.assertEqual(module.expectedLabelPosition(peptide1, self.labelDescriptor.labels[2]), {0: ['188', '36'], 6: ['188', '36']})
        self.assertEqual(module.expectedLabelPosition(peptide1, self.labelDescriptor.labels[3]), {0: ['188', '199'], 6: ['188', '199']})
