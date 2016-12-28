import sys
sys.path.append('D:/Dropbox/python/maspy')
import unittest
import maspy.ontology as module

class TestOntology(unittest.TestCase):
    def test_DefaultTranslator(self):
        translator = module.DefaultTranslator()
        testTermIdToName = {'MS:1000011': 'mass resolution',
                            'MS:1002289': 'ProteinProphet',
                            'MS:1002461': 'protein group-level confidence',
                            'MS:1002089': 'ProteomeDiscoverer:Peptide Without Protein XCorr Threshold',
                            'UO:0000010': 'second',
                            'UO:0000052': 'mass density unit',
                            'UO:0000031': 'minute',
                            'UO:0000111': 'energy unit'
                            }
        for termId in testTermIdToName:
            termName = testTermIdToName[termId]
            self.assertEqual(translator.getNameWithId(termId), termName)