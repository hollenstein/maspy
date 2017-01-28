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
