"""
TODO
"""

#  Copyright 2015-2017 David M. Hollenstein, Jakob J. Hollenstein
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.

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

import maspy.inference as INFERENCE
import maspy._proteindb_refactoring as PROTEINDB

#Import your proteindb file
def the_magic_mapping_function(peptides, fastaPath, importAttributes=None):
    """ TODO

    :param peptides: a set of peptide sequences
    :parm fastaPath: FASTA file path
    :param importAttributes: dict, can be used to override default parameters
        passed to the function maspy.proteindb.importProteinDatabase().
    :returns: TODO
    """

    missedCleavage = max([p.count('K') + p.count('R') for p in peptides]) - 1
    minLength = min([len(p) for p in peptides])
    maxLength = max([len(p) for p in peptides])
    defaultAttributes = {
        'cleavageRule': '[KR]', 'minLength': minLength, 'maxLength': maxLength,
        'removeNtermM': True,  'ignoreIsoleucine': True,
        'missedCleavage': missedCleavage, 'forceId': True,
        'headerParser': PROTEINDB.fastaParserSpectraClusterPy,
    }
    if importAttributes is not None:
        defaultAttributes.update(importAttributes)

    proteindb = PROTEINDB.importProteinDatabase(fastaPath, **defaultAttributes)

    #This could be automated by adding a function to the inference module
    proteinToPeptides = ddict(set)
    for peptide in peptides:
        proteins = proteindb.peptides[peptide].proteins
        for protein in proteins:
            proteinToPeptides[protein].add(peptide)

    #Generate the ProteinInference instance
    inference = INFERENCE.mappingBasedGrouping(proteinToPeptides)

    peptideGroupMapping = dict()
    for peptide in peptides:
        groupLeaders = set()
        for proteinId in inference.pepToProts[peptide]:
            for proteinGroup in inference.getGroups(proteinId):
                groupLeaders.add(';'.join(sorted(proteinGroup.leading)))
        peptideGroupMapping[peptide] = groupLeaders

    return peptideGroupMapping
