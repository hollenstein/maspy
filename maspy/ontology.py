"""
This module provides access to OBO ontology files. The class DefaulTranslator
already contains OBO Terms used in mzML files and other file formats specified
by the Proteomics Standards Initiative (PSI).

Refer to http://pythonhosted.org/Orange-Bioinformatics/index.html and the OBO
ontology module for a much more comprehensive OBO library.
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

import io
import os

import maspy.auxiliary as aux
import maspy.errors


def oboTermParser(filepath):
    """Read a obo file and yield '[Term]' entries.

    :param filepath: file path of the .obo file

    :yields: lists containing all lines from a obo '[Term]' entry. Lines are not
        processed and still contain the newline character.
    """
    with io.open(filepath) as openfile:
        lineIter = iter([i.rstrip() for i in openfile.readlines()])

    #Iterate through lines until the first obo "[Term]" is encountered
    try:
        line = lineIter.next()
        while line != '[Term]':
            line = lineIter.next()
        header = line #Remove
        entryLines = list()
    except StopIteration:
        errorText = 'File does not contain obo "[Term]" entries.'
        raise maspy.errors.FileFormatError(errorText)

    for line in lineIter:
        #Skip empty lines between entries
        if not line:
            continue
        if line == '[Term]':
            yield entryLines
            header = line #Remove
            entryLines = list()
        else:
            entryLines.append(line)

    #Yield last entry
    if entryLines:
        yield entryLines


def _attributeLinesToDict(attributeLines):
    """Converts a list of obo 'Term' lines to a dictionary.

    :param attributeLines: a list of obo 'Term' lines. Each line contains a key
        and a value part which are separated by a ':'.

    :return: a dictionary containing the attributes of an obo 'Term' entry.

    NOTE: Some attributes can occur multiple times in one single term, for 
          example 'is_a' or 'relationship'. However, currently only the last
          occurence is stored.
    """
    attributes = dict()
    for line in attributeLines:
        attributeId, attributeValue = line.split(':', 1)
        attributes[attributeId.strip()] = attributeValue.strip()
    return attributes


def _termIsObsolete(oboTerm):
    """Determine wheter an obo 'Term' entry is marked as obsolete.

    :param oboTerm: a dictionary as return by
        :func:`maspy.ontology._attributeLinesToDict()`

    :return: bool
    """
    isObsolete = False
    if u'is_obsolete' in oboTerm:
        if oboTerm[u'is_obsolete'].lower() == u'true':
            isObsolete = True
    return isObsolete


class OboTranslator(object):
    """A class that provides access to OBO ontology terms.

    Use OboTranslator.load(filepath) to import '[Term]' entries from an ontology
    file. Loaded terms are stored in self.oboTerms by using their 'id'.  

    :ivar oboTerms: a dictionary that stores imported obo terms in the form
        {termId: {attributeName: attributeValue, ...}, ...}

    NOTE: Some attributes can occur multiple times in one single term, for 
          example 'is_a' or 'relationship'. However, currently only the last
          occurence is stored.
    """
    def __init__(self):
        self.oboTerms = dict()

    def load(self, filepath):
        """Import '[Term]' entries from an .obo file."""
        for attributeLines in oboTermParser(filepath):
            oboTerm = _attributeLinesToDict(attributeLines)
            if oboTerm['id'] not in self.oboTerms:
                self.oboTerms[oboTerm['id']] = oboTerm            
            else:
                oldOboTerm = self.oboTerms[oboTerm['id']]
                oldTermIsObsolete = _termIsObsolete(oldOboTerm)
                newTermIsObsolete = _termIsObsolete(oboTerm)
                if oldTermIsObsolete and not newTermIsObsolete:
                    self.oboTerms[oboTerm['id']] = oboTerm
                else:
                    #At least one of two terms with identical id must be obsolete
                    assert oldTermIsObsolete or newTermIsObsolete

    def __getitem__(self, oboTermId):
        return self.oboTerms[oboTermId]

    def getNameWithId(self, oboTermId):
        """Return the term name by using the term id."""
        return self.oboTerms[oboTermId]['name']

    def getIdWithName(self, oboTermName):
        """Return the term id by using the term name."""
        raise NotImplementedError


class DefaultTranslator(OboTranslator):
    """A class that provides access to OBO ontology terms used in mzML files
    and other file formats specified by the Proteomics Standards Initiative.
    Imports the files 'psi-ms.obo' and 'unit.obo' which are distributed together
    with maspy.

    Use :func:`maspy.ontology.DefaultTranslator.getNameWithId()` to get a terms
    'name' by its 'id'.

    :ivar oboTerms: a dictionary that stores imported obo terms in the form
        {termId: {attributeName: attributeValue, ...}, ...}

    NOTE: Some attributes can occur multiple times in one single term, for 
          example 'is_a' or 'relationship'. However, currently only the last
          occurence is stored.
    """
    def __init__(self):
        super(DefaultTranslator, self).__init__()
        psiOntologyPath = aux.joinpath(os.path.dirname(aux.__file__),
                                       'ontologies', 'psi-ms.obo'
                                       )
        unitOntologyPath = aux.joinpath(os.path.dirname(aux.__file__),
                                        'ontologies', 'unit.obo'
                                        )
        self.load(psiOntologyPath)
        self.load(unitOntologyPath)
