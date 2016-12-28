"""
This module provides access to OBO ontology files. The class DefaulTranslator
already contains OBO Terms used in mzML files and other file formats specified
by the Proteomics Standards Initiative (PSI).
"""
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
#Refer to http://pythonhosted.org/Orange-Bioinformatics/index.html and the 
#OBO ontology module for a much more comprehensive OBO library.
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
                oldTermIsObsolete = termIsObsolete(oldOboTerm)
                newTermIsObsolete = termIsObsolete(oboTerm)
                if oldTermIsObsolete and not newTermIsObsolete:
                    self.oboTerms[oboTerm['id']] = oboTerm
                else:
                    #At least one of two terms with identical id must be obsolete
                    assert oldTermIsObsolete or newTermIsObsolete

    def __getitem__(self, oboTermId):
        return self.oboTerms[oboTermId]

    def getNameFromId(self, oboTermId):
        """Return the term name by using the term id."""
        return self.oboTerms[oboTermId]['name']

    def getIdFromName(self, oboTermName):
        """Return the term id by using the term name."""
        raise NotImplementedError


class DefaultTranslator(OboTranslator):
    """A class that provides access to OBO ontology terms used in mzML files
    and other file formats specified by the Proteomics Standards Initiative.
    Imports the files 'psi-ms.obo' and 'unit.obo' which are distributed together
    with maspy.

    Use :func:`maspy.ontology.DefaultTranslator.getNameFromId()` to get a terms
    'name' by its 'id'.
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
