"""
The module "inference" is still in development. The interface of high and low
level functions is not yet stable!

This module will allow a simple protein inference; reporting a minimal set of
proteins that can explain a set of peptides; defining shared und unique
peptides; defining proteins that share equal evidence and protein subgroups;
defining protein redundancy groups.
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

import warnings
warnings.warn('Module maspy.inference is currently in development.',
              ImportWarning)

from collections import Counter
from collections import defaultdict as ddict
import operator

"""
Mapping objects
---------------

:param pepToProts: dict, for each peptide (=key) contains a set of parent
    proteins (=value). For Example {peptide: {protein, ...}, ...}
:param protToPeps: dict, for each protein (=key) contains a set of
    associated peptides (=value). For Example {protein: {peptide, ...}, ...}
"""

def _groupConnectedPeptides(pepToProts, protToPeps):
    """ #TODO

    :param pepToProts: dict, for each peptide (=key) contains a set of parent
        proteins (=value). For Example {peptide: {protein, ...}, ...}
    :param protToPeps: dict, for each protein (=key) contains a set of
        associated peptides (=value). For Example {protein: {peptide, ...}, ...}
    :returns: a list of peptide sets
    """
    connectedGroups = list()
    resolvingPeptides = set(pepToProts)
    while resolvingPeptides:
        peptide = resolvingPeptides.pop()
        peptideGroup = set([peptide])
        proteins = set(pepToProts[peptide])
        parsedProteins = set()

        while len(proteins) != len(parsedProteins):
            for protein in proteins:
                peptideGroup.update(protToPeps[protein])
            parsedProteins.update(proteins)

            for peptide in peptideGroup:
                proteins.update(pepToProts[peptide])
        connectedGroups.append(peptideGroup)
        resolvingPeptides = resolvingPeptides.difference(peptideGroup)
    return connectedGroups


def _findUniqueProteins(pepToProts):
    """ #TODO

    :param pepToProts: dict, for each peptide (=key) contains a set of parent
        proteins (=value). For Example {peptide: {protein, ...}, ...}
    :returns: #TODO
    """
    uniqueProteins = list()
    for proteins in viewvalues(pepToProts):
        if len(proteins) == 1:
            uniqueProteins.extend(list(proteins))
    return set(uniqueProteins)


def _findEqualEvidenceProteins(proteins, protToPeps):
    """#TODO

    :param proteins: iterable, proteins that are tested for having equal
        evidence.
    :param protToPeps: dict, for each protein (=key) contains a set of
        associated peptides (=value). For Example {protein: {peptide, ...}, ...}
    :returns: a list of protein sets that share equal peptide evidence
    """
    equalEvidence = ddict(set)
    for protein in proteins:
        peptides = protToPeps[protein]
        equalEvidence[tuple(sorted(peptides))].add(protein)
    equalProteins = list()
    for proteins in viewvalues(equalEvidence):
        if len(proteins) > 1:
            equalProteins.append(proteins)
    return equalProteins


def _findSubsetProteins(proteins, protToPeps, pepToProts):
    """Find proteins which peptides are a subset of or equal to other proteins.

    :param proteins: iterable, proteins that are tested for being a subset
    :param pepToProts: dict, for each peptide (=key) contains a set of parent
        proteins (=value). For Example {peptide: {protein, ...}, ...}
    :param protToPeps: dict, for each protein (=key) contains a set of
        associated peptides (=value). For Example {protein: {peptide, ...}, ...}
    :returns: a list of pairs of protein and their superset proteins.
        [(protein, {superset protein, ...}), ...]
    """
    subGroups = list()
    for protein in proteins:
        peptideCounts = Counter()
        for peptide in protToPeps[protein]:
            proteins = pepToProts[peptide]
            peptideCounts.update(proteins)
        peptideCount = peptideCounts.pop(protein)

        superGroups = set()
        for sharingProtein, sharedPeptides in peptideCounts.most_common():
            if sharedPeptides == peptideCount:
                superGroups.add(sharingProtein)
            else:
                break
        if superGroups:
            subGroups.append((protein, superGroups))
    return subGroups


def _findRedundantProteins(protToPeps, pepToProts, proteins=None):
    """Returns a set of proteins with redundant peptide evidence.

    After removing the redundant proteins from the "protToPeps" and "pepToProts"
    mapping, all remaining proteins have at least one unique peptide. The
    remaining proteins are a "minimal" set of proteins that are able to explain
    all peptides. However, this is not guaranteed to be the optimal solution
    with the least number of proteins. In addition it is possible that multiple
    solutions with the same number of "minimal" proteins exist.

    Procedure for finding the redundant proteins:
    1.  Generate a list of proteins that do not contain any unique peptides, a
        unique peptide has exactly one protein entry in "pepToProts".
    2.  Proteins are first sorted in ascending order of the number of peptides.
        Proteins with an equal number of peptides are sorted in descending order
            of their sorted peptide frequencies (= proteins per peptide).
        If two proteins are still equal, they are sorted alpha numerical in
        ascending order according to their protein names.
    3.  Parse this list of sorted non unique proteins;
        If all its peptides have a frequency value of greater 1;
        mark the protein as redundant; remove its peptides from the peptide
        frequency count, continue with the next entry.
    4.  Return the set of proteins marked as redundant.

    :param pepToProts: dict, for each peptide (=key) contains a set of parent
        proteins (=value). For Example {peptide: {protein, ...}, ...}
    :param protToPeps: dict, for each protein (=key) contains a set of
        associated peptides (=value). For Example {protein: {peptide, ...}, ...}
    :param proteins: iterable, proteins that are tested for being redundant. If
        None all proteins in "protToPeps" are parsed.
    :returns: a set of redundant proteins, i.e. proteins that are not necessary
        to explain all peptides
    """
    if proteins is None:
        proteins = viewkeys(protToPeps)

    pepFrequency = _getValueCounts(pepToProts)
    protPepCounts = _getValueCounts(protToPeps)

    getCount = operator.itemgetter(1)
    getProt = operator.itemgetter(0)

    sort = list()
    for protein in sorted(proteins):
        protPepFreq = [pepFrequency[pep] for pep in protToPeps[protein]]
        if min(protPepFreq) > 1:
            sortValue = (len(protPepFreq)*-1, sorted(protPepFreq, reverse=True))
            sort.append((protein, sortValue))
    sortedProteins = map(getProt, sorted(sort, key=getCount, reverse=True))

    redundantProteins = set()
    for protein in sortedProteins:
        for pep in protToPeps[protein]:
            if pepFrequency[pep] <= 1:
                break
        else:
            protPepFrequency = Counter(protToPeps[protein])
            pepFrequency.subtract(protPepFrequency)
            redundantProteins.add(protein)
    return redundantProteins


def _mergeProteinEntries(proteinGroups, protToPeps):
    """Returns a new "protToPeps" dictionary with entries merged that are
    present in proteinGroups.

    NOTE:
        The key of the merged entry is a tuple of the sorted protein keys. This
        behaviour might change in the future; the tuple might be replaced by
        simply one of the protein entries which is then representative for all.

    :param proteinGroups: a list of protein groups that will be merged
        [{protein, ...}, ...]
    :param protToPeps: dict, for each protein (=key) contains a set of
        associated peptides (=value). For Example {protein: {peptide, ...}, ...}
    :returns: dict, {protein: set([peptid, ...])}
    """
    mergedProtToPeps = dict(protToPeps)
    for proteinGroup in proteinGroups:
        for protein in proteinGroup:
            peptides = mergedProtToPeps.pop(protein)
        mergedProtein = tuple(sorted(proteinGroup))
        mergedProtToPeps[mergedProtein] = peptides
    return mergedProtToPeps


def _reducedProtToPeps(protToPeps, proteins):
    """Returns a new, reduced "protToPeps" dictionary that does not contain
    entries present in "proteins".

    :param protToPeps: dict, for each protein (=key) contains a set of
        associated peptides (=value). For Example {protein: {peptide, ...}, ...}
    :param proteins: a list of proteinSet
    :returns: dict, protToPeps not containing entries from "proteins"
    """
    return {k: v for k, v in viewitems(protToPeps) if k not in proteins}


def _invertMapping(mapping):
    """Converts a protein to peptide or peptide to protein mapping.

    :param protToPeps: dict, for each protein (=key) contains a set of
        associated peptides (=value). For Example {protein: {peptide, ...}, ...}

    :returns: "pepToProts" ddict, for each peptide (=key) contains a set of
        parent proteins (=value). For Example {peptide: {protein, ...}, ...}
    """
    invertedMapping = ddict(set)
    for key, values in viewitems(mapping):
        for value in values:
            invertedMapping[value].add(key)
    return invertedMapping


def _getValueCounts(mapping):
    """ #TODO:

    :param mapping: dict, for each key contains a set of values. Can be either a
        protToPeps or pepToProts dictionary.
    :returns: #TODO:
    """
    return Counter({k: len(v) for k, v in viewitems(mapping)})

