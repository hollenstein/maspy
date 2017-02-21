"""
The module "inference" is still in development. The interface of high and low
level functions is not yet stable!

This module will allow a simple protein inference; reporting a minimal set of
proteins that can explain a set of peptides; defining shared und unique
peptides; defining proteins that share equal evidence and protein subgroups;
defining protein redundancy groups.

Protein group
    A protein group is a collection of proteins that is necessary

.. Note:
    Whenever possible, we use the PSI-MS ontology terms to describe protein
    groups and protein relationships. However, in order to be less ambiguous
    some of the terms are more strictly defined in maspy or used in a slightly
    different context. The details of used terms are described below.

    Protein cluster (proteinCluster)
        A protein cluster comprises all proteins that are somehow linked by
        shared peptide evidence. (related to PSI-MS MS:1002407)

    Protein group (proteinGroup)
        A group of proteins that are essential to explain at least one peptide.
        A protein group consists of its leading, subset and subsumable proteins.
        However, subsumable proteins are more loosely associated with the group.
        The common peptide evidence of a group is the sum of all its leading and
        subset proteins. As a consequence, all leading and subset proteins must
        be a sameset of subset of the protein group common peptide evidence.

    Representative protein (representative)
        Each protein group must have exactle one representative protein. This
        protein must be amongst the leading proteins of the group. (identical to
        PSI-MS MS:1002403)

    Leading protein (leading)
        Each protein group must have at least one leading protein, to indicate
        proteins with the strongest evidence. (identical to PSI-MS MS:1002401)
        All proteins present in a protein group that are not leading proteins
        can be considered to be non-leading proteins (PSI-MS MS:1002402).

    Same-set protein (sameset)
        A protein that shares exactly the same set of peptide evidence with one
        or multiple other proteins. (identical to PSI-MS MS:1001594)

    Sub-set protein (subset)
        A protein with a sub-set of the peptide evidence of another protein. It
        can be a sub-set of multiple proteins, and those don't have to be in the
        same protein group. (identical to PSI-MS MS:1001596)

    Subsumable protein (subsumable)
        A protein which peptide evidence is distributed across multiple proteins
        in at least two different protein groups. This term is mutually exclusiv
        with the term sub-set protein, and hierarchically below sub-set. That
        means that a protein can only be a subsumable if it is not already a
        subset of another protein. Also the protein is not allowed to be a
        "leading" protein in any protein group. (is related to the PSI-MS term
        MS:1002570 "multiply subsumable protein", but defined more strictly)


Relevant terms from the PSI-MS ontology
---------------------------------------

Group representative (MS:1002403)
    An arbitrary and optional flag applied to exactly one protein per group to
    indicate it can serve as the representative of the group, amongst leading
    proteins, in effect serving as a tiebreaker for approaches that require
    exactly one group representative.

Leading protein (MS:1002401)
    At least one protein within each group should be annotated as a leading
    protein to indicate it has the strongest evidence, or approximately equal
    evidence as other group members.

Non-leading protein (MS:1002402)
    Zero to many proteins within each group should be annotated as non-leading
    to indicate that other proteins have stronger evidence.

Cluster identifier (MS:1002407)
    An identifier applied to protein groups to indicate that they are linked by
    shared peptides.

Same-set protein (MS:1001594)
    A protein which is indistinguishable or equivalent to another protein,
    having matches to an identical set of peptide sequences.

Sequence sub-set protein (MS:1001596)
    A protein with a sub-set of the peptide sequence matches for another
    protein, and no distinguishing peptide matches.

Sequence subsumable protein (MS:1001598)
    A sequence same-set or sequence sub-set protein where the matches are
    distributed across two or more proteins.

Sequence multiply subsumable protein (MS:1002570)
    A protein for which the matched peptide sequences are the same, or a
    subset of, the matched peptide sequences for two or more other proteins
    combined. These other proteins need not all be in the same group.

Peptide unique to one protein (MS:1001363)
    A peptide matching only one.

Peptide shared in multiple proteins (MS:1001175)
    A peptide matching multiple proteins.


Additional terms used in maspy
------------------------------

Unique protein (uniqueProtein)
    A protein mapped to at least one unique peptide.

Super-set protein:
    A protein with a set of peptide matches, that is a strict superset of other
    proteins (strict means that the peptide matches are not equal).
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

import maspy.auxiliary as AUX


"""
Mapping objects
---------------

:param pepToProts: dict, for each peptide (=key) contains a set of parent
    proteins (=value). For Example {peptide: {protein, ...}, ...}
:param protToPeps: dict, for each protein (=key) contains a set of
    associated peptides (=value). For Example {protein: {peptide, ...}, ...}
"""

class ProteinInference(object):
    """Contains the result of a protein inference analysis.

    :ivar protToPeps: dict, for each protein (=key) contains a set of
        associated peptides (=value). For Example {protein: {peptide, ...}, ...}
    :ivar pepToProts: dict, for each peptide (=key) contains a set of parent
        proteins (=value). For Example {peptide: {protein, ...}, ...}
    :ivar proteins: dict, protein ids pointing to ProteinObjects
    :ivar groups: dict, protein group ids pointing to ProteinGroup objects
    :ivar clusters: dict, cluster ids pointing to a list of protein group ids
    """
    def __init__(self, protToPeps):
        self.protToPeps = protToPeps
        self.pepToProts = _invertMapping(protToPeps)
        self.proteins = dict()
        self.groups = dict()
        self.clusters = dict()
        self._proteinToGroupIds = ddict(set)
        self._nextGroupId = 1

    def getGroups(self, proteinId):
        """Return a list of protein groups a protein is associated with."""
        return [self.groups[gId] for gId in self._proteinToGroupIds[proteinId]]

    def addProteinGroup(self, groupRepresentative):
        """Adds a new protein group and returns the groupId.

        The groupId is defined using an internal counter, which is incremented
        every time a protein group is added. The groupRepresentative is added
        as a leading protein.

        :param groupRepresentative: the protein group representing protein
        :returns: the protein groups groupId
        """
        groupId = self._getNextGroupId()
        self.groups[groupId] = ProteinGroup(groupId, groupRepresentative)
        self.addLeadingToGroups(groupRepresentative, groupId)
        return groupId

    def addLeadingToGroups(self, proteinIds, groupIds):
        """Add one or multiple leading proteins to one or multiple protein
        groups.

        :param proteinIds: a proteinId or a list of proteinIds, a proteinId
            must be a string.
        :param groupIds: a groupId or a list of groupIds, a groupId
            must be a string.
        """
        for groupId in AUX.toList(groupIds):
            self.groups[groupId].addLeadingProteins(proteinIds)
            self._addProteinIdsToGroupMapping(proteinIds, groupId)

    def addSubsetToGroups(self, proteinIds, groupIds):
        """Add one or multiple subset proteins to one or multiple protein
        groups.

        :param proteinIds: a proteinId or a list of proteinIds, a proteinId
            must be a string.
        :param groupIds: a groupId or a list of groupIds, a groupId
            must be a string.
        """
        for groupId in AUX.toList(groupIds):
            self.groups[groupId].addSubsetProteins(proteinIds)
            self._addProteinIdsToGroupMapping(proteinIds, groupId)

    def addSubsumableToGroups(self, proteinIds, groupIds):
        """Add one or multiple subsumable proteins to one or multiple protein
        groups.

        :param proteinIds: a proteinId or a list of proteinIds, a proteinId
            must be a string.
        :param groupIds: a groupId or a list of groupIds, a groupId
            must be a string.
        """
        for groupId in AUX.toList(groupIds):
            self.groups[groupId].addSubsumableProteins(proteinIds)
            self._addProteinIdsToGroupMapping(proteinIds, groupId)

    def _getNextGroupId(self):
        """Returns the next internal groupId and increments the counter."""
        groupId = self._nextGroupId
        self._nextGroupId += 1
        return str(groupId)

    def _addProteinIdsToGroupMapping(self, proteinIds, groupId):
        """Add a groupId to one or multiple entries of the internal
        proteinToGroupId mapping.

        :param proteinIds: a proteinId or a list of proteinIds, a proteinId
            must be a string.
        :param groupId: str, a groupId
        """
        for proteinId in AUX.toList(proteinIds):
            self._proteinToGroupIds[proteinId].add(groupId)


class ProteinGroup(object):
    """A protein amiguity group.

    :ivar proteins: set, all proteins that are associated with this group
    :ivar representative: a single protein that represents the group
    :ivar leading: set, leading proteins (MS:1002401)
    :ivar subset: a list of sub-set proteins
    :ivar subsumable: a list of sub-sumable proteins
    """
    def __init__(self, groupId, representative):
        self.id = groupId
        self.representative = representative
        self.proteins = set()
        self.leading = set()
        self.subset = set()
        self.subsumable = set()

    def _addProteins(self, proteinIds, containerName):
        """Add one or multiple proteinIds to the respective container.

        :param proteinIds: a proteinId or a list of proteinIds, a proteinId
            must be a string.
        :param containerName: must be either 'leading', 'subset' or 'subsumable'
        """
        proteinIds = AUX.toList(proteinIds)
        proteinContainer = getattr(self, containerName)
        proteinContainer.update(proteinIds)
        self.proteins.update(proteinIds)

    def addLeadingProteins(self, proteinIds):
        """Add one or multiple proteinIds as leading proteins.

        :param proteinIds: a proteinId or a list of proteinIds, a proteinId
            must be a string.
        """
        self._addProteins(proteinIds, 'leading')

    def addSubsetProteins(self, proteinIds):
        """Add one or multiple proteinIds as subset proteins.

        :param proteinIds: a proteinId or a list of proteinIds, a proteinId
            must be a string.
        """
        self._addProteins(proteinIds, 'subset')

    def addSubsumableProteins(self, proteinIds):
        """Add one or multiple proteinIds as subsumable proteins.

        :param proteinIds: a proteinId or a list of proteinIds, a proteinId
            must be a string.
        """
        self._addProteins(proteinIds, 'subsumable')


def mappingBasedGrouping(protToPeps):
    """Performs protein grouping based only on protein to peptide mappings.

    :param protToPeps: dict, for each protein (=key) contains a set of
        associated peptides (=value). For Example {protein: {peptide, ...}, ...}

    :Todo:
        Refactoring
            Processing of individual clusters can be put into a function
            Processing of groupInitiatingProteins, redundantProteins and
            subsetProteins can be put into function
            The First part that defines various sets of proteins can be put into
            a function

        Proteins
            Generate protein objects and characterize peptides

    returns a ProteinInference object
    """
    inference = ProteinInference(protToPeps)
    pepToProts = inference.pepToProts

    proteinClusters = _findProteinClusters(protToPeps, pepToProts)
    proteins = {}
    for clusterId, proteinCluster in enumerate(proteinClusters, 1):
        clusterProtToPeps = {p: protToPeps[p] for p in proteinCluster}

        #Find sameset proteins, define unique and non unique sameset proteins
        #NOTE: already unique proteins could be excluded to find sameset proteins
        samesetProteins = _findSamesetProteins(clusterProtToPeps)
        mergedProtToPeps = _mergeProteinEntries(samesetProteins,
                                                clusterProtToPeps)
        mergedPepToProts = _invertMapping(mergedProtToPeps)
        uniqueProteins = _findUniqueProteins(mergedPepToProts)
        remainingProteins = set(mergedProtToPeps).difference(uniqueProteins)

        # Remove subset proteins and check if remaining proteins become unique
        subsetProteinInfo = _findSubsetProteins(remainingProteins,
                                                mergedProtToPeps,
                                                mergedPepToProts)
        subsetProteins = [p for p, _ in subsetProteinInfo]
        subsetRemovedProtToPeps = _reducedProtToPeps(mergedProtToPeps,
                                                     subsetProteins)
        subsetRemovedPepToProts = _invertMapping(subsetRemovedProtToPeps)
        uniqueSubsetRemoved = _findUniqueProteins(subsetRemovedPepToProts)
        remainingProteins = remainingProteins.difference(subsetProteins)
        remainingProteins = remainingProteins.difference(uniqueSubsetRemoved)

        # Find redundant proteins #
        redundantProteins = _findRedundantProteins(subsetRemovedProtToPeps,
                                                   subsetRemovedPepToProts)
        remainingNonRedundant = remainingProteins.difference(redundantProteins)


        #Generate protein groups
        groupInitiatingProteins = uniqueSubsetRemoved.union(remainingNonRedundant)
        for protein in groupInitiatingProteins:
            proteinIds = AUX.toList(protein)

            groupId = inference.addProteinGroup(proteinIds[0])
            inference.addLeadingToGroups(proteinIds, groupId)

        #Add redundant proteins here (must be subsumable I guess)
        for protein in redundantProteins:
            proteinIds = AUX.toList(protein)

            connectedProteins = _mappingGetValueSet(
                mergedPepToProts, mergedProtToPeps[protein]
            )
            flatConnectedProteins = _flattenMergedProteins(connectedProteins)
            groupIds = _mappingGetValueSet(
                inference._proteinToGroupIds, flatConnectedProteins
            )
            inference.addSubsumableToGroups(proteinIds, groupIds)
            assert len(groupIds) > 1

        #Add subgroup proteins to the respective groups
        for protein, supersetProteins in subsetProteinInfo:
            proteinIds = AUX.toList(protein)

            flatSupersetProteins = _flattenMergedProteins(supersetProteins)
            superGroupIds = _mappingGetValueSet(
                inference._proteinToGroupIds, flatSupersetProteins
            )
            inference.addSubsetToGroups(proteinIds, superGroupIds)
            assert superGroupIds

    allProteins = reduce(lambda i,j: i.union(j),
                         [pg.proteins for pg in viewvalues(inference.groups)])
    assert len(allProteins) == len(protToPeps)
    return inference


def _findProteinClusters(protToPeps, pepToProts):
    """Find protein clusters in the specified protein to peptide mappings.

    A protein cluster is a group of proteins that are somehow directly or
    indirectly connected by shared peptides.

    :param protToPeps: dict, for each protein (=key) contains a set of
        associated peptides (=value). For Example {protein: {peptide, ...}, ...}
    :param pepToProts: dict, for each peptide (=key) contains a set of parent
        proteins (=value). For Example {peptide: {protein, ...}, ...}
    :returns: a list of protein clusters, each cluster is a set of proteins
    """
    clusters = list()
    resolvingProteins = set(protToPeps)
    while resolvingProteins:
        protein = resolvingProteins.pop()
        proteinCluster = set([protein])

        peptides = set(protToPeps[protein])
        parsedPeptides = set()

        while len(peptides) != len(parsedPeptides):
            for peptide in peptides:
                proteinCluster.update(pepToProts[peptide])
            parsedPeptides.update(peptides)

            for protein in proteinCluster:
                peptides.update(protToPeps[protein])
        clusters.append(proteinCluster)
        resolvingProteins = resolvingProteins.difference(proteinCluster)
    return clusters


def _findUniqueProteins(pepToProts):
    """Find proteins that are mapped to at least one unique peptide.

    :param pepToProts: dict, for each peptide (=key) contains a set of parent
        proteins (=value). For Example {peptide: {protein, ...}, ...}
    :returns: a set of unique proteins
    """
    uniqueProteins = list()
    for proteins in viewvalues(pepToProts):
        if len(proteins) == 1:
            uniqueProteins.extend(list(proteins))
    return set(uniqueProteins)


def _findSamesetProteins(protToPeps, proteins=None):
    """Find proteins that are mapped to an identical set of peptides.

    :param protToPeps: dict, for each protein (=key) contains a set of
        associated peptides (=value). For Example {protein: {peptide, ...}, ...}
    :param proteins: iterable, proteins that are tested for having equal
        evidence. If not specified all proteins are tested
    :returns: a list of sorted protein tuples that share equal peptide evidence
    """
    proteins = viewkeys(protToPeps) if proteins is None else proteins

    equalEvidence = ddict(set)
    for protein in proteins:
        peptides = protToPeps[protein]
        equalEvidence[tuple(sorted(peptides))].add(protein)
    equalProteins = list()
    for proteins in viewvalues(equalEvidence):
        if len(proteins) > 1:
            equalProteins.append(tuple(sorted(proteins)))
    return equalProteins


def _findSubsetProteins(proteins, protToPeps, pepToProts):
    """Find proteins which peptides are a sub-set, but not a same-set to other
    proteins.

    :param proteins: iterable, proteins that are tested for being a subset
    :param pepToProts: dict, for each peptide (=key) contains a set of parent
        proteins (=value). For Example {peptide: {protein, ...}, ...}
    :param protToPeps: dict, for each protein (=key) contains a set of
        associated peptides (=value). For Example {protein: {peptide, ...}, ...}
    :returns: a list of pairs of protein and their superset proteins.
        [(protein, {superset protein, ...}), ...]
    """
    proteinsEqual = lambda prot1, prot2: protToPeps[prot1] == protToPeps[prot2]

    subGroups = list()
    for protein in proteins:
        peptideCounts = Counter()
        for peptide in protToPeps[protein]:
            proteins = pepToProts[peptide]
            peptideCounts.update(proteins)
        peptideCount = peptideCounts.pop(protein)

        superGroups = set()
        for sharingProtein, sharedPeptides in peptideCounts.most_common():
            if peptideCount == sharedPeptides:
                if not proteinsEqual(protein, sharingProtein):
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
        descending order according to their protein names. For example in the
        case of a tie between proteins "A" and "B", protein "B" would be
        removed.
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
    for protein in sorted(proteins, reverse=True):
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


def _mergeProteinEntries(proteinLists, protToPeps):
    """Returns a new "protToPeps" dictionary with entries merged that are
    present in proteinLists.

    NOTE:
        The key of the merged entry is a tuple of the sorted protein keys. This
        behaviour might change in the future; the tuple might be replaced by
        simply one of the protein entries which is then representative for all.

    :param proteinLists: a list of protein groups that will be merged
        [{protein, ...}, ...]
    :param protToPeps: dict, for each protein (=key) contains a set of
        associated peptides (=value). For Example {protein: {peptide, ...}, ...}
    :returns: dict, {protein: set([peptid, ...])}
    """
    mergedProtToPeps = dict(protToPeps)
    for proteins in proteinLists:
        for protein in proteins:
            peptides = mergedProtToPeps.pop(protein)
        mergedProtein = tuple(sorted(proteins))
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

    :param mapping: dict, for each key contains a set of entries

    :returns: an inverted mapping that each entry of the values points to a set
        of initial keys.
    """
    invertedMapping = ddict(set)
    for key, values in viewitems(mapping):
        for value in values:
            invertedMapping[value].add(key)
    return invertedMapping


def _getValueCounts(mapping):
    """Returns a counter object; contains for each key of the mapping the counts
    of the respective value element (= set length).

    :param mapping: dict, for each key contains a set of entries.
    :returns: a counter
    """
    return Counter({k: len(v) for k, v in viewitems(mapping)})


def _mappingGetValueSet(mapping, keys):
    """Return a combined set of values from the mapping.

    :param mapping: dict, for each key contains a set of entries

    returns a set of combined entries
    """
    setUnion = lambda i,j: i.union(j)
    return reduce(setUnion, [mapping[k] for k in keys])


def _flattenMergedProteins(proteins):
    """Return a set where merged protein entries in proteins are flattened.

    :param proteins: an iterable of proteins, can contain merged protein entries
        in the form of tuple([protein1, protein2]).
    returns a set of protein entries, where all entries are strings
    """
    proteinSet = set()
    for protein in proteins:
        if isinstance(protein, tuple):
            proteinSet.update(protein)
        else:
            proteinSet.add(protein)
    return proteinSet


# - Deprecated early version of grouping - #
"""
def makeInference(protToPeps):
    proteinInference = ProteinInference(protToPeps)
    summary = doTheGrouping(protToPeps, proteinInference.pepToProts)
    proteinInference.groups = summary['proteinGroups']
    proteinInference._proteinToGroupIds = summary['proteinToGroupIds']
    return proteinInference


def getGroupIds(idLookup, entries):
    """ """
    setUnion = lambda i,j: i.union(j)
    return reduce(setUnion, [idLookup[e] for e in entries])


def doTheGrouping(protToPeps, pepToProts):
    proteinClusters = _findProteinClusters(protToPeps, pepToProts)
    proteinGroups = {}
    proteins = {}
    proteinToGroupIds = ddict(set)
    currGroupId = 1
    for clusterId, proteinCluster in enumerate(proteinClusters, 1):
        clusterProtToPeps = {p: protToPeps[p] for p in proteinCluster}

        #Find sameset proteins, define unique and non unique sameset proteins
        #NOTE: already unique proteins could be excluded to find sameset proteins
        samesetProteins = _findSamesetProteins(clusterProtToPeps)
        mergedProtToPeps = _mergeProteinEntries(samesetProteins,
                                                clusterProtToPeps)
        mergedPepToProts = _invertMapping(mergedProtToPeps)
        uniqueProteins = _findUniqueProteins(mergedPepToProts)
        remainingProteins = set(mergedProtToPeps).difference(uniqueProteins)

        # Remove subset proteins and check if remaining proteins become unique
        subsetProteinInfo = _findSubsetProteins(remainingProteins,
                                                mergedProtToPeps,
                                                mergedPepToProts)
        subsetProteins = [p for p, _ in subsetProteinInfo]
        subsetRemovedProtToPeps = _reducedProtToPeps(mergedProtToPeps,
                                                     subsetProteins)
        subsetRemovedPepToProts = _invertMapping(subsetRemovedProtToPeps)
        uniqueSubsetRemoved = _findUniqueProteins(subsetRemovedPepToProts)
        remainingProteins = remainingProteins.difference(subsetProteins)
        remainingProteins = remainingProteins.difference(uniqueSubsetRemoved)

        # Find redundant proteins #
        redundantProteins = _findRedundantProteins(subsetRemovedProtToPeps,
                                                   subsetRemovedPepToProts)
        remainingNonRedundant = remainingProteins.difference(redundantProteins)


        # - Generate protein groups and sssign protein to groups - #
        groupInitiatingProteins = uniqueSubsetRemoved.union(remainingNonRedundant)
        for protein in groupInitiatingProteins:
            proteinIds = AUX.toList(protein)

            proteinGroup = ProteinGroup(currGroupId, proteinIds[0])
            proteinGroup.addLeadingProteins(proteinIds)
            proteinGroups[proteinGroup.id] = proteinGroup

            #Code duplication, nearly#
            for proteinId in proteinIds:
                proteinToGroupIds[proteinId].add(currGroupId)
                proteinToGroupIds[protein].add(currGroupId)
            #-------#
            currGroupId += 1

        #Add redundant proteins here (must be subsumable I guess)
        for protein in redundantProteins:
            proteinIds = AUX.toList(protein)

            connectedProteins = set()
            for peptide in mergedProtToPeps[protein]:
                connectedProteins.update(mergedPepToProts[peptide])
            superGroupIds = getGroupIds(proteinToGroupIds, connectedProteins)
            #-Tests-#
            assert superGroupIds
            assert len(superGroupIds) > 1
            #-------#

            #Code duplication, nearly#
            for groupId in superGroupIds:
                proteinGroups[groupId].addSubsumableProteins(proteinIds)
            for proteinId in proteinIds:
                proteinToGroupIds[proteinId].update(superGroupIds)
            proteinToGroupIds[protein].update(superGroupIds)
            #-------#

        #Add subgroup proteins to the respective groups
        for protein, supersetProteins in subsetProteinInfo:
            proteinIds = AUX.toList(protein)

            superGroupIds = getGroupIds(proteinToGroupIds, supersetProteins)
            #-Tests-#
            assert superGroupIds
            #-------#

            #Code duplication, nearly#
            for groupId in superGroupIds:
                proteinGroups[groupId].addSubsetProteins(proteinIds)
            for proteinId in proteinIds:
                proteinToGroupIds[proteinId].update(superGroupIds)
            #-------#

    allProteins = set()
    for proteinGroup in viewvalues(proteinGroups):
        allProteins.update(proteinGroup.proteins)
    #assert len(allProteins) == len(protToPeps)

    summary = {'proteinGroups': proteinGroups, 'proteins': proteins,
               'proteinToGroupIds': proteinToGroupIds}
    return summary
"""
#"""
