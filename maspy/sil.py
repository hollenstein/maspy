from __future__ import print_function, division
from future.utils import viewkeys, viewvalues, viewitems, listvalues, listitems

import itertools

import maspy.auxiliary as aux
import maspy.constants
import maspy.peptidemethods

# --- Functions dealing with isotopic labels --- #
class LabelDescriptor(object):
    """Describes a MS1 stable isotope label setup for quantification.

    :ivar labels: Contains a dictionary with all possible label states, keys (=labelStates) are increasing integers starting from 0
    :ivar excludingModifictions: bool, set to True if any label has specified excludingModifications
    """
    def __init__(self):
        self.labels = dict()
        self.excludingModifictions = False
        self._labelCounter = 0

    def addLabel(self, aminoAcidLabels, excludingModifications=None):
        """Adds a new labelstate.

        :ivar aminoAcidsLabels: Describes which amino acids can bear which labels
        possible keys amino acids in one letter code and ('nTerm', 'cTerm')
        possible values are keys from :var:`maspy.constants.aaModMass` as strings or list of strings
        eg. {'nTerm':'u:188', 'K':['u:188', 'u:188']} for one expected label at the nterminus and two expected labels at Lysine
        :ivar excludingModifications: Describes which modifications can prevent the addition of labels
        keys and values have to be keys from :var:`maspy.constants.aaModMass` written as a string.
        eg. {'u:1':'u:188'} For each modification 'u:1' that is present at an amino acid or terminus of a peptide
        the number of expected labels at this position is reduced by one
        """
        if excludingModifications is not None:
            self.excludingModifictions = True

        self.labels[self._labelCounter] = dict()
        self.labels[self._labelCounter]['aminoAcidLabels'] = aminoAcidLabels
        self.labels[self._labelCounter]['excludingModifications'] = excludingModifications
        self._labelCounter += 1


def returnLabelStateMassDifferences(peptide, labelDescriptor, labelState=None, sequence=None):
    """Calculates the mass difference for alternative possible label states of a given peptide.

    :ivar peptide: Peptide to calculate alternative label states
    :ivar labelDescriptor: :class:`LabelDescriptor` describes the label setup of an experiment
    :ivar labelState: label state of the peptide, if None it is calculated by :func:`returnLabelState`
    :ivar sequence: unmodified amino acid sequence of :var:`peptide`, if None it is calculated with :func:`removeModifications`

    :return: {alternativeLabelSate: massDifference, ...} or {} if the peptide label state is -1
    (massDifference + peptide mass = expected mass of alternatively labeled peptide)

    See also :class:`LabelDescriptor`, :func:`returnLabelState`
    """
    labelState = returnLabelState(peptide, labelDescriptor) if labelState is None else labelState
    sequence = maspy.peptidemethods.removeModifications(peptide) if sequence is None else sequence

    if labelState < 0:
        # special case for mixed label... #
        return dict()

    # define type and number of labels of the peptide
    labelModNumbers = dict()
    for labelStateModList in viewvalues(expectedLabelPosition(peptide, labelDescriptor.labels[labelState], sequence=sequence)):
        for labelMod in labelStateModList:
            labelModNumbers.setdefault(labelMod, int())
            labelModNumbers[labelMod] += 1

    # calculate the combined labels mass of the peptide
    labelMass = int()
    for labelMod, modCounts in viewitems(labelModNumbers):
        labelMass += maspy.constants.aaModMass[labelMod] * modCounts

    # calculate mass differences to all other possible label states
    labelStateMassDifferences = dict()
    for possibleLabelState in viewkeys(labelDescriptor.labels):
        if possibleLabelState == labelState:
            continue

        labelModNumbers = dict()
        for labelStateModList in viewvalues(expectedLabelPosition(peptide, labelDescriptor.labels[possibleLabelState], sequence=sequence)):
            for labelMod in labelStateModList:
                labelModNumbers.setdefault(labelMod, int())
                labelModNumbers[labelMod] += 1

        possibleLabelMass = int()
        for labelMod, modCounts in viewitems(labelModNumbers):
            possibleLabelMass += maspy.constants.aaModMass[labelMod] * modCounts

        possibleLabelMassDifference = possibleLabelMass - labelMass
        labelStateMassDifferences[possibleLabelState] = possibleLabelMassDifference
    return labelStateMassDifferences


def returnLabelState(peptide, labelDescriptor, labelSymbols=None, labelAminoacids=None):
    """Calculates the label state of a given peptide for the label setup described in labelDescriptor

    :ivar peptide: peptide which label state should be calcualted
    :ivar labelDescriptor: :class:`LabelDescriptor` describes the label setup of an experiment
    :ivar labelSymbols: modifications that show a label, calculated by :func:`modSymbolsFromLabelInfo`
    :ivar labelAminoacids: amino acids that can bear a label, calculated by :func:`modAminoacidsFromLabelInfo`

    :return: integer that shows the label state
    >=0: predicted label state of the peptide
     -1: peptide sequence can't bear any labelState modifications
     -2: peptide modifications don't fit to any predicted labelState
     -3: peptide modifications fit to a predicted labelState, but not all predicted labelStates are distinguishable
    """
    labelSymbols = modSymbolsFromLabelInfo(labelDescriptor) if labelSymbols is None else labelSymbols
    labelAminoacids = modAminoacidsFromLabelInfo(labelDescriptor) if labelAminoacids is None else labelAminoacids

    sequence = maspy.peptidemethods.removeModifications(peptide)
    modPositions = maspy.peptidemethods.returnModPositions(peptide, indexStart=0, removeModString=False)

    labelState = None

    # No amino acids in sequence which can bear a label modification (ignores presence of excluding modifications)
    if all([(True if sequence.find(labelAminoacid) == -1 else False) for labelAminoacid in labelAminoacids]):
        # No terminal label modifications specified by labelDescriptor
        if 'nTerm' not in labelAminoacids and 'cTerm' not in labelAminoacids:
            labelState = -1

    # Check if the peptide mofidifcations fit to any predicted label state
    if labelState is None:
        peptideLabelPositions = dict()
        for labelSymbol in labelSymbols:
            if labelSymbol in viewkeys(modPositions):
                for sequencePosition in modPositions[labelSymbol]:
                    peptideLabelPositions.setdefault(sequencePosition, list())
                    peptideLabelPositions[sequencePosition].append(labelSymbol)
        for sequencePosition in list(viewkeys(peptideLabelPositions)):
            peptideLabelPositions[sequencePosition] = sorted(peptideLabelPositions[sequencePosition])

        predictedLabelStates = dict()
        for predictedLabelState, labelStateInfo in viewitems(labelDescriptor.labels):
            expectedLabelMods = expectedLabelPosition(peptide, labelStateInfo, sequence=sequence, modPositions=modPositions)
            predictedLabelStates[predictedLabelState] = expectedLabelMods
            if peptideLabelPositions == expectedLabelMods:
                # If another expectedLabel state has already been matched, there is an ambiguity between label states...
                labelState = predictedLabelState

    if labelState is None:
        # Peptide mofidifcations don't fit to any predicted label state
        labelState = -2
    elif labelState != -1:
        # Check if all predicted label states are distinguishable
        for labelState1, labelState2 in set(itertools.combinations(range(len(predictedLabelStates)), 2)):
            if predictedLabelStates[labelState1] == predictedLabelStates[labelState2]:
                labelState = -3
                break

    return labelState


def modSymbolsFromLabelInfo(labelDescriptor):
    """Returns a set of all modiciation symbols which were used in the labelDescriptor

    :ivar labelDescriptor: :class:`LabelDescriptor` describes the label setup of an experiment
    """
    modSymbols = set()
    for labelStateEntry in viewvalues(labelDescriptor.labels):
        for labelPositionEntry in viewvalues(labelStateEntry['aminoAcidLabels']):
            for modSymbol in aux.toList(labelPositionEntry):
                if modSymbol != '':
                    modSymbols.add(modSymbol)
    return modSymbols


def modAminoacidsFromLabelInfo(labelDescriptor):
    """Returns a set of all amino acids and termini which can bear a label, as described in labelDescriptor

    :ivar labelDescriptor: :class:`LabelDescriptor` describes the label setup of an experiment
    """
    modAminoacids = set()
    for labelStateEntry in viewvalues(labelDescriptor.labels):
        for labelPositionEntry in viewkeys(labelStateEntry['aminoAcidLabels']):
            for modAminoacid in aux.toList(labelPositionEntry):
                if modAminoacid != '':
                    modAminoacids.add(modAminoacid)
    return modAminoacids


def expectedLabelPosition(peptide, labelStateInfo, sequence=None, modPositions=None):
    """Returns a modification description of a certain label state of a peptide.

    :ivar peptide: Peptide sequence used to calculat the expected label state modifications
    :ivar labelStateInfo: An entry of :attr:`LabelDescriptor.labels` that describes a label state
    :ivar sequence: unmodified amino acid sequence of :var:`peptide`, if None it is calculated with :func:`removeModifications`
    :ivar modPositions: dictionary describing the modification state of :var:`peptide`,
    if None it is calculated with :func:`returnModPositions`

    :return: {sequence position: sorted list of expected label modifications on that position, ...}
    """
    modPositions = maspy.peptidemethods.returnModPositions(peptide, indexStart=0) if modPositions is None else modPositions
    sequence = maspy.peptidemethods.removeModifications(peptide) if sequence is None else sequence

    currLabelMods = dict()
    for labelPosition, labelSymbols in viewitems(labelStateInfo['aminoAcidLabels']):
        labelSymbols = aux.toList(labelSymbols)
        if labelSymbols == ['']:
            pass
        elif labelPosition == 'nTerm':
            currLabelMods.setdefault(0, list())
            currLabelMods[0].extend(labelSymbols)
        else:
            for sequencePosition in aux.findAllSubstrings(sequence, labelPosition):
                currLabelMods.setdefault(sequencePosition, list())
                currLabelMods[sequencePosition].extend(labelSymbols)

    if labelStateInfo['excludingModifications'] is not None:
        for excludingModification, excludedLabelSymbol in viewitems(labelStateInfo['excludingModifications']):
            if excludingModification in modPositions:
                for excludingModPosition in modPositions[excludingModification]:
                    if excludingModPosition in currLabelMods:
                        if excludedLabelSymbol in currLabelMods[excludingModPosition]:
                            if len(currLabelMods[excludingModPosition]) == 1:
                                del(currLabelMods[excludingModPosition])
                            else:
                                excludedModIndex = currLabelMods[excludingModPosition].index(excludedLabelSymbol)
                                currLabelMods[excludingModPosition].pop(excludedModIndex)

    for sequencePosition in list(viewkeys(currLabelMods)):
        currLabelMods[sequencePosition] = sorted(currLabelMods[sequencePosition])
    return currLabelMods
