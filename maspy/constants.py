from __future__ import print_function, division, unicode_literals
from future.utils import viewkeys, viewvalues, viewitems, listvalues, listitems

try: # python 2.7
    from itertools import izip as zip
except ImportError: # python 3.x series
    pass
###############################################################################
import copy

import pyteomics.mass
""" for details see http://pythonhosted.org/pyteomics/ """

from pyteomics.mass import Composition as COMPOSITION
""" A Composition object stores a chemical composition of a
    substance. Basically, it is a dict object, with the names
    of chemical elements as keys and values equal to an integer number of
    atoms of the corresponding element in a substance.

    The main improvement over dict is that Composition objects allow
    adding and subtraction. -> For details see pyteomics.mass.Composition
"""

aaComp = dict()
""" A dictionary with elemental compositions of the twenty standard amino acid residues.
    This concept was inherited from pyteomics.mass.std_aa_comp
"""
aaComp.update({'A':COMPOSITION({'H': 5, 'C': 3, 'O': 1, 'N': 1}),
               'C':COMPOSITION({'H': 5, 'C': 3, 'S': 1, 'O': 1, 'N': 1}),
               'D':COMPOSITION({'H': 5, 'C': 4, 'O': 3, 'N': 1}),
               'E':COMPOSITION({'H': 7, 'C': 5, 'O': 3, 'N': 1}),
               'F':COMPOSITION({'H': 9, 'C': 9, 'O': 1, 'N': 1}),
               'G':COMPOSITION({'H': 3, 'C': 2, 'O': 1, 'N': 1}),
               'H':COMPOSITION({'H': 7, 'C': 6, 'N': 3, 'O': 1}),
               'I':COMPOSITION({'H': 11, 'C': 6, 'O': 1, 'N': 1}),
               'K':COMPOSITION({'H': 12, 'C': 6, 'N': 2, 'O': 1}),
               'L':COMPOSITION({'H': 11, 'C': 6, 'O': 1, 'N': 1}),
               'M':COMPOSITION({'H': 9, 'C': 5, 'S': 1, 'O': 1, 'N': 1}),
               'N':COMPOSITION({'H': 6, 'C': 4, 'O': 2, 'N': 2}),
               'P':COMPOSITION({'H': 7, 'C': 5, 'O': 1, 'N': 1}),
               'Q':COMPOSITION({'H': 8, 'C': 5, 'O': 2, 'N': 2}),
               'R':COMPOSITION({'H': 12, 'C': 6, 'N': 4, 'O': 1}),
               'S':COMPOSITION({'H': 5, 'C': 3, 'O': 2, 'N': 1}),
               'T':COMPOSITION({'H': 7, 'C': 4, 'O': 2, 'N': 1}),
               'V':COMPOSITION({'H': 9, 'C': 5, 'O': 1, 'N': 1}),
               'W':COMPOSITION({'C': 11, 'H': 10, 'N': 2, 'O': 1}),
               'Y':COMPOSITION({'H': 9, 'C': 9, 'O': 2, 'N': 1})
              })
aaMass = dict([(name, pyteomics.mass.calculate_mass(comp)) for name, comp in viewitems(aaComp)])

aaModComp = dict()
""" A dictionary with elemental compositions of the peptide modifications.
    Modifications present at "www.unimod.org" should be written as "u:X", where
    X is the unimod accession number. If a modification is not present an text
    abbriviation should be used. This concept was inherited from
    ``pyteomics.mass.std_aa_comp``.

    TODO: in the future this table should be imported from two external files.
    The first is directly obtained from www.unimod.org, the second contains
    user specified entries. It is also possible to specify a modification folder
    where multiple user specified files can be deposited for importing.
"""
aaModComp.update({'u:1':COMPOSITION({'C': 2, 'H': 2, 'O': 1}),
                  'u:3':COMPOSITION({'C': 10, 'H': 14, 'N': 2, 'O': 2, 'S': 1}),
                  'u:4':COMPOSITION({'C': 2, 'H': 3, 'N': 1, 'O': 1}),
                  'u:5':COMPOSITION({'C': 1, 'H': 1, 'N': 1, 'O': 1}),
                  'u:7':COMPOSITION({'H': -1, 'N': -1, 'O': 1}),
                  'u:21':COMPOSITION({'H': 1, 'O': 3, 'P': 1}),
                  'u:27':COMPOSITION({'H': -2, 'O': -1}),
                  'u:28':COMPOSITION({'H': -3, 'N': -1}),
                  'u:34':COMPOSITION({'C': 1, 'H': 2}),
                  'u:35':COMPOSITION({'O': 1}),
                  'u:36':COMPOSITION({'C': 2, 'H': 4}),
                  'u:121':COMPOSITION({'C': 4, 'H': 6, 'N': 2, 'O': 2}),
                  'u:188':COMPOSITION({'C': -6, 'C[13]': 6}),
                  'u:199':COMPOSITION({'H[2]': 4, 'C': 2}),
                  'u:374':COMPOSITION({'H': -1}),
                  'u:1020':COMPOSITION({'C': 8, 'H': 12, 'O': 3}),
                  'u:1356':COMPOSITION({'C': 5, 'H': 9, 'O': 7, 'P': 1}),
                  'DSS':COMPOSITION({'C': 8, 'H': 10, 'O': 2}),
                  '*':COMPOSITION({})
                 })
aaModMass = dict([(name, pyteomics.mass.calculate_mass(comp)) for name, comp in viewitems(aaModComp)])
#TODO change all modification instances from "UNIMOD:X" to "u:X"
for accession, composition in listitems(aaModMass):
    if accession.startswith('u:'):
        aaModMass[accession.replace('u:', 'UNIMOD:')] = composition

# Define constants #
atomicMassH = 1.00782504
atomicMassProton = 1.00727646677

""" TODO: Maybe this is still needed to import xtandem results and convert modification mass strings to "aaModMass" keys...
unimodToMassDict = dict()
unimodToMassDict['34'] = 14.015650 # Methylation
unimodToMassDict['36'] = 28.031300 # Dimethyl light label
unimodToMassDict['199'] = 32.056407 # Dimethyl medium label
unimodToMassDict['4'] = 57.021464 # Carbamidomethylation
unimodToMassDict['374'] = -1.007825 # Half of a disulfide bridge
unimodToMassDict['7'] = 0.984016 # Deamidated
unimodToMassDict['188'] = 6.020129 # Label:13C(6)
unimodToMassDict['35'] = 15.994915 # Oxidation
unimodToMassDict['21'] = 79.966331 # Phospho
unimodToMassDict['1'] = 42.010565 # Acetyl
unimodToMassDict['27'] = -18.010565 # Glu->pyro-Glu
unimodToMassDict['28'] = -17.026549 # Gln->pyro-Glu
unimodToMassDict['121'] = 114.042927 # GG, ubiquitinlyation residue
unimodToMassDict['DSS'] = 138.068 # Xlink:DSS / BS3
unimodToMassDict['1020'] = 156.078644 # Xlink:DSS, Water-quenched monolink of of DSS/BS3 crosslinker
unimodToMassDict['1356'] = 212.008590 # phosphate-ribosylation: R, (D, E)
unimodToMassDict['213'] = 541.061110 # ADP-Ribosyl, R
unimodToMassDict['5'] = 43.005814 # Carbamyl, pep-n / K / R
unimodToMassDict['3'] = 226.077598 # Biotin K
unimodToMassDict['*'] = 0.0 #Place holder for the second position of a dipeptide modification like a crosslink.

xTandemMassToUniModDict = copy.deepcopy(unimodToMassDict)
xTandemMassToUniModDict[4] = 57.02147
xTandemMassToUniModDict[374] = -1.00783
xTandemMassToUniModDict[1] = 42.01057
xTandemMassToUniModDict = dict([(round(mass, 5), unimod) for unimod, mass in viewitems(xTandemMassToUniModDict)])
"""

expasy_rules = dict()
""" The dictionary expasy_rules contains regular expressions for cleavage rules
    of the most popular proteolytic enzymes. The rules were copied from
    `Pyteomics <http://pythonhosted.org/pyteomics/>`_ and initially taken from
    the PeptideCutter tool at `Expasy
    <http://ca.expasy.org/tools/peptidecutter/peptidecutter_enzymes.html>`_
"""
expasy_rules = {'arg-c': 'R',
                'asp-n': '\\w(?=D)',
                'bnps-skatole': 'W',
                'caspase 1': '(?<=[FWYL]\\w[HAT])D(?=[^PEDQKR])',
                'caspase 10': '(?<=IEA)D',
                'caspase 2': '(?<=DVA)D(?=[^PEDQKR])',
                'caspase 3': '(?<=DMQ)D(?=[^PEDQKR])',
                'caspase 4': '(?<=LEV)D(?=[^PEDQKR])',
                'caspase 5': '(?<=[LW]EH)D',
                'caspase 6': '(?<=VE[HI])D(?=[^PEDQKR])',
                'caspase 7': '(?<=DEV)D(?=[^PEDQKR])',
                'caspase 8': '(?<=[IL]ET)D(?=[^PEDQKR])',
                'caspase 9': '(?<=LEH)D',
                'chymotrypsin high specificity': '([FY](?=[^P]))|(W(?=[^MP]))',
                'chymotrypsin low specificity': '([FLY](?=[^P]))|(W(?=[^MP]))|(M(?=[^PY]))|(H(?=[^DMPW]))',
                'clostripain': 'R',
                'cnbr': 'M',
                'enterokinase': '(?<=[DE]{3})K',
                'factor xa': '(?<=[AFGILTVM][DE]G)R',
                'formic acid': 'D',
                'glutamyl endopeptidase': 'E',
                'granzyme b': '(?<=IEP)D',
                'hydroxylamine': 'N(?=G)',
                'iodosobenzoic acid': 'W',
                'lysc': 'K',
                'ntcb': '\\w(?=C)',
                'pepsin ph1.3': '((?<=[^HKR][^P])[^R](?=[FLWY][^P]))|((?<=[^HKR][^P])[FLWY](?=\\w[^P]))',
                'pepsin ph2.0': '((?<=[^HKR][^P])[^R](?=[FL][^P]))|((?<=[^HKR][^P])[FL](?=\\w[^P]))',
                'proline endopeptidase': '(?<=[HKR])P(?=[^P])',
                'proteinase k': '[AEFILTVWY]',
                'staphylococcal peptidase i': '(?<=[^E])E',
                'thermolysin': '[^DE](?=[AFILMV])',
                'thrombin': '((?<=G)R(?=G))|((?<=[AFGILTVM][AFGILTVWA]P)R(?=[^DE][^DE]))',
                'trypsin': '([KR](?=[^P]))|((?<=W)K(?=P))|((?<=M)R(?=P))',
                'trypsin simple': '[KR]'
                }
