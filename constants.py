import copy

# Define constants #
atomicMassH = 1.00782504
atomicMassProton = 1.00727646677

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
xTandemMassToUniModDict = dict([(round(mass, 5), unimod) for unimod, mass in xTandemMassToUniModDict.items()])

"""
This dict contains regular expressions for cleavage rules of the most popular proteolytic enzymes.
The rules were copied from pyteomics http://pythonhosted.org/pyteomics/ and originally taken from
`PeptideCutter tool <http://ca.expasy.org/tools/peptidecutter/peptidecutter_enzymes.html>`_ at Expasy.
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

mzmlAccessions = dict()
mzmlAccessions['MS:1000927'] = {'name':'iit', 'msLevel':None}
mzmlAccessions['MS:1000285'] = {'name':'tic', 'msLevel':None}
mzmlAccessions['MS:1000016'] = {'name':'rt', 'msLevel':None}
mzmlAccessions['MS:1000827'] = {'name':'targetWindowMz', 'msLevel':2}
mzmlAccessions['MS:1000828'] = {'name':'lowWindowOffset', 'msLevel':2}
mzmlAccessions['MS:1000829'] = {'name':'highWindowOffset', 'msLevel':2}
mzmlAccessions['MS:1000744'] = {'name':'obsMz', 'msLevel':2}
mzmlAccessions['MS:1000042'] = {'name':'obsI', 'msLevel':2}
mzmlAccessions['MS:1000041'] = {'name':'charge', 'msLevel':2}
mzmlAccessions['MS:1000501'] = {'name':'lowWindowLimit', 'msLevel':2}
mzmlAccessions['MS:1000500'] = {'name':'highWindowLimit', 'msLevel':2}
mzmlAccessions['MS:1000504']  = {'name':'basePeakMz', 'msLevel':None}
mzmlAccessions['MS:1000505']  = {'name':'basePeakI', 'msLevel':None}

#mzmlAccessions['MS:1000512'] = {'name':'filterString', 'msLevel':None}

