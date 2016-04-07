from __future__ import print_function, division, unicode_literals
from future.utils import viewkeys, viewvalues, viewitems, listvalues, listitems

try: # python 2.7
    from itertools import izip as zip
except ImportError: # python 3.x series
    pass
###############################################################################

import os
import sys
sys.path.append('D:/Dropbox/python/maspy')
sys.path.append('C:/Users/David/Dropbox/python/maspy')

import maspy.new.auxiliary as aux
import maspy.new
import maspy.new.import_export
import maspy.new.core
import maspy.new.xml


# --- Testing mzml related functionalities (core, import_export, xml) --- #
import io
from lxml import etree as etree

def compareParams(paramList, referenceParamList):
    referenceDict = {'refParamGroups': list()}
    for referenceParam in referenceParamList:
        if len(referenceParam) < 3:
            referenceDict['refParamGroups'].append(referenceParam[1])
        else:
            referenceDict[referenceParam[0]] = referenceParam[1:]
    observedRefParamGroups = list()

    for param in paramList:
        if len(param) == 3:
            paramType = 'cvParam'
        elif len(param) == 4:
            paramType = 'userParam'
        else:
            paramType = 'refParamGroup'

        if paramType in ['cvParam', 'userParam']:
            if param[0] not in referenceDict:
                raise Exception(param)
            try:
                if round(float(param[1]), 6) != round(float(referenceDict[param[0]][0]), 6):
                    raise Exception(param)
            except ValueError:
                if param[1] != referenceDict[param[0]][0]:
                    raise Exception(param)
            if param[2] != referenceDict[param[0]][1]:
                raise Exception(param)
        if paramType == 'userParam':
            if param[3] != referenceDict[param[0]][2]:
                raise Exception(param)
        if paramType == 'refParamGroup':
            observedRefParamGroups.append(param[1])
    if observedRefParamGroups != referenceDict['refParamGroups']:
        raise Exception(observedRefParamGroups)
    return True


testfilepath = aux.joinpath(os.path.dirname(aux.__file__), 'testdata', 'spectrum.xml')
with io.open(testfilepath, 'r', encoding='utf-8') as openfile:
    root = etree.XML(openfile.read())

#TESTING mzml.smiFromXmlSpectrum(), mzml.extractBinaries()
spectra = list()
for xmlSpectrum in root.getchildren():
    smi, binaryDataArrayList = maspy.new.import_export.smiFromXmlSpectrum(xmlSpectrum, 'test')
    sai = maspy.new.core.Sai(smi.id, smi.specfile)
    sai.arrays, sai.arrayInfo = maspy.new.xml.extractBinaries(binaryDataArrayList,
                                                     smi.attributes['defaultArrayLength'])
    spectra.append({'smi': smi, 'sai': sai})


#TESTING Smi, MzmlScan, MzmlPrecursor, TODO: add MzmlProduct
smi = spectra[0]['smi']
smi.attributes == {'dataProcessingRef': 'adataprocessingref',
                   'defaultArrayLength': '870',
                   'id': 'controllerType=0 controllerNumber=1 scan=4',
                   'index': '3',
                   'sourceFileRef': 'aref',
                   'spotID': 'aspotid'
                   }
smi.id == '4'
smi.specfile == 'test'
len(smi.scanList) == 1
len(smi.scanList[0].params) == 4

smi.scanListParams == [('MS:1000795', '', None)]
scan = smi.scanList[0]
scan.scanWindowList == ((('MS:1000501', '300.0', 'MS:1000040'), ('MS:1000500', '1650.0', 'MS:1000040')),)

scanReferenceParams = [('MS:1000016', '0.018546918', 'UO:0000031'),
                       ('MS:1000512', 'FTMS + p NSI Full ms [300.00-1650.00]', None),
                       ('MS:1000616', '1', None),
                       ('MS:1000927', '19.999999552965', 'UO:0000028')
                       ]
compareParams(scan.params, scanReferenceParams)

smi = spectra[1]['smi']
len(smi.precursorList) == 1
smiReferenceParams = [('MS:1000580', '', None),
                      ('MS:1000511', '2', None),
                      ('MS:1000130', '', None),
                      ('MS:1000127', '', None),
                      ('MS:1000504', '251.1017004', 'MS:1000040'),
                      ('MS:1000505', '1.7075384e05', 'MS:1000131'),
                      ('MS:1000285', '6.3746325e05', None),
                      ('MS:1000528', '61.927135467529', 'MS:1000040'),
                      ('MS:1000527', '543.267944335938', 'MS:1000040'),
                      ('userParam1', 'a simple user param test case', None, None),
                      ('userParam2', '100', 'UO:0000029', 'xsd:float'),
                      ('ref', 'fakeparamgroup')
                      ]
compareParams(smi.params, smiReferenceParams)

precursor = smi.precursorList[0]
precursor.activation == (('MS:1000422', '', None), ('MS:1000045', '27.0', 'UO:0000266'))
precursor.spectrumRef == 'controllerType=0 controllerNumber=1 scan=4'
isolationWindowReferenceParams = (('MS:1000827', '410.68', 'MS:1000040'),
                                  ('MS:1000828', '0.800000011921', 'MS:1000040'),
                                  ('MS:1000829', '0.800000011921', 'MS:1000040')
                                  )
compareParams(precursor.isolationWindow, isolationWindowReferenceParams)
len(precursor.selectedIonList) == 1
selectedIon = precursor.selectedIonList[0]
selectedIonReferenceParams = (('MS:1000744', '410.681187454293', 'MS:1000040'),
                              ('MS:1000041', '2', None),
                              ('MS:1000042', '7.2570925e05', 'MS:1000131')
                              )
compareParams(selectedIon, selectedIonReferenceParams)


#TESTING Si, defaultFetchSiAttrFromSmi(), fetchParentIon(), fetchScanInfo(), fetchSpectrumInfo()
smi = spectra[1]['smi']
si = maspy.new.core.Si(smi.id, smi.specfile)
maspy.new.import_export.defaultFetchSiAttrFromSmi(smi, si)
siAttrReference = {'basepeakI': 170753.84, 'basepeakMz': 251.1017004,'i': 725709.25, 'id': '5',
                   'iit': 59.999998658895,'isValid': None, 'msLevel': 2, 'mz': 410.681187454293,
                   'rt': 1.39947804, 'specfile': 'test', 'tic': 637463.25, 'charge': 2
                   }
for key, refValue in viewitems(siAttrReference):
    siAttrValue = getattr(si, key)
    try:
        print(round(siAttrValue, 6) == round(refValue, 6))
        round(siAttrValue, 6) == round(refValue, 6)
    except TypeError:
        print(siAttrValue == refValue)
        siAttrValue == refValue

parentIonAttr = maspy.new.import_export.fetchParentIon(smi)
parentIonAttrReference = {'i': 7.2570925e05, 'mz': 410.681187454293, 'charge': 2, 'precursorId': 'controllerType=0 controllerNumber=1 scan=4'}
for key, value in viewitems(parentIonAttr):
    try:
        round(value, 6) == round(parentIonAttrReference[key], 6)
        print(round(value, 6) == round(parentIonAttrReference[key], 6))
    except TypeError:
        value == parentIonAttrReference[key]
        print(value == parentIonAttrReference[key])

scanAttr = maspy.new.import_export.fetchScanInfo(smi)
scanAttrReference = {'iit': 59.999998658895, 'rt': 0.023324634*60}
for key, value in viewitems(scanAttr):
    round(value, 6) == round(scanAttrReference[key], 6)
    print(round(value, 6) == round(scanAttrReference[key], 6))

specAttr = maspy.new.import_export.fetchSpectrumInfo(smi)
specAttrReference = {'basepeakI': 1.7075384e05, 'basepeakMz': 251.1017004, 'msLevel': 2, 'tic': 6.3746325e05}
for key, value in viewitems(specAttr):
    round(value, 6) == round(specAttrReference[key], 6)
    print(round(value, 6) == round(specAttrReference[key], 6))


#TESTING the results of mzml.extractBinaries()
sai = spectra[0]['sai']
smi = spectra[0]['smi']
saiArrayInfoReference = {'i': {'dataProcessingRef': None,
                               'params': [('MS:1000523', '', None),
                                          ('MS:1000574', '', None),
                                          ('MS:1000515', '', 'MS:1000131')
                                          ]
                                },
                         'mz': {'dataProcessingRef': None,
                                'params': [('MS:1000523', '', None),
                                           ('MS:1000574', '', None),
                                           ('MS:1000514', '', 'MS:1000040')
                                           ]
                                }
                          }
for key, value in viewitems(sai.arrayInfo):
    value['dataProcessingRef'] == saiArrayInfoReference[key]['dataProcessingRef']
    print(value['dataProcessingRef'] == saiArrayInfoReference[key]['dataProcessingRef'])
    compareParams(value['params'] , saiArrayInfoReference[key]['params'])
    print(compareParams(value['params'] , saiArrayInfoReference[key]['params']))

for key, value in viewitems(sai.arrays):
    value.size == int(smi.attributes['defaultArrayLength'])
    print(value.size == int(smi.attributes['defaultArrayLength']))

#TODO: ciFromXml
#TODO: xmlGenScanList, xmlGenPrecursorList, xmlGenProductList, xmlGenBinaryDataArrayList, xmlSpectrumFromSmi, xmlChromatogramFromCi


# --- Testing import / load / save MsrunContainer --- #
specfilename = 'JD_06232014_sample1_A'
specfiles = ['JD_06232014_sample1_A']

mzmlPath = 'J:/lab/maspy_testfiles/' + specfilename + '.mzML'
msrunContainer = maspy.new.import_export.importMzml(mzmlPath)

msrunContainer.setPath('D:/maspy_test')
msrunContainer.save()

newMsrunContainer = maspy.new.core.MsrunContainer()
newMsrunContainer.addSpecfile(specfiles, 'D:/maspy_test')
newMsrunContainer.load()

# --- Testing writing an mzml file and reimportint it --- #
maspy.new.import_export.writeMzml(specfilename, newMsrunContainer, 'D:/maspy_test', spectrumIds=None, chromatogramIds=None)
newMzmlMsrunContainer = maspy.new.import_export.importMzml('D:/maspy_test/' + specfilename + '.mzML')


# --- Testing import sii (SpectrumIdentificationItem) and SiiContainer --- #
psmFolder = 'J:/lab/maspy_testfiles'
siiContainer = maspy.new.core.SiiContainer()
for specfile in specfiles:
    print('Importing PSMs from: ', specfile)
    currSearchFolder = os.path.join(psmFolder, specfile).replace('\\','/')
    percolatorTargetOutputLocation = os.path.join(currSearchFolder, 'target_percolator_output.tsv').replace('\\','/')
    maspy.new.import_export.importPsmResults(siiContainer, percolatorTargetOutputLocation, specfile, psmType='percolator', psmEngine='comet', qValue=0.01)
siiContainer.setPath('D:/maspy_test')
siiContainer.save()

newSiiContainer = maspy.new.core.SiiContainer()
newSiiContainer.addSpecfile(specfile, 'D:/maspy_test')
newSiiContainer.load()

newSiiContainer.addSiInfo(msrunContainer, attributes=['mz', 'rt', 'charge'])
newSiiContainer.calcMz()


# --- Testing import fi (FeatureItem) and SiContainer --- #
specfiles = ['JD_06232014_sample1_A']
folder = 'C:/Users/David/Dropbox/python/maspy_workspace/maspy_testfiles'
folder = 'J:/lab/maspy_testfiles'
fiContainer = maspy.new.core.FiContainer()
for specfile in specfiles:
    print('Importing Fi from: ', specfile)
    currSearchFolder = os.path.join(folder, specfile).replace('\\','/')
    filelocation = os.path.join(folder, specfile+'.featureXML').replace('\\','/')
    maspy.new.import_export.importPeptideFeatures(fiContainer, filelocation, specfile)
fiContainer.setPath('D:/maspy_test')
fiContainer.save()

newFiContainer = maspy.new.core.FiContainer()
newFiContainer.addSpecfile(specfile, 'D:/maspy_test')
newFiContainer.load()


items = [item for item in maspy.new.core._getListItems(siiContainer.container, sort='score', reverse=True, selector=lambda item: item.isValid)]
siItems = [_ for _ in newMsrunContainer.getItems(selector=lambda si: si.msLevel==2 and si.mz > 1000)]
siArrays = newMsrunContainer.getArrays(['rt', 'mz', 'i'], selector= lambda si: si.msLevel==2)

siiItems = [item for item in newSiiContainer.getItems(sort='score', reverse=True, selector=lambda item: item.isValid)]
siiArrays = newSiiContainer.getArrays(['rt', 'charge', 'mz'])

fiItems = [item for item in newFiContainer.getItems(sort='mz', reverse=False, selector=lambda item: True)]
fiArrays = newFiContainer.getArrays(['rt', 'mz'])

from matplotlib import pyplot as plt
fig, ax = plt.subplots(4, sharex=True, sharey=True)
ax[0].scatter(siArrays['rt'], siArrays['mz'], marker='.', alpha=0.1, color='grey')
ax[1].scatter(siiArrays['rt'], siiArrays['mz'], marker='.', alpha=0.1, color='red')
ax[2].scatter(fiArrays['rt'], fiArrays['mz'], marker='.', alpha=0.1, color='green')
ax[3].scatter(siArrays['rt'], siArrays['mz'], marker='.', alpha=0.5, color='grey')
ax[3].scatter(siiArrays['rt'], siiArrays['mz'], marker='.', alpha=0.5, color='red')
ax[3].scatter(fiArrays['rt'], fiArrays['mz'], marker='.', alpha=0.5, color='green')
plt.show()



# --- Testing importing a protein database --- #
