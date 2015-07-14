from __future__ import print_function
import numpy
import os
from collections import defaultdict as ddict
from matplotlib import pyplot as plt

try:
    import seaborn
except:
    pass

import pyms.core as pyms
import pyms.auxiliary as aux

import timeit
globalStart = timeit.time.time()

specfiles = list()
specfiles.append('JD_06232014_sample1_A')
specfiles.append('JD_06232014_sample2_B')
specfiles.append('JD_06232014_sample3_C')

specFileDir = 'D:/massData/specFiles/pymzml'

# Import spectrum files
siContainer = pyms.importSpecfiles(specfiles, specFileDir)

# Generate QC data for plotting
targetSpecfiles = specfiles #[specfiles[0], specfiles[1], specfiles[2]]
specfileNumber = len(targetSpecfiles)

arrays = siContainer.getArrays(['rt', 'tic', 'iit', 'msLevel', 'msnIdList'])
qcReport = ddict(dict)

for specfileCounter, specfile in enumerate(targetSpecfiles):
    specfileMask = (arrays['specfile'] == specfile)
    rtMask = (arrays['rt'] > 0)
    ms1Mask = (arrays['msLevel'] == 1) & rtMask & specfileMask
    ms2Mask = (arrays['msLevel'] == 2) & rtMask & specfileMask

    qcReport[specfile]['reportPos'] = specfileCounter
    qcReport[specfile]['rtMs1'] = arrays['rt'][ms1Mask]
    qcReport[specfile]['rtMs2'] = arrays['rt'][ms2Mask]

    qcReport[specfile]['numMs1'] = ms1Mask.sum()
    qcReport[specfile]['numMs2'] = ms2Mask.sum()

    qcReport[specfile]['ms2PerMs1'] = numpy.array([len(x) for x in arrays['msnIdList'][ms1Mask]])

    qcReport[specfile]['ticListMs1'] = arrays['tic'][ms1Mask]
    qcReport[specfile]['ticListMs2'] = arrays['tic'][ms2Mask]

    qcReport[specfile]['iitListMs1'] = arrays['iit'][ms1Mask]
    qcReport[specfile]['iitListMs2'] = arrays['iit'][ms2Mask]


# Generate plot
plotTitle = 'Comparison of '+str(specfileNumber)+' spectrum files:'
colorList = ['#73d216','#3465a4','#f57900','#75507b','#cc0000', 'grey']
generalAlpha = 0.75

plt.clf()
f, axarr = plt.subplots(3, 2,figsize=(12,15))
f.suptitle(plotTitle, fontsize=14, weight = 'bold')
f.subplots_adjust(hspace=0.4)


# Histogram - Number of MS2 scans per MS1 scan
dataList = list()
for specfile in targetSpecfiles:
    report = qcReport[specfile]
    ms2PerMs1 = report['ms2PerMs1']
    dataList.append(ms2PerMs1)
maxValue = max([max(y) for y in dataList])
binList = [x-0.5 for x in range(0,maxValue+2)]
tickList = range(0,maxValue+1)

axarr[0, 0].set_title('Number of MS2 scans per MS1 full scan', fontsize=13)
axarr[0, 0].hist(dataList, bins=binList, alpha=generalAlpha, color=colorList[0:specfileNumber], normed=False)
axarr[0, 0].set_xlabel('MS2 scans after MS1 full scan')
axarr[0, 0].set_xticks(tickList)
axarr[0, 0].set_xlim(-0.5,maxValue+0.5)
axarr[0, 0].set_ylabel('MS1 counts')
axarr[0, 0].legend(targetSpecfiles, prop={'size':7}, bbox_to_anchor=(-0.2, 1.04, 2.6, .102), loc=3,ncol=2, mode="expand", borderaxespad=1.5)


# Number of MS1 and MS2 scans
xTickLabelList = list()
for specfile in targetSpecfiles:
    report = qcReport[specfile]
    xAxis = (report['reportPos'],report['reportPos']+0.5)
    yAxis = (report['numMs1'],report['numMs2'])
    axarr[0, 1].bar(xAxis, yAxis, 0.5, alpha=generalAlpha, color=colorList[report['reportPos']])
    xTickLabelList.append('# MS1')
    xTickLabelList.append('# MS2')

axarr[0, 1].set_title('Number of MS1 and MS2 scans ', fontsize=13)
axarr[0, 1].set_xticklabels(xTickLabelList)
axarr[0, 1].set_xticks([x/2.0 +0.25 for x in range(0,specfileNumber*2)])
axarr[0, 1].set_ylabel('Counts')


# Histogram - TotalIonCurrent (log10) MS1
dataList = list()
for specfile in targetSpecfiles:
    report = qcReport[specfile]
    dataList.append([numpy.log10(x) for x in report['ticListMs1']])
minBin = int(min([min(y) for y in dataList]))
maxBin = int(max([max(y) for y in dataList]))+1
binList = [x/5.0 for x in range(minBin*5,maxBin*5)]
tickList = [x/1.0 for x in range(minBin*1,maxBin*1)]

axarr[1, 0].set_title('Total Ion Current distribution of all MS1 scans', fontsize=13)
axarr[1, 0].hist(dataList, bins=binList, alpha=generalAlpha, color=colorList[0:specfileNumber], normed=False)
axarr[1, 0].set_xlabel('TIC (log10)')
axarr[1, 0].set_xticks(tickList)
axarr[1, 0].set_xticklabels(['E'+str(x) for x in tickList])
axarr[1, 0].set_ylabel('MS1 counts')


# Histogram - TotalIonCurrent (log10) MS2
dataList = list()
for specfile in targetSpecfiles:
    report = qcReport[specfile]
    dataList.append([numpy.log10(x) for x in report['ticListMs2']])
minBin = int(min([min(y) for y in dataList]))
maxBin = int(max([max(y) for y in dataList]))+1
binList = [x/5.0 for x in range(minBin*5,maxBin*5)]
tickList = [x/1.0 for x in range(minBin*1,maxBin*1)]

axarr[1, 1].set_title('Total Ion Current distribution of all MS2 scans', fontsize=13)
axarr[1, 1].hist(dataList,bins=binList, alpha=generalAlpha, color=colorList[0:specfileNumber], normed=False)
axarr[1, 1].set_xlabel('TIC (log10)')
axarr[1, 1].set_xticks(tickList)
axarr[1, 1].set_xticklabels(['E'+str(x) for x in tickList])
axarr[1, 1].set_ylabel('MS2 counts')


# Histogram - IonInjectionTime MS1
dataList = list()
for specfile in targetSpecfiles:
    report = qcReport[specfile]
    dataList.append(report['iitListMs1'])

binNumber = 15
binStep = int( numpy.ceil( max([max(y) for y in dataList]) / binNumber ) )
binList = [x*binStep for x in range(0,binNumber+1)]
tickNumber = 5
tickStep = int( numpy.ceil( max([max(y) for y in dataList]) / tickNumber ) )
tickList = [x*tickStep for x in range(0,tickNumber+1)]

axarr[2, 0].set_title('Ion Injection Time distribution of all MS1 scans', fontsize=13)
axarr[2, 0].hist(dataList,bins=binList, alpha=generalAlpha, color=colorList[0:specfileNumber], normed=False)
axarr[2, 0].set_xlabel('Ion Injection Time (milliseconds)')
axarr[2, 0].set_xticks(tickList)
axarr[2, 0].set_xticklabels([str(x)+' ms' for x in tickList])
axarr[2, 0].set_ylabel('MS1 counts')


# Histogram - IonInjectionTime MS2
dataList = list()
for specfile in targetSpecfiles:
    report = qcReport[specfile]
    dataList.append(report['iitListMs2'])

binNumber = 15
binStep = int( numpy.ceil( max([max(y) for y in dataList]) / binNumber ) )
binList = [x*binStep for x in range(0,binNumber+1)]
tickNumber = 5
tickStep = int( numpy.ceil( max([max(y) for y in dataList]) / tickNumber ) )
tickList = [x*tickStep for x in range(0,tickNumber+1)]

axarr[2, 1].set_title('Ion Injection Time distribution of all MS2 scans', fontsize=13)
axarr[2, 1].hist(dataList,bins=binList, alpha=generalAlpha, color=colorList[0:specfileNumber], normed=False)
axarr[2, 1].set_xlabel('Ion Injection Time (milliseconds)')
axarr[2, 1].set_xticks(tickList)
axarr[2, 1].set_xticklabels([str(x)+' ms' for x in tickList])
axarr[2, 1].set_ylabel('MS2 counts')

globalEnd = timeit.time.time()
print(globalEnd - globalStart, 'processing time in seconds.')

plt.show()
plt.clf()




for specfile in specfiles:
    plt.scatter(qcReport[specfile]['rtMs1'], qcReport[specfile]['iitListMs1'], color='grey', alpha=0.2)
    plt.show()




""" This is something different as the rawMeat like report
for specfile in specfiles[0:2]:
    values = qcReport[specfile]['iitListMs2']
    binStep = 20
    bins = [x * binStep for x in range(0, int(max(values)/binStep) +2)]
    plt.hist(values, bins=bins, alpha=0.3, label=specfile)
plt.legend()
plt.show()


for specfile in specfiles[0:2]:
    medians, lowerPercentile, upperPercentile, indVar = recalFU.returnRunningMedian(qcReport[specfile]['rtMs1'],
                                                                                    qcReport[specfile]['ms2PerMs1'],
                                                                                    windowSize=50
                                                                                    )
    plt.scatter(qcReport[specfile]['rtMs1'], qcReport[specfile]['ms2PerMs1'], alpha=0.1, color='grey')
    plt.plot(indVar, medians, alpha=0.5)
    plt.show()


    medians, lowerPercentile, upperPercentile, indVar = recalFU.returnRunningMedian(qcReport[specfile]['rtMs1'],
                                                                                    qcReport[specfile]['iitListMs1'],
                                                                                    windowSize=200
                                                                                    )
    plt.plot(indVar, medians)
    plt.plot(indVar, lowerPercentile, color='grey')
    plt.ylim(0,max(medians)+5)
    plt.show()

"""



