from __future__ import print_function

import os
import subprocess
import time

# --- external non PSMatching tools (openMS, msConvert) --- #
def startMsConvert(sourceSpecFilePath, outputFileFormat, outputSpecFilePath, **kwargs):
    msConvertPath = kwargs.get('msConvertPath', 'msconvert')
    # startMsConvert, return 0 (success) or 1 (failed)
    if outputFileFormat not in ['mzML','mgf']:
        print('Wrong output file format specified: '+outputFileFormat+' (has to be \'mzML\' or \'mgf\')')
        return 1

    # add filters which were specified with additional arguments #
    filterDict = dict()
    filterList = list()

    filterDict['turbocharger'] = kwargs.get('turbocharger', False)
    if filterDict['turbocharger'] != False:
        filterList.append('\"' + 'turbocharger minCharge=2 maxCharge=6 precursorsBefore=2 precursorsAfter=2 halfIsoWidth=1.5' + '\"')

    filterDict['ms2Deisotope'] = kwargs.get('ms2Deisotope', False) #should be something like '0.005Da'
    if filterDict['ms2Deisotope'] != False:
        filterList.append('\"' + 'MS2Deisotope true ' + filterDict['ms2Deisotope'] + '\"')

    filterDict['peakPicking'] = kwargs.get('peakPicking', False) #should be something like '1' or '1,2'
    if filterDict['peakPicking'] != False:
        filterList.append('\"' + 'peakPicking [' + filterDict['peakPicking'] + ']\"')

    filterDict['peakPickingVendor'] = kwargs.get('peakPickingVendor', False) #should be something like '1' or '1,2'
    if filterDict['peakPickingVendor'] != False:
        filterList.append('\"' + 'peakPicking true [' + filterDict['peakPickingVendor'] + ']\"')

    processList = list()
    processList.append( msConvertPath )
    processList.append( sourceSpecFilePath )
    processList.append( '--'+outputFileFormat )

    if outputFileFormat == 'mzML':
        processList.append( '--zlib' )

    # add filter parameters to the processList if any filters where specified #
    for filterEntry in filterList:
        processList.extend(('--filter', filterEntry))

    processList.extend(('-o', outputSpecFilePath))
    print(processList)

    # execute the proceeList code to start msConvert #
    p = subprocess.Popen(processList, bufsize=4096 , stdout = subprocess.PIPE)
    stdout, stderr = p.communicate()
    if p.returncode == 0:# is 0 if success
        print('msConvert: succesfully executed')
        return 0
    else:
        print('msConvert: execution failed')
        print(stdout, stderr)
        return 1


def createOpenMsTrfFile(mzMLFileLocationList,trfFileLocation):
    with open(trfFileLocation, 'wb') as openFile:
        openFile.write('<?xml version="1.0" encoding="ISO-8859-1"?>'+'\n')
        openFile.write('<PARAMETERS version="1.3" xsi:noNamespaceSchemaLocation="http://open-ms.sourceforge.net/schemas/Param_1_3.xsd" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">'+'\n')
        ## Write NODE1 input List ##
        openFile.write('<NODE name="1" description="">'+'\n')
        openFile.write('<ITEMLIST name="url_list" type="string" description="">'+'\n')
        for mzMLFileLocation in mzMLFileLocationList:
            openFile.write('<LISTITEM value="file:///'+mzMLFileLocation+'"/>'+'\n')
        openFile.write('</ITEMLIST>'+'\n')
        openFile.write('</NODE>'+'\n')
        openFile.write('</PARAMETERS>'+'\n')


def runOpenMsWorkflow(workFlowPath,outputDir,fileLocationList, **kwargs):
    numJobs     = kwargs.get('numJobs', '1')
    numThreads  = kwargs.get('numThreads', '6')
    trfFilePath = kwargs.get('trfFilePath', os.path.join(outputDir,time.strftime('TEMP_%Y%m%d%H%M%S.trf') ).replace('\\','/') )
    createOpenMsTrfFile(fileLocationList,trfFilePath)

    processList = ['ExecutePipeline','-in',workFlowPath,'-out_dir',outputDir,'-resource_file',trfFilePath,'-num_jobs',str(numJobs)]#,'-threads',str(numThreads)]
    p = subprocess.Popen(processList, bufsize=4096 , stdout = subprocess.PIPE)
    stdout, stderr = p.communicate()
    if p.returncode == 0:
        print('Toppas workflow successfully executed: ' + workFlowPath)
    else:
        print('Toppas workflow aborted')
        print(stdout, stderr)

    os.remove(trfFilePath)
