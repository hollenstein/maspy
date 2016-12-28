# This is an outdated version that uses the module maspy_resources.external
# 07.09.2016
from __future__ import print_function, division, unicode_literals
from future.utils import viewkeys, viewvalues, viewitems, listvalues, listitems

try: # python 2.7
    from itertools import izip as zip
except ImportError: # python 3.x series
    pass
###############################################################################


import argparse
import glob
import os
import subprocess
import sys
sys.path.append('D:/Dropbox/python/maspy')
sys.path.append('C:/Users/David/Dropbox/python/maspy')

import pyteomics.mgf

import maspy.auxiliary as aux
import maspy.core
import maspy.reader
import maspy.writer
import maspy.xml

import maspy_resources.external as EXTERNAL


# --- Global parameters --- #
#pparsePath = 'pParse.exe'
#msconvertPath = 'msconvert.exe'
pparsePath = 'D:/shared/pParse_2.0_2016.02.03/pParse.exe'
msconvertPath = 'D:/programs/ProteoWizard 3.0.7320/msconvert.exe'
msconvertParams = '--zlib --filter "peakPicking true [1,2]"'


def extractIsolationWindows(msrunContainer):
    """Try to extract the isolation windows for all MS2 scans present in the
    msrunContainer

    :param msrunContainer: :class:`maspy.core.MsrunContaine`

    :returns: a list of extracted isolation window values

    """
    isolationWindows = set()
    for si in msrunContainer.getItems(selector=lambda si: si.msLevel == 2):
        smi = msrunContainer.smic[si.specfile][si.id]
        for prec in smi.precursorList:
            lowerLimit = maspy.xml.findParam(prec.isolationWindow,
                                             'MS:1000828'
                                             )[1]
            upperLimit = maspy.xml.findParam(prec.isolationWindow,
                                             'MS:1000829'
                                             )[1]
            isolationWindow = round(float(lowerLimit) + float(upperLimit), 3)
            isolationWindows.add(isolationWindow)
    return list(isolationWindows)


def pParseToMzml(rawfilepath, mzmlfilepath, pParseExePath, outputpath,
                 isolationWindow=None, removeMgf=True):
    """Execute pParse, import the mzML into MasPy, read the MGF file generated
    by pParse and replace precursor information in the msrunContainer. Export
    the msrunContainer to generate a new mzML. Remove files generated by pParse.

    :param rawfilepath: location of the thermo ".raw" file
    :param mzmlfilepath: location of the ".mzML" specfile, generated from the
        raw file specified in "rawfilepath".
    :param pParseExePath: location of the pParse executable
    :param outputpath: path to the output directory of pParse
    :param isolationWindow: MSn isolation window that was used for the
        aquisition of the specified thermo raw file
    :param removeMgf: bool, True if the ".mgf" file generated by pParse should
        be removed after the precursor information has been written to the mzML
        file.
    """
    #Read mzml file
    print('\nReading mzml file')
    print('---------------------------------')
    specfile = os.path.basename(mzmlfilepath)
    msrunContainer = maspy.reader.importMzml(mzmlfilepath,
                                             maspy.core.MsrunContainer()
                                             )
    if isolationWindow is None:
        windows = extractIsolationWindows(msrunContainer)
        if len(windows) > 1:
            raise Exception('Multiple different isolation windows present in' +
                            'mzML file: '+mzmlfilepath )
        else:
            isolationWindow = windows[0]
            print('MS2 isolation window = ', isolationWindow)

    #Execute pParse
    print('\nExecuting pParse')
    print('---------------------------------')
    paramFilepath = EXTERNAL.preparePparse(rawfilepath, outputpath,
                                           isolationWindow
                                           )
    pParseReturnCode = EXTERNAL.executePparse(paramFilepath,
                                              executable=pParseExePath
                                              )
    if pParseReturnCode == 0:# is 0 if success
        print('pParse: succesfully executed')
    else:
        raise Exception('pParse: execution failed')

    #Read mgf file
    print('\nReading mgf file')
    print('---------------------------------')
    basename, fileext = os.path.splitext(os.path.basename(rawfilepath))
    for _filename in os.listdir(outputpath):
        if _filename.find(basename) != -1 and _filename.endswith('.mgf'):
            filename = _filename
            break
    mgfFilepath = aux.joinpath(outputpath, filename)
    mgfRead = pyteomics.mgf.read(mgfFilepath, convert_arrays=1,
                                 read_charges=False
                                 )

    #Transfer precursor information from mgf to mzml
    print('\nTransfering precursor information')
    print('---------------------------------')
    for spectrum in mgfRead:
        title = os.path.basename(spectrum['params']['title'])
        specfile = title.split('.')[0]
        scanNr = title.split('.')[1]
        obsMz = str(spectrum['params']['pepmass'][0])
        charge = str(int(spectrum['params']['charge'][0]))

        # - replace the precursor m/z and charge value - #
        #Note: what happens in case of no charge? How is no charge handled by
        #   pyteomics? (eg 0 or None or a default charge value)
        smi = msrunContainer.smic[specfile][scanNr]
        if len(smi.precursorList) > 1:
            text = ' '.join(['Multiple precursorList entries in', str(smi)])
            raise Exception(text)
        if len(smi.precursorList[0].selectedIonList) > 1:
            text = ' '.join(['Multiple selected precursor ions in', str(smi)])
            raise Exception(text)
        selectedIonList = [(('MS:1000744', obsMz, 'MS:1000040'),
                            ('MS:1000041', charge, None)
                            )]
        smi.precursorList[0].selectedIonList = selectedIonList

    #Generate a new mzML file
    print('\nWriting new mzML file')
    print('---------------------------------')
    maspy.writer.writeMzml(specfile, msrunContainer, outputpath)

    #Cleaning up unnecessary pParse files
    print('\nCleaning up pParse files')
    print('---------------------------------')
    EXTERNAL.cleanUpPparse(outputpath, specfile, mgf=removeMgf)
    #Remove .ms1, .ms2, .xtract file fram raw file folder
    print('Done!')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Short sample app')
    parser.add_argument('-r', '--rawfile', action="store",
                        default=False, help='TODO')
    parser.add_argument('-i', '--isolationWindow', action="store",
                        default=None, help='TODO')

    arguments = parser.parse_args()
    currScriptPath = os.path.dirname(sys.argv[0])
    isolationWindow = arguments.isolationWindow

    if arguments.rawfile:
        paths = glob.glob(arguments.rawfile)
    else:
        paths = glob.glob(aux.joinpath(currScriptPath, '*.raw'))

    for rawfilepath in paths:
        if not rawfilepath.lower().endswith('.raw'):
            continue

        print()
        print('---------------------------------------------------------------')
        print('Processing ', rawfilepath)
        print('---------------------------------------------------------------')

        mzmlfilepath = rawfilepath.replace('.raw', '.mzML')
        outputpath = os.path.dirname(rawfilepath)

        # - use msConvert to generate a new mzML file from the raw file - #
        EXTERNAL.runMsConvert(rawfilepath, 'mzML', outputpath, msconvertParams,
                              executable=msconvertPath
                              )

        # - Import precursor information, update msrunContainer and write mzML - #
        pParseToMzml(rawfilepath, mzmlfilepath, pparsePath, outputpath,
                     isolationWindow
                     )

        # - Generate mgf and mzXML files - #
        EXTERNAL.runMsConvert(mzmlfilepath, 'mgf', outputpath,
                              executable=msconvertPath
                              )
        EXTERNAL.runMsConvert(mzmlfilepath, 'mzXML', outputpath, '--zlib --32',
                              executable=msconvertPath
                              )
        print('---------------------------------------------------------------')
    raw_input('Workflow complete, press any key to quit.')