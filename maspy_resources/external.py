from __future__ import print_function, division, unicode_literals
from future.utils import viewkeys, viewvalues, viewitems, listvalues, listitems

try: # python 2.7
    from itertools import izip as zip
except ImportError: # python 3.x series
    pass
################################################################################

import os
import subprocess
import sys

import maspy.auxiliary as aux

# --- MsConvert --- #
def runMsConvert(filelocation, outformat, outdir, params='',
                 executable='msConvert.exe'):
    """Execute the msConvert tool on Windows operating systems.

    :param filelocation: input file path
    :param outformat: output format, must be one of foloowing: "mzML, mzXML,
        mz5, mgf, text, ms1, cms1, ms2"
    :param params: str(), specify additional parameters and filters, for details
        see the msConvert.exe help.
    :param outdir: path of the output directory
    :param executable: must specify the complete file path of the msConvert.exe
        if its location is not in the ``PATH`` environment variable.
    """

    processList = ['"'+executable+'"', '"'+filelocation+'"', '--'+outformat,
                   params, '-o', '"'+outdir+'"'
                   ]
    command = ' '.join(processList)

    ## run it ##
    p = subprocess.Popen(command, shell=True, stderr=subprocess.PIPE)

    ## But do not wait till netstat finish, start displaying output immediately ##
    while True:
        out = p.stderr.read(1)
        if out == '' and p.poll() != None:
            break
        if out != '':
            sys.stdout.write(out)
            sys.stdout.flush()

"""
outdir = 'D:/maspy_test/pparse'
filelocation = 'D:/maspy_test/pparse/20160608_QexHF4_RSLCbeta_25ng_HeLa_short_1305_02.raw'
msConvertParams = '--zlib --64 --filter "peakPicking true [1,2]"'
runMsConvert(filelocation, 'mzML', outdir, msConvertParams,
             executable='D:/programs/ProteoWizard 3.0.7320/msconvert.exe')
"""

# --- pParse ---#
def generatePparseParams(rawfilepath, outputpath, isolationWindow, coElute=0):
    """Generates a string containing the parameters for a pParse parameter file
    but doesn't write any file yet.

    :param rawfilepath: location of the thermo ".raw" file
    :param outputpath: path to the output directory of pParse
    :param isolationWindow: MSn isolation window that was used for the
        aquisition of the specified thermo raw file
    :param coElute: 0 or 1, see "[Advanced Options]" below

    :returns: string containing pParse parameters

    .. note:
        # pParse.para params template
        # For help: mail to tuhuijun@ict.ac.cn
        # Time: 2014.12.08

        [Basic Options]
        datapath = C:\filedirectory\filename
        logfilepath = C:\filedirectory
        outputpath = C:\filedirectory

        [Advanced Options]
        co-elute = 1
        # 0, output single precursor for single scan;
        # 1, output all co-eluted precursors.
        input_format = raw
        # raw / ms1
        isolation_width = 1.6
        # 2 / 2.5 / 3 / 4
        mars_threshold = -0.5
        ipv_file = .\IPV.txt
        trainingset = EmptyPath

        [Internal Switches]
        output_mars_y = 0
        delete_msn = 0
        output_mgf = 1
        output_pf = 1
        debug_mode = 0
        check_activationcenter = 1
        output_all_mars_y = 0
        rewrite_files = 0
        export_unchecked_mono = 0
        cut_similiar_mono = 1
        mars_model = 4
        output_trainingdata = 0
    """
    output = str()
    #Basic options
    output = '\n'.join([output, ' = '.join(['datapath', rawfilepath])])
    output = '\n'.join([output, ' = '.join(['logfilepath', outputpath])])
    output = '\n'.join([output, ' = '.join(['outputpath', outputpath])])
    #Advanced options
    output = '\n'.join([output, ' = '.join(['co-elute', str(coElute)])])
    output = '\n'.join([output, ' = '.join(['input_format', 'raw'])])
    output = '\n'.join([output, ' = '.join(['isolation_width',
                                            str(isolationWindow)]
                                            )])
    output = '\n'.join([output, ' = '.join(['mars_threshold', '-0.5'])])
    output = '\n'.join([output, ' = '.join(['ipv_file', '.\IPV.txt'])])
    output = '\n'.join([output, ' = '.join(['trainingset', 'EmptyPath'])])
    #Internal Switches
    output = '\n'.join([output, ' = '.join(['output_mars_y', '0'])])
    output = '\n'.join([output, ' = '.join(['delete_msn', '0'])])
    output = '\n'.join([output, ' = '.join(['output_mgf', '1'])])
    output = '\n'.join([output, ' = '.join(['output_pf', '0'])])
    output = '\n'.join([output, ' = '.join(['debug_mode', '0'])])
    output = '\n'.join([output, ' = '.join(['check_activationcenter', '1'])])
    output = '\n'.join([output, ' = '.join(['output_all_mars_y', '0'])])
    output = '\n'.join([output, ' = '.join(['rewrite_files', '0'])])
    output = '\n'.join([output, ' = '.join(['export_unchecked_mono', '0'])])
    output = '\n'.join([output, ' = '.join(['cut_similiar_mono', '1'])])
    output = '\n'.join([output, ' = '.join(['mars_model', '4'])])
    output = '\n'.join([output, ' = '.join(['output_trainingdata', '0'])])

    return output


def preparePparse(rawfilepath, outputpath, isolationWindow, coElute=0):
    """Generate and write a pParse parameter file.

    :param rawfilepath: location of the thermo ".raw" file
    :param outputpath: path to the output directory of pParse
    :param isolationWindow: MSn isolation window that was used for the
        aquisition of the specified thermo raw file
    :param coElute:

    :returns: file path of the pParse parameter file
    """
    paramText = generatePparseParams(rawfilepath, outputpath, isolationWindow,
                                     coElute=coElute)
    paramPath = aux.joinpath(outputpath, 'pParse.para')
    with open(paramPath, 'wb') as openfile:
        openfile.write(paramText)
    return paramPath


def executePparse(paramPath, executable='pParse.exe'):
    """Execute pParse with the specified parameter file.

    :param paramPath: location of the pParse parameter file
    :param executable: must specify the complete file path of the pParse.exe
        if its location is not in the ``PATH`` environment variable.

    :returns: :func:`subprocess.Popen` return code, 0 if pParse was executed
        successful
    """
    processList = ['"'+executable+'"', '"'+paramPath+'"']
    command = ' '.join(processList)

    ## run it ##
    p = subprocess.Popen(command, shell=True, stderr=subprocess.PIPE)

    ## But do not wait till netstat finish, start displaying output immediately ##
    while True:
        out = p.stderr.read(1)
        if out == '' and p.poll() != None:
            break
        if out != '':
            sys.stdout.write(out)
            sys.stdout.flush()

    return p.returncode


def cleanUpPparse(outputpath, filename, mgf=False):
    """Delete temporary files generated by pparse, including the filetypes
    ".csv", ".ms1", ".ms2", ".xtract",  the files "pParsePlusLog.txt" and
    "pParse.para" and optionally also the ".mgf" file generated by pParse.

    .. warning:
        When the parameter "mgf" is set to "True" all files ending with ".mgf"
        and containing the specified "filename" are deleted. This could
        potentially also affect MGF files not generated by pParse.

    :param outputpath: path to the output directory of pParse
    :param filename: specfile name
    :param mgf: bool, if True the ".mgf" file generated by pParse is also
        removed
    """
    extensions = ['csv', 'ms1', 'ms2', 'xtract']
    basename, fileext = os.path.splitext(os.path.basename(filename))
    additionalFiles = [aux.joinpath(outputpath, 'pParsePlusLog.txt'),
                       aux.joinpath(outputpath, 'pParse.para'),
                       ]

    for ext in extensions:
        filepath = aux.joinpath(outputpath, '.'.join([basename, ext]))
        if os.path.isfile(filepath):
            print('Removing file: ', filepath)
            os.remove(filepath)
    for filepath in additionalFiles:
        if os.path.isfile(filepath):
            print('Removing file: ', filepath)
            os.remove(filepath)
    if mgf:
        for _filename in os.listdir(outputpath):
            _basename, _fileext = os.path.splitext(_filename)
            if _fileext.lower() != '.mgf':
                continue
            if _basename.find(basename) != -1 and _basename != basename:
                filepath = aux.joinpath(outputpath, _filename)
                print('Removing file: ', filepath)
                os.remove(filepath)
