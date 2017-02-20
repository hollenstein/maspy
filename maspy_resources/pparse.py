"""
This module is intended for usage of pParse version 2.0 from 20160203 at
http://pfind.ict.ac.cn/software/pParse/index.html

For the download see the section "Supplemental Files" on the webpage.
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

import os
import subprocess
import sys

import maspy.auxiliary as aux

def generateParams(rawfilepath, outputpath, isolationWindow, coElute):
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


def writeParams(rawfilepath, outputpath, isolationWindow, coElute=0):
    """Generate and write a pParse parameter file.

    :param rawfilepath: location of the thermo ".raw" file
    :param outputpath: path to the output directory of pParse
    :param isolationWindow: MSn isolation window that was used for the
        aquisition of the specified thermo raw file
    :param coElute:

    :returns: file path of the pParse parameter file
    """
    paramText = generateParams(rawfilepath, outputpath, isolationWindow,
                               coElute)
    filename, fileext = os.path.splitext(os.path.basename(rawfilepath))
    paramPath = aux.joinpath(outputpath, filename+'.pparse.para')
    with open(paramPath, 'wb') as openfile:
        openfile.write(paramText)
    return paramPath


def execute(paramPath, executable='pParse.exe'):
    """Execute pParse with the specified parameter file.

    :param paramPath: location of the pParse parameter file
    :param executable: must specify the complete file path of the pParse.exe
        if its location is not in the ``PATH`` environment variable.

    :returns: :func:`subprocess.Popen` return code, 0 if pParse was executed
        successful
    """
    procArgs = [executable, paramPath]

    ## run it ##
    proc = subprocess.Popen(procArgs, stderr=subprocess.PIPE)

    ## But do not wait till netstat finish, start displaying output immediately ##
    while True:
        out = proc.stderr.read(1)
        if out == '' and proc.poll() != None:
            break
        if out != '':
            sys.stdout.write(out)
            sys.stdout.flush()

    return proc.returncode


def cleanUpPparse(outputpath, rawfilename, mgf=False):
    """Delete temporary files generated by pparse, including the filetypes
    ".csv", ".ms1", ".ms2", ".xtract",  the files "pParsePlusLog.txt" and
    "pParse.para" and optionally also the ".mgf" file generated by pParse.

    .. warning:
        When the parameter "mgf" is set to "True" all files ending with ".mgf"
        and containing the specified "filename" are deleted. This could
        potentially also affect MGF files not generated by pParse.

    :param outputpath: path to the output directory of pParse
    :param rawfilename: filename of the thermo ".raw" file
    :param mgf: bool, if True the ".mgf" file generated by pParse is also
        removed
    """
    extensions = ['csv', 'ms1', 'ms2', 'xtract']
    filename, fileext = os.path.splitext(os.path.basename(rawfilename))
    additionalFiles = [aux.joinpath(outputpath, 'pParsePlusLog.txt'),
                       aux.joinpath(outputpath, filename+'.pparse.para'),
                       ]

    for ext in extensions:
        filepath = aux.joinpath(outputpath, '.'.join([filename, ext]))
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
