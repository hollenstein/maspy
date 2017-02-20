"""
This module is intended for usage of RawConverter version 1.1.0.18 at
http://fields.scripps.edu/rawconv/
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


def execute(filelocation, outformat, outdir, log=False,
            executable='RawConverter.exe'):
    """Execute the msConvert tool on Windows operating systems.

    :param filelocation: input file path
    :param outformat: output format, must be one of the following: ms1, ms2, ms3, mgf
    :param outdir: path of the output directory
    :param log: #TODO
    :param executable: must specify the complete file path of the RawConverter.exe
        if its location is not in the ``PATH`` environment variable.

    .. note:
        Specifying the complete path to the executable is probably always
        necessary because RawConverter looks for the file "AveragineTable.txt"
        in the working directory.
    """
    assert outformat in ['ms1', 'ms2', 'ms3', 'mgf']

    args = [executable, filelocation, '--'+outformat, '--out_folder', outdir,
            '--select_mono_prec']

    ## run it ##
    proc = subprocess.Popen(args, cwd=os.path.dirname(executable),
                            stderr=subprocess.PIPE)

    ## But do not wait till netstat finish, start displaying output immediately ##
    while True:
        out = proc.stderr.read(1)
        if out == '' and proc.poll() != None:
            break
        if out != '':
            sys.stdout.write(out)
            sys.stdout.flush()


"""
Usage:
        RawXtract <input_file> [options]        run the commandline version.
        Options:
                --ms1   output MS1 file.
                --ms2   output MS2 file.
                --ms3   output MS3 file.
                --mgf   output MGF file.
                --log   output log file.
                --out_folder    output folder.
                --select_mono_prec      select the monoisotopic m/z values of precursors in DDA data.
                --predict_precursors    predict the precursors for DIA data.
"""

