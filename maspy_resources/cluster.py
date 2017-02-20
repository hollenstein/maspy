"""
This module is intended for usage of spectra-cluster-cli version 1.0.2 at
https://github.com/spectra-cluster/spectra-cluster-cli

To install download and extract the release 1.0.2 BETA from
https://github.com/spectra-cluster/spectra-cluster-cli/releases
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

import subprocess
import sys

import maspy.auxiliary as aux

def execute(filelocation, outpath, executable, args=None, switchArgs=None):
    """Executes the dinosaur tool on Windows operating systems.

    :param filelocation: either a single mgf file path or a list of file paths.
    :param outpath: path of the output file, file must not exist
    :param executable: must specify the complete file path of the
        spectra-cluster-cli.jar file, supported version is 1.0.2 BETA.
    :param args: list of arguments containing a value, for details see the
        spectra-cluster-cli help. Arguments should be added as tuples or a list.
        For example: [('precursor_tolerance', '0.5'), ('rounds', '3')]
    :param switchArgs: list of arguments not containing a value, for details see
        the spectra-cluster-cli help. Arguments should be added as strings.
        For example: ['fast_mode', 'keep_binary_files']
    """
    procArgs = ['java', '-jar', executable]
    procArgs.extend(['-output_path', outpath])
    if args is not None:
        for arg in args:
            procArgs.extend(['-'+arg[0], arg[1]])
    if switchArgs is not None:
        procArgs.extend(['-'+arg for arg in switchArgs])

    procArgs.extend(aux.toList(filelocation))

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

"""
execute(['D:/massData/cluster_test/JD_06232014_sample1_A.mgf',
         'D:/massData/cluster_test/JD_06232014_sample1_B.mgf',
         'D:/massData/cluster_test/JD_06232014_sample1_C.mgf'
         ], 'D:/massData/cluster_test/cluster_03.cluster',
        'D:/programs/spectra_cluster_cli_1_0_2/spectra-cluster-cli-1.0.2.jar',
        args=None, switchArgs=None
        )
"""


""" spectra-cluster-cli help
usage: Spectra Cluster - Clusterer [-binary_directory <arg>] [-fast_mode]
       [-fragment_tolerance] [-help] [-keep_binary_files]
       [-major_peak_jobs <arg>] [-output_path <arg>] [-precursor_tolerance
       <arg>] [-remove_reporters <QUANTITATION TYPE>]
       [-reuse_binary_files] [-rounds <arg>] [-threshold_end <arg>]
       [-threshold_start <arg>] [-verbose] [-x_disable_mgf_comments]
       [-x_learn_cdf <output filename>] [-x_load_cdf <CDF filename>]
       [-x_min_comparisons <arg>] [-x_n_prefiltered_peaks <number peaks>]
Clusters the spectra found in an MGF file and writes the results in a
text-based file.
 -binary_directory <arg>                 path to the directory to
                                         (temporarily) store the binary
                                         files. By default a temporary
                                         directory is being created
 -fast_mode                              if this option is set the 'fast
                                         mode' is enabled. In this mode,
                                         the radical peak filtering used
                                         for the comparison function is
                                         already applied during spectrum
                                         conversion. Thereby, the
                                         clustering and consensus spectrum
                                         quality is slightly decreased but
                                         speed increases 2-3 fold.
 -fragment_tolerance                     fragment ion tolerance in m/z to
                                         use for fragment peak matching
 -help                                   print this message.
 -keep_binary_files                      if this options is set, the
                                         binary files are not deleted
                                         after clustering.
 -major_peak_jobs <arg>                  number of threads to use for
                                         major peak clustering.
 -output_path <arg>                      path to the outputfile.
                                         Outputfile must not exist.
 -precursor_tolerance <arg>              precursor tolerance (clustering
                                         window size) in m/z used during
                                         matching.
 -remove_reporters <QUANTITATION TYPE>   remove reporter ion peaks in
                                         quantitation experiments.
                                         Possible QUANTIATION TYPES are
                                         'ITRAQ', 'TMT' and 'ALL' ('TMT'
                                         and 'ITRAQ' peaks are removed.
 -reuse_binary_files                     if this option is set, the binary
                                         files found in the binary file
                                         directory will be used for
                                         clustering.
 -rounds <arg>                           number of clustering rounds to
                                         use.
 -threshold_end <arg>                    (lowest) final clustering
                                         threshold
 -threshold_start <arg>                  (highest) starting threshold
 -verbose                                if set additional status
                                         information is printed.
 -x_disable_mgf_comments                 (Advanced option) If set, MGF
                                         comment strings are NOT
                                         supported. This will increase
                                         performance but only works for
                                         MGF files that do not contain any
                                         comments
 -x_learn_cdf <output filename>          (Experimental option) Learn the
                                         used cumulative distribution
                                         function directly from the
                                         processed data. This is only
                                         recommended for high-resolution
                                         data. The result will be written
                                         to the defined file.
 -x_load_cdf <CDF filename>              (Experimental option) Loads the
                                         cumulative distribution function
                                         to use from the specified file.
                                         These files can be created using
                                         the x_learn_cdf parameter
 -x_min_comparisons <arg>                (Experimental option) Sets the
                                         minimum number of comparisons
                                         used to calculate the probability
                                         that incorrect spectra are
                                         clustered.
 -x_n_prefiltered_peaks <number peaks>   (Experimental option) Set the
                                         number of highest peaks that are
                                         kept per spectrum during loading.
"""
