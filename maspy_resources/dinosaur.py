"""
This module is intended for usage of Dinosaur version 1.1.3 at
https://github.com/fickludd/dinosaur
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

def execute(filelocation, outdir, args=None,
            executable='Dinosaur-1.1.3.free.jar', xmx='2G'):
    """Executes the dinosaur tool on Windows operating systems.

    :param filelocation: mzml input file path
    :param outdir: path of the output directory
    :param args: list of dinosaur arguments, for details see the dinosaur help.
        Arguments with a has to be added as "argument=value".
    :param executable: must specify the complete file path of the msConvert.exe
        if its location is not in the ``PATH`` environment variable.
    :param xmx: java maximum heap size
    """

    procArgs = ['java', '-jar', '-Xmx'+xmx, executable]
    if args is not None:
        procArgs.extend(['--'+arg for arg in args])

    procArgs.extend(['--outDir='+outdir])
    procArgs.append(filelocation)

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
outdir = 'D:/maspy_test/dinosaur'
filelocation = 'D:/maspy_test/test/20160608_QexHF4_RSLCbeta_25ng_HeLa_short_1305_02.mzML'
dinosaurArguments = [()]
executable = 'D:/programs/dinosaur/Dinosaur-1.1.3.free.jar'
execute(filelocation, outdir, executable=executable,
        args=['concurrency=4'], xmx='16G')
"""

"""Dinosaur-1.1.3.free.jar, OPTIONS:
          PARAMETER   DEFAULT   DESCRIPTION
            advHelp   false     set to output adv param file help and quit
          advParams             path to adv param file
        concurrency   2         the number of assays to analyze in parallel
              force   false     ignore missing mzML params
          maxCharge   6         max searched ion charge
          minCharge   1         min searched ion charge
               mode   global    analysis mode: global or target. Global mode reports all isotope patterns, targeted only those matching targets.
               mzML   -         The shotgun MzML file to analyze
            nReport   10        number of random assay to export control figure for
             outDir             output directory (by default same as input mzML)
            outName             basename for output files (by default same as input mzML)
          profiling   false     set to enable CPU profiling
reportDeisoMzHeight   15.0      mz range in deisotoper reports
      reportHighRes   false     generate high-resolution plot trail when supported (for print)
         reportSeed   -1        seed to use for report assay selection (<0 means random)
      reportTargets   false     set to create a special report figure for each target
   targetPreference   rt        if multiple isotope patterns fit target, take the closest rt apex (rt) or the most intense (intensity)

            targets             path to isotope patterns target file (not used by default)
            verbose   false     increase details in output
        writeBinary   false     set to output binary MSFeatureProtocol file
         writeHills   false     set to output csv file with all hills assigned to isotope patterns
     writeMsInspect   false     set to output MsInspect feature csv file
       writeQuantML   false     set to output mzQuantML file
        zipQcFolder   false     set to zip the entire qc folder on algorithm completion

              PARAMETER   DEFAULT   DESCRIPTION
          averagineCorr   0.6       Required correlation with averagine
     averagineExplained   0.5       Required probability mass of averagine explained
     backgroundQuantile   0.0       (NOT IMPL) low abund quantile to remove
         chargePairCorr   0.6       required cosine correlation between intensity profiles of charge pairs
          chargePairPPM   7.0       Max ppm diff to allow when assigning charge pairs
         deisoBatchSize   500       Number of hill clusters to batch into non-parallel computation unit
              deisoCorr   0.6       required cosine correlation between intensity profiles of hills in feature
    deisoDecompMaxSeeds   100       max seeds checked for feature decomposition per round
           deisoOverlap   3         minimum overlap required when using the averagine-overlap model
            deisoSigmas   5.0       number of m/z std-devs from expected to tolerate in feature isotopes
      deisoValleyFactor   1.3       split feature if smallest surrounding local max is larger than valleyFactor * local min
          hillBatchSize   100       Number of consecutive spectra to batch into non-parallel computation unit
         hillMaxMissing   1         Max number of consecutive centroids that can be missing in a hill
          hillMinLength   3         Min number of centroid in hill
      hillMzGuessLength   3         number of recent mz values to use for live-guessing hill mz
             hillNBoots   150       n bootstrap calculations for mz determination
                hillPPM   8.0       Max ppm diff to allow in consecutive centroids for joining as a hill
         hillPeakFactor   2         require hill smoothed endpoints to be < PeakFactor * apexIntensity
hillPeakFactorMinLength   40        only use hillPeakFactor on hills at least this long
   hillSmoothMeanWindow   1         number of +-data points for the hill smooth mean filter (applied after median)
 hillSmoothMedianWindow   1         number of +-data points for the hill smooth median filter (applied before mean)
       hillValleyFactor   1.3       Split hill if smallest surrounding local max is larger than valleyFactor * local min
     intensityThreshold   0.0       (NOT IMPL) hard intensity threshold below which peaks are removed
         massCalcNBoots   150       n bootstrap calculations for mass calcuation
          massEstPoints   3         Number of most abundant peaks to use for centroid determination
            maxBootSize   300       max number of samples to use in each bootstrap iteration during mz determination
           maxIntensity   false     Use max intensity for centroided peaks instead of the sum
            noHillSplit   false     Don't split hills based on intensity profile
     reportTargetIntMax   2000.0    Max value for spectral heatmaps unless bigger real max
     subtractBackground   false     (NOT IMPL) Remove spectral background
"""
