"""
This module is intended for usage of MsConvert from the ProteoWizard package
version 3.0.9992 (should also work with new versions) at
http://proteowizard.sourceforge.net/
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

def execute(filelocation, args, outdir, filters=None,
            executable='msConvert.exe'):
    """Execute the msConvert tool on Windows operating systems.

    :param filelocation: input file path
    :param args: str() or list(), msConvert arguments for details see the
        msConvert help below.
    :param outdir: path of the output directory
    :param filters: str() or list(), specify additional parameters and filters,
        for details see the msConvert help below.
    :param executable: must specify the complete file path of the msConvert.exe
        if its location is not in the ``PATH`` environment variable.
    """

    procArgs = [executable, filelocation]
    procArgs.extend(aux.toList(args))
    if filters is not None:
        for arg in aux.toList(filters):
            procArgs.extend(['--filter', arg])
    procArgs.extend(['-o', outdir])

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
    :param outformat: output format, must be one of foloowing: "mzML, mzXML,
        mz5, mgf, text, ms1, cms1, ms2"

    assert outformat in ['mzML', 'mzXML', 'mz5', 'mgf', 'text', 'ms1', 'cms1', 'ms2']
"""

"""
outdir = 'D:/maspy_test/pparse'
filelocation = 'D:/maspy_test/pparse/20160608_QexHF4_RSLCbeta_25ng_HeLa_short_1305_02.raw'
args = ['--mzML', --zlib', '--64']
filters = '"peakPicking true [1,2]"'
execute(filelocation, args, outdir, filters,
             executable='D:/programs/ProteoWizard 3.0.9992/msconvert.exe')
"""


""" MSCONVERT HELP
Usage: msconvert [options] [filemasks]
Convert mass spec data file formats.

Return value: # of failed files.

Options:
  -f [ --filelist ] arg : specify text file containing filenames
  -o [ --outdir ] arg (=.) : set output directory ('-' for stdout) [.]
  -c [ --config ] arg : configuration file (optionName=value)
  --outfile arg : Override the name of output file.
  -e [ --ext ] arg : set extension for output files
  [mzML|mzXML|mgf|txt|mz5]
  --mzML : write mzML format [default]
  --mzXML : write mzXML format
  --mz5 : write mz5 format
  --mgf : write Mascot generic format
  --text : write ProteoWizard internal text format
  --ms1 : write MS1 format
  --cms1 : write CMS1 format
  --ms2 : write MS2 format
  --cms2 : write CMS2 format
  -v [ --verbose ] : display detailed progress information
  --64 : set default binary encoding to 64-bit precision
  [default]
  --32 : set default binary encoding to 32-bit precision
  --mz64 : encode m/z values in 64-bit precision [default]
  --mz32 : encode m/z values in 32-bit precision
  --inten64 : encode intensity values in 64-bit precision
  --inten32 : encode intensity values in 32-bit precision
  [default]
  --noindex : do not write index
  -i [ --contactInfo ] arg : filename for contact info
  -z [ --zlib ] : use zlib compression for binary data
  --numpressLinear [toler] : use numpress linear prediction lossy compression for binary mz and rt data (relative error guaranteed less than given tolerance, default is 2e-009)
  --numpressPic : use numpress positive integer lossy compression for binary intensities (maximum 0.5 absolute error guaranteed)
  --numpressSlof [toler] : use numpress short logged float lossy compression for binary intensities (relative error guaranteed less than given tolerance, default is 0.0002)
  -n [ --numpressAll] : same as --numpressLinear --numpressSlof (see https://github.com/fickludd/ms-numpress for more info)
  -g [ --gzip ] : gzip entire output file (adds .gz to filename)
  --filter arg : add a spectrum list filter
  --merge : create a single output file from multiple input
  files by merging file-level metadata and
  concatenating spectrum lists
  --simAsSpectra : write selected ion monitoring as spectra, not
  chromatograms
  --srmAsSpectra : write selected reaction monitoring as spectra, not
  chromatograms
  --combineIonMobilitySpectra : write all drift bins/scans in a frame/block as one spectrum instead of individual spectra
  --acceptZeroLengthSpectra : some vendor readers have an efficient way of filtering out empty spectra, but it takes more time to open the file
  --ignoreUnknownInstrumentError : if true, if an instrument cannot be determined from a vendor file, it will not be an error
  --help : show this message, with extra detail on filter options

FILTER OPTIONS
run this command with --help to see more detail
index <index_value_set>
msLevel <mslevels>
chargeState <charge_states>
precursorRecalculation
mzRefiner input1.pepXML input2.mzid [msLevels=<1->] [thresholdScore=<CV_Score_Name>] [thresholdValue=<floatset>] [thresholdStep=<float>] [maxSteps=<count>]
lockmassRefiner mz=<real> mzNegIons=<real (mz)> tol=<real (1.0 Daltons)>
precursorRefine
peakPicking [<PickerType> [snr=<minimum signal-to-noise ratio>] [peakSpace=<minimum peak spacing>] [msLevel=<ms_levels>]]
scanNumber <scan_numbers>
scanEvent <scan_event_set>
scanTime <scan_time_range>
sortByScanTime
stripIT
metadataFixer
titleMaker <format_string>
threshold <type> <threshold> <orientation> [<mslevels>]
mzWindow <mzrange>
mzPrecursors <precursor_mz_list>
defaultArrayLength <peak_count_range>
zeroSamples <mode> [<MS_levels>]
mzPresent <tolerance> <type> <threshold> <orientation> <mz_list> [<include_or_exclude>]
scanSumming [precursorTol=<precursor tolerance>] [scanTimeTol=<scan time tolerance>]
MS2Denoise [<peaks_in_window> [<window_width_Da> [multicharge_fragment_relaxation]]]
MS2Deisotope [hi_res [mzTol=<mzTol>]] [Poisson [minCharge=<minCharge>] [maxCharge=<maxCharge>]]
ETDFilter [<removePrecursor> [<removeChargeReduced> [<removeNeutralLoss> [<blanketRemoval> [<matchingTolerance> ]]]]]
chargeStatePredictor [overrideExistingCharge=<true|false (false)>] [maxMultipleCharge=<int (3)>] [minMultipleCharge=<int (2)>] [singleChargeFractionTIC=<real (0.9)>] [maxKnownCharge=<int (0)>] [makeMS2=<true|false (false)>]
turbocharger [minCharge=<minCharge>] [maxCharge=<maxCharge>] [precursorsBefore=<before>] [precursorsAfter=<after>] [halfIsoWidth=<half-width of isolation window>] [defaultMinCharge=<defaultMinCharge>] [defaultMaxCharge=<defaultMaxCharge>] [useVendorPeaks=<useVendorPeaks>]
activation <precursor_activation_type>
analyzer <analyzer>
analyzerType <analyzer>
polarity <polarity>


Examples:

# convert data.RAW to data.mzML
msconvert data.RAW

# convert data.RAW to data.mzXML
msconvert data.RAW --mzXML

# put output file in my_output_dir
msconvert data.RAW -o my_output_dir

# extract scan indices 5...10 and 20...25
msconvert data.RAW --filter "index [5,10] [20,25]"

# extract MS1 scans only
msconvert data.RAW --filter "msLevel 1"

# extract MS2 and MS3 scans only
msconvert data.RAW --filter "msLevel 2-3"

# extract MSn scans for n>1
msconvert data.RAW --filter "msLevel 2-"

# apply ETD precursor mass filter
msconvert data.RAW --filter ETDFilter

# remove non-flanking zero value samples
msconvert data.RAW --filter "zeroSamples removeExtra"

# remove non-flanking zero value samples in MS2 and MS3 only
msconvert data.RAW --filter "zeroSamples removeExtra 2 3"

# add missing zero value samples (with 5 flanking zeros) in MS2 and MS3 only
msconvert data.RAW --filter "zeroSamples addMissing=5 2 3"

# keep only HCD spectra from a decision tree data file
msconvert data.RAW --filter "activation HCD"

# keep the top 42 peaks or samples (depending on whether spectra are centroid or profile):
msconvert data.RAW --filter "threshold count 42 most-intense"

# multiple filters: select scan numbers and recalculate precursors
msconvert data.RAW --filter "scanNumber [500,1000]" --filter "precursorRecalculation"

# multiple filters: apply peak picking and then keep the bottom 100 peaks:
msconvert data.RAW --filter "peakPicking true 1-" --filter "threshold count 100 least-intense"

# multiple filters: apply peak picking and then keep all peaks that are at least 50% of the intensity of the base peak:
msconvert data.RAW --filter "peakPicking true 1-" --filter "threshold bpi-relative .5 most-intense"

# use a configuration file
msconvert data.RAW -c config.txt

# example configuration file
mzXML=true
zlib=true
filter="index [3,7]"
filter="precursorRecalculation"


Questions, comments, and bug reports:
http://proteowizard.sourceforge.net
support@proteowizard.org

ProteoWizard release: 3.0.9134 (2015-11-11)
ProteoWizard MSData: 3.0.9134 (2015-11-11)
ProteoWizard Analysis: 3.0.9098 (2015-11-1)
Build date: Nov 13 2015 11:57:52
"""
