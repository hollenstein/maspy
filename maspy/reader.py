from __future__ import print_function, division, unicode_literals
from future.utils import viewkeys, viewvalues, viewitems, listvalues, listitems

try: # python 2.7
    from itertools import izip as zip
except ImportError: # python 3.x series
    pass
################################################################################
from collections import defaultdict as ddict
import io
from lxml import etree as ETREE
import numpy
from operator import itemgetter as ITEMGETTER
import os
import warnings

import pyteomics.mzid
import pyteomics.mass

import maspy.auxiliary as aux
import maspy.core
import maspy.peptidemethods
import maspy.xml


#####################################################################
### mzml import methods                       #######################
#####################################################################
def importMzml(filepath, msrunContainer=None, siAttrFromSmi=None, specfilename=None):
    """Performs a complete import of a mzml file into a maspy MsrunContainer.

    :paramsiAttrFromSmi: allow here to specify a custom function that extracts params a from spectrumMetadataItem
    :param specfilename: by default the filename will be used as the specfilename in the MsrunContainer and all
        mzML item instances, specify here an alternative specfilename to override the default one
    """
    #TODO: docstring
    siAttrFromSmi = defaultFetchSiAttrFromSmi if siAttrFromSmi is None else siAttrFromSmi
    if msrunContainer is None:
        msrunContainer = maspy.core.MsrunContainer()

    basename = os.path.basename(filepath)
    dirname = os.path.dirname(filepath)
    filename, extension = os.path.splitext(basename)
    specfilename = filename if specfilename is None else specfilename

    #Check if the specified file is valid for an import
    if not os.path.isfile(filepath):
        raise IOError('File does not exist: %s' % filepath)
    elif extension.lower() != '.mzml':
        raise IOError('Filetype is not "mzml": %s' % filepath)
    elif specfilename in msrunContainer.info:
        print(specfilename, 'already present in the msrunContainer, aborting import.')
        return None

    mzmlReader = maspy.xml.MzmlReader(filepath)
    masterContainer = {'rm': str(), 'ci': {}, 'si': {}, 'sai': {}, 'smi': {}}
    #Dictionary recording which MS2 scans follow a MS1 scan
    ms1Record = ddict(list)

    for xmlSpectrum in mzmlReader.parseSpectra():
        smi, binaryDataArrayList = smiFromXmlSpectrum(xmlSpectrum, specfilename)
        #Generate SpectrumItem
        si = maspy.core.Si(smi.id, smi.specfile)
        si.isValid = True
        siAttrFromSmi(smi, si)
        if si.msLevel > 1:
            si.precursorId = si.precursorId.split('scan=')[1] #TODO: change to use regex to extract from known vendor format
            ms1Record[si.precursorId].append(si.id)
        else:
            ms1Record[si.id] #Touch the ddict to add the MS1 id, if it is not already present
        #Generate SpectrumArrayItem
        sai = maspy.core.Sai(smi.id, smi.specfile)
        sai.arrays, sai.arrayInfo = maspy.xml.extractBinaries(binaryDataArrayList,
                                                              smi.attributes['defaultArrayLength'])
        #Store all items in the appropriate containers
        masterContainer['smi'][smi.id] = smi
        masterContainer['si'][smi.id] = si
        masterContainer['sai'][smi.id] = sai

    for siId, msnIdList in viewitems(ms1Record):
        #Ignore KeyError if the spectrum is not present in the mzML file for whatever reason
        try:
            setattr(masterContainer['si'][siId], 'msnIdList', msnIdList)
        except KeyError:
            pass

    for xmlChromatogram in mzmlReader.chromatogramList:
        ci = ciFromXml(xmlChromatogram, specfilename)
        masterContainer['ci'][ci.id] = ci
    masterContainer['rm'] = mzmlReader.metadataNode

    msrunContainer._addSpecfile(specfilename, dirname)
    msrunContainer.rmc[specfilename] = masterContainer['rm']
    msrunContainer.info[specfilename]['status']['rm'] = True
    msrunContainer.smic[specfilename] = masterContainer['smi']
    msrunContainer.info[specfilename]['status']['smi'] = True
    msrunContainer.sic[specfilename] = masterContainer['si']
    msrunContainer.info[specfilename]['status']['si'] = True
    msrunContainer.saic[specfilename] = masterContainer['sai']
    msrunContainer.info[specfilename]['status']['sai'] = True
    msrunContainer.cic[specfilename] = masterContainer['ci']
    msrunContainer.info[specfilename]['status']['ci'] = True

    return msrunContainer


# --- generate python objects from mzML xml elements --- #
def ciFromXml(xmlelement, specfile):
    #TODO: docstring
    ciId = xmlelement.attrib['id']
    ci = maspy.core.Ci(ciId, specfile)
    for key in ['dataProcessingRef']:
        if key in xmlelement.attrib:
            setattr(ci, key, xmlelement.attrib[key])
    ci.params, children = maspy.xml.extractParams(xmlelement)
    arrayLength = xmlelement.attrib['defaultArrayLength']

    for child in children:
        childElements, childParams = maspy.xml.sublistReader(child)
        childTag = maspy.xml.clearTag(child.tag)
        if childTag == 'precursor':
            #TODO: THIS HAS TO BE ADAPTED, REPLACE "**element", ADD A LOOP
            newElement = MzmlPrecursor(**element)
        elif childTag == 'product':
            #TODO: THIS HAS TO BE ADAPTED
            newElement = MzmlProduct(**element)
        elif childTag == 'binaryDataArrayList':
            for childElement in childElements:
                dataType, dataTypeParam = maspy.xml.findBinaryDataType(childElement['params'])
                ci.arrayInfo[dataType] = {'dataProcessingRef': None, 'params': childElement['params']}
                if 'dataProcessingRef' in childElement:
                    ci.arrayInfo[dataType]['dataProcessingRef'] = childElement['dataProcessingRef']
            ci.arrays, ci.arrayInfo = maspy.xml.extractBinaries(childElements, arrayLength)
    return ci


def smiFromXmlSpectrum(xmlelement, specfile):
    #TODO: docstring
    scanId = xmlelement.attrib['id'].split('scan=')[1] #TODO: change to use regex to extract from known vendor format
    smi = maspy.core.Smi(scanId, specfile)
    smi.params, spectrumChildren = maspy.xml.extractParams(xmlelement)
    smi.attributes.update(xmlelement.attrib)
    binaryDataArrayList = list()

    for spectrumChild in spectrumChildren:
        spectrumChildTag = maspy.xml.clearTag(spectrumChild.tag)
        elements, params = maspy.xml.sublistReader(spectrumChild)
        if params:
            #Define scanListParams here
            setattr(smi, spectrumChildTag+'Params', params)
        for element in elements:
            if spectrumChildTag == 'scanList':
                newElement = maspy.core.MzmlScan(**element)
            elif spectrumChildTag == 'precursorList':
                newElement = maspy.core.MzmlPrecursor(**element)
            elif spectrumChildTag == 'productList':
                newElement = maspy.core.MzmlProduct(**element)
            elif spectrumChildTag == 'binaryDataArrayList':
                binaryDataArrayList.append(element)
                continue
            else:
                raise Exception('unexpected spectrum xml child tag!')
            getattr(smi, spectrumChildTag).append(newElement)
    return smi, binaryDataArrayList


# --- extract spectrum attributes from the mzml like structure of a Smi object --- #
def fetchSpectrumInfo(smi):
    #TODO: docstring
    attributes = dict()
    for key, accession, dtype in [('msLevel', 'MS:1000511', int),
                                  ('tic', 'MS:1000285', float),
                                  ('basepeakMz', 'MS:1000504', float),
                                  ('basepeakI', 'MS:1000505', float)
                                  ]:
        param = maspy.xml.findParam(smi.params, accession)
        if param is not None:
            attributes[key] = dtype(param[1])
        else:
            #Or don't add key to attributes?
            attributes[key] = None
    return attributes


def fetchScanInfo(smi):
    #TODO: docstring
    if smi.scanList:
        attributes = dict()
        for key, accession, dtype in [('rt', 'MS:1000016', float),
                                      ('iit', 'MS:1000927', float),
                                      ]:
            param = maspy.xml.findParam(smi.scanList[0]['params'], accession)
            if param is not None:
                attributes[key] = dtype(param[1])
                #TODO: raise a warning if time is neither in minutes nor in seconds
                if key == 'rt' and param[2] == 'UO:0000031':
                    #convert retention time from minutes to seconds
                    attributes['rt'] *= 60
            else:
                attributes[key] = None
    else:
        attributes = None
    return attributes


def fetchParentIon(smi):
    #TODO: docstring
    if smi.precursorList:
        attributes = dict()
        attributes['precursorId'] = smi.precursorList[0].spectrumRef
        selectedIon = smi.precursorList[0]['selectedIonList'][0]
        for key, accession, dtype in [('obsMz', 'MS:1000744', float),
                                      ('i', 'MS:1000042', float),
                                      ('charge', 'MS:1000041', int)
                                      ]:
            param = maspy.xml.findParam(selectedIon, accession)
            if param is not None:
                attributes[key] = dtype(param[1])
            else:
                attributes[key] = None
    else:
        attributes = None
    return attributes


def defaultFetchSiAttrFromSmi(smi, si):
    """Default method to extract attributes from a spectrum metadata item (sai) and
    adding them to a spectrum item (si)."""
    for key, value in viewitems(fetchSpectrumInfo(smi)):
        setattr(si, key, value)
    for key, value in viewitems(fetchScanInfo(smi)):
        setattr(si, key, value)
    if si.msLevel > 1:
        for key, value in viewitems(fetchParentIon(smi)):
            setattr(si, key, value)


def convertMzml(mzmlPath, outputDirectory=None):
    """Imports an mzml file and converts it to a MsrunContainer file

    :param mzmlPath: path of the mzml file
    :param outputDirectory: directory where the MsrunContainer file should be written
    if it is not specified, the output directory is set to the mzml files directory.
    """
    outputDirectory = outputDirectory if outputDirectory is not None else os.path.dirname(mzmlPath)
    msrunContainer = importMzml(mzmlPath)
    msrunContainer.setPath(outputDirectory)
    msrunContainer.save()


#####################################################################
### import for SiiContainer class             #######################
#####################################################################
def prepareSiiImport(siiContainer, specfile, path, qcAttr, qcLargerBetter,
                     qcCutoff, rankAttr, rankLargerBetter):
    """Prepares the ``siiContainer`` for the import of peptide spectrum matching
    results. Adds entries to ``siiContainer.container`` and to
    ``siiContainer.info``.

    :param siiContainer: instance of :class:`maspy.core.SiiContainer`
    :param specfile: unambiguous identifier of a ms-run file. Is also used as
        a reference to other MasPy file containers.
    :param path: folder location used by the ``SiiContainer`` to save and load
        data to the hard disk.
    :param qcAttr: name of the parameter to define a ``Sii`` quality cut off.
        Typically this is some sort of a global false positive estimator,
        for example a 'false discovery rate' (FDR).
    :param qcLargerBetter: bool, True if a large value for the ``.qcAttr`` means
        a higher confidence.
    :param qcCutOff: float, the quality threshold for the specifed ``.qcAttr``
    :param rankAttr: name of the parameter used for ranking ``Sii`` according
        to how well they match to a fragment ion spectrum, in the case when
        their are multiple ``Sii`` present for the same spectrum.
    :param rankLargerBetter: bool, True if a large value for the ``.rankAttr``
        means a better match to the fragment ion spectrum.

    For details on ``Sii`` ranking see :func:`applySiiRanking()`

    For details on ``Sii`` quality validation see :func:`applySiiQcValidation()`
    """
    if specfile not in siiContainer.info:
        siiContainer.addSpecfile(specfile, path)
    else:
        raise Exception('...')

    siiContainer.info[specfile]['qcAttr'] = qcAttr
    siiContainer.info[specfile]['qcLargerBetter'] = qcLargerBetter
    siiContainer.info[specfile]['qcCutoff'] = qcCutoff
    siiContainer.info[specfile]['rankAttr'] = rankAttr
    siiContainer.info[specfile]['rankLargerBetter'] = rankLargerBetter


def addSiiToContainer(siiContainer, specfile, siiList):
    """Adds the ``Sii`` elements contained in the siiList to the appropriate
    list in ``siiContainer.container[specfile]``.

    :param siiContainer: instance of :class:`maspy.core.SiiContainer`
    :param specfile: unambiguous identifier of a ms-run file. Is also used as
        a reference to other MasPy file containers.
    :param siiList: a list of ``Sii`` elements imported from any PSM search
        engine results
    """
    for sii in siiList:
        if sii.id not in siiContainer.container[specfile]:
            siiContainer.container[specfile][sii.id] = list()
        siiContainer.container[specfile][sii.id].append(sii)


def applySiiRanking(siiContainer, specfile):
    """Iterates over all Sii entries of a specfile in siiContainer and sorts Sii
    elements of the same spectrum according to the score attribute specified in
    ``siiContainer.info[specfile]['rankAttr']``. Sorted Sii elements are then
    ranked  according to their sorted position, if multiple Sii have the same
    score, all get the same rank and the next entries rank is its list position.

    :param siiContainer: instance of :class:`maspy.core.SiiContainer`
    :param specfile: unambiguous identifier of a ms-run file. Is also used as
        a reference to other MasPy file containers.
    """
    attr = siiContainer.info[specfile]['rankAttr']
    reverse = siiContainer.info[specfile]['rankLargerBetter']
    for itemList in listvalues(siiContainer.container[specfile]):
        sortList = [(getattr(sii, attr), sii) for sii in itemList]
        itemList = [sii for score, sii in sorted(sortList, reverse=reverse)]

        #Rank Sii according to their position
        lastValue = None
        for itemPosition, item in enumerate(itemList, 1):
            if getattr(item, attr) != lastValue:
                rank = itemPosition
            item.rank = rank
            lastValue = getattr(item, attr)


def applySiiQcValidation(siiContainer, specfile):
    """Iterates over all Sii entries of a specfile in siiContainer and validates
    if they surpass a user defined quality threshold. The parameters for
    validation are defined in ``siiContainer.info[specfile]``:

        - ``qcAttr``, ``qcCutoff`` and ``qcLargerBetter``

    In addition to passing this validation a ``Sii`` has also to be at the first
    list position in the ``siiContainer.container``. If both criteria are met
    the attribute ``Sii.isValid`` is set to ``True``.

    :param siiContainer: instance of :class:`maspy.core.SiiContainer`
    :param specfile: unambiguous identifier of a ms-run file. Is also used as
        a reference to other MasPy file containers.
    """
    attr = siiContainer.info[specfile]['qcAttr']
    cutOff = siiContainer.info[specfile]['qcCutoff']
    if siiContainer.info[specfile]['qcLargerBetter']:
        evaluator = lambda sii: getattr(sii, attr) >= cutOff and sii.rank == 1
    else:
        evaluator = lambda sii: getattr(sii, attr) <= cutOff and sii.rank == 1

    for itemList in listvalues(siiContainer.container[specfile]):
        #Set the .isValid attribute of all Sii to False
        for sii in itemList:
            sii.isValid = False

        #Validate the first Sii
        sii = itemList[0]
        if evaluator(sii):
            sii.isValid = True


def readPercolatorResults(filelocation, specfile, psmEngine):
    """Reads percolator PSM results from a txt file and returns a list of
    :class:`Sii <maspy.core.Sii>` elements.

    :param filelocation: file path of the percolator result file
    :param specfile: unambiguous identifier of a ms-run file. Is also used as
        a reference to other MasPy file containers.
    :param psmEngine: PSM PSM search engine used for peptide spectrum matching
        before percolator. This is important to specify, since the scanNr
        information is written in a different format by some engines. It might
        be necessary to adjust the settings for different versions of percolator
        or the PSM search engines used.

        Possible values are 'comet', 'xtandem', 'msgf'.

    :returns: [sii, sii, sii, ...]
    """
    if psmEngine not in ['comet', 'msgf', 'xtandem']:
        raise Exception('PSM search engine not supported: ', psmEngine)
    itemList = list()

    #Note: regarding headerline, xtandem seperates proteins with ';',
    #msgf separates proteins with a tab
    with io.open(filelocation, 'r', encoding='utf-8') as openfile:
        lines = openfile.readlines()

    headerDict = dict([[y,x] for (x,y) in
                       enumerate(lines[0].strip().split('\t'))
                       ])
    scanEntryList = list()
    for line in lines[1:]:
        if len(line.strip()) == 0:
            continue
        fields = line.strip().split('\t')

        if psmEngine in ['comet', 'msgf']:
            scanNr = fields[headerDict['PSMId']].split('_')[-3]
        elif psmEngine in ['xtandem']:
            scanNr = fields[headerDict['PSMId']].split('_')[-2]

        peptide = fields[headerDict['peptide']]
        if peptide.find('.') != -1:
            peptide = peptide.split('.')[1]
        #Change to the new unimod syntax
        peptide = peptide.replace('[UNIMOD:', '[u:')
        sequence = maspy.peptidemethods.removeModifications(peptide)

        qValue = fields[headerDict['q-value']]
        score = fields[headerDict['score']]
        pep = fields[headerDict['posterior_error_prob']]

        sii = maspy.core.Sii(scanNr, specfile)
        sii.peptide = peptide
        sii.sequence = sequence
        sii.qValue = float(qValue)
        sii.score = float(score)
        sii.pep = float(pep)
        sii.isValid = False

        itemList.append(sii)
    return itemList


def importPercolatorResults(siiContainer, filelocation, specfile, psmEngine,
                            qcAttr='qValue', qcLargerBetter=False,
                            qcCutoff=0.01, rankAttr='score',
                            rankLargerBetter=True):
    """Import peptide spectrum matches (PSMs) from a percolator result file,
    generate :class:`Sii <maspy.core.Sii>` elements and store them in the
    specified :class:`siiContainer <maspy.core.SiiContainer>`. Imported ``Sii``
    are ranked according to a specified attribute and validated if they surpass
    a specified quality threshold.

    :param siiContainer: imported PSM results are added to this instance of
        :class:`siiContainer <maspy.core.SiiContainer>`
    :param filelocation: file path of the percolator result file
    :param specfile: unambiguous identifier of a ms-run file. Is also used as
        a reference to other MasPy file containers.
    :param psmEngine: PSM search engine used for peptide spectrum matching
        before percolator. For details see :func:`readPercolatorResults()`.
        Possible values are 'comet', 'xtandem', 'msgf'.
    :param qcAttr: name of the parameter to define a quality cut off. Typically
        this is some sort of a global false positive estimator (eg FDR)
    :param qcLargerBetter: bool, True if a large value for the ``.qcAttr`` means
        a higher confidence.
    :param qcCutOff: float, the quality threshold for the specifed ``.qcAttr``
    :param rankAttr: name of the parameter used for ranking ``Sii`` according
        to how well they match to a fragment ion spectrum, in the case when
        their are multiple ``Sii`` present for the same spectrum.
    :param rankLargerBetter: bool, True if a large value for the ``.rankAttr``
        means a better match to the fragment ion spectrum

    For details on ``Sii`` ranking see :func:`applySiiRanking()`

    For details on ``Sii`` quality validation see :func:`applySiiQcValidation()`
    """

    path = os.path.dirname(filelocation)
    siiList = readPercolatorResults(filelocation, specfile, psmEngine)
    prepareSiiImport(siiContainer, specfile, path, qcAttr, qcLargerBetter,
                     qcCutoff, rankAttr, rankLargerBetter)
    addSiiToContainer(siiContainer, specfile, siiList)
    applySiiRanking(siiContainer, specfile)
    applySiiQcValidation(siiContainer, specfile)


def readMsgfMzidResults(filelocation, specfile=None):
    """Reads MS-GF+ PSM results from a mzIdentML file and returns a list of
    :class:`Sii <maspy.core.Sii>` elements.

    :param filelocation: file path of the percolator result file
    :param specfile: optional, unambiguous identifier of a ms-run file. Is also
        used as a reference to other MasPy file containers. If specified all
        the ``.specfile`` attribute of all ``Sii`` are set to this value, else
        it is read from the mzIdentML file.

    :returns: [sii, sii, sii, ...]
    """
    readSpecfile = True if specfile is None else False

    unimod = pyteomics.mass.mass.Unimod()
    _tempMods = dict()
    mzid_refs = pyteomics.mzid.read(filelocation, retrieve_refs=True,
                                    iterative=False)

    siiList = list()
    for mzidEntry in mzid_refs:
        mzidSii = mzidEntry['SpectrumIdentificationItem'][0]
        scanNumber = str(int(mzidEntry['scan number(s)']))
        if readSpecfile:
            specfile = os.path.splitext(mzidEntry['name'])[0]

        sii = maspy.core.Sii(scanNumber, specfile)
        sii.isValid = mzidSii['passThreshold']
        sii.rank = mzidSii['rank']
        sii.eValue = mzidSii['MS-GF:EValue']
        sii.charge = mzidSii['chargeState']
        sii.sequence = mzidSii['PeptideSequence']
        sii.specEValue = mzidSii['MS-GF:SpecEValue']
        sii.score = numpy.log10(sii.eValue)*-1

        if 'Modification' in mzidSii:
            modifications = list()
            for modEntry in mzidSii['Modification']:
                try:
                    modSymbolMaspy = _tempMods[modEntry['name']]
                except KeyError:
                    unimodEntry = unimod.by_title(modEntry['name'])
                    if len(unimodEntry) != 0:
                        modSymbol = 'u:'+str(unimodEntry['record_id'])
                    else:
                        modSymbol = modEntry['name']
                    modSymbolMaspy = '[' + modSymbol + ']'
                    _tempMods[modEntry['name']] = modSymbolMaspy
                modifications.append((modEntry['location'], modSymbolMaspy))
            modifications.sort(key=ITEMGETTER(0))

            _lastPos = 0
            _peptide = list()
            for pos, mod in modifications:
                 _peptide.extend((sii.sequence[_lastPos:pos], mod))
                 _lastPos = pos
            _peptide.append(sii.sequence[_lastPos:])

            sii.peptide = ''.join(_peptide)
        else:
            sii.peptide = sii.sequence
        siiList.append(sii)
    return siiList


def importMsgfMzidResults(siiContainer, filelocation, specfile=None,
                          qcAttr='eValue', qcLargerBetter=False, qcCutoff=0.01,
                          rankAttr='score', rankLargerBetter=True):
    """Import peptide spectrum matches (PSMs) from a MS-GF+ mzIdentML file,
    generate :class:`Sii <maspy.core.Sii>` elements and store them in the
    specified :class:`siiContainer <maspy.core.SiiContainer>`. Imported ``Sii``
    are ranked according to a specified attribute and validated if they surpass
    a specified quality threshold.

    :param siiContainer: imported PSM results are added to this instance of
        :class:`siiContainer <maspy.core.SiiContainer>`
    :param filelocation: file path of the percolator result file
    :param specfile: optional, unambiguous identifier of a ms-run file. Is also
        used as a reference to other MasPy file containers. If specified the
        attribute ``.specfile`` of all ``Sii`` is set to this value, else
        it is read from the mzIdentML file.
    :param qcAttr: name of the parameter to define a quality cut off. Typically
        this is some sort of a global false positive estimator (eg FDR)
    :param qcLargerBetter: bool, True if a large value for the ``.qcAttr`` means
        a higher confidence.
    :param qcCutOff: float, the quality threshold for the specifed ``.qcAttr``
    :param rankAttr: name of the parameter used for ranking ``Sii`` according
        to how well they match to a fragment ion spectrum, in the case when
        their are multiple ``Sii`` present for the same spectrum.
    :param rankLargerBetter: bool, True if a large value for the ``.rankAttr``
        means a better match to the fragment ion spectrum

    For details on ``Sii`` ranking see :func:`applySiiRanking()`

    For details on ``Sii`` quality validation see :func:`applySiiQcValidation()`
    """
    path = os.path.dirname(filelocation)
    siiList = readMsgfMzidResults(filelocation, specfile)

    #If the mzIdentML file contains multiple specfiles, split the sii elements
    # up according to their specified "specfile" attribute.
    specfiles = ddict(list)
    if specfile is None:
        for sii in siiList:
            specfiles[sii.specfile].append(sii)
    else:
        specfiles[specfile] = siiList

    for specfile in specfiles:
        _siiList = specfiles[specfile]

        prepareSiiImport(siiContainer, specfile, path, qcAttr, qcLargerBetter,
                         qcCutoff, rankAttr, rankLargerBetter)
        addSiiToContainer(siiContainer, specfile, _siiList)
        applySiiRanking(siiContainer, specfile)
        applySiiQcValidation(siiContainer, specfile)


#####################################################################
### import for FeatureContainer class         #######################
#####################################################################
def importPeptideFeatures(fiContainer, filelocation, specfile):
    """ Import peptide features from a featureXml file, as generated for example
    by the OpenMS node featureFinderCentroided, or a features.tsv file by the
    Dinosaur command line tool.

    :param fiContainer: imported features are added to this instance of
        :class:`FeatureContainer <maspy.core.FeatureContainer>`.
    :param filelocation: Actual file path
    :param specfile: Keyword (filename) to represent file in the
        :class:`FeatureContainer`. Each filename can only occure once, therefore
        importing the same filename again is prevented.
    """
    if not os.path.isfile(filelocation):
        warnings.warn('The specified file does not exist %s' %(filelocation, ))
        return None
    elif (not filelocation.lower().endswith('.featurexml') and
          not filelocation.lower().endswith('.features.tsv')
          ):
        #TODO: this is depricated as importPeptideFeatues
        #is not longer be used solely for featurexml
        print('Wrong file extension, %s' %(filelocation, ))
    elif specfile in fiContainer.info:
        print('%s is already present in the SiContainer, import interrupted.'
              %(specfile, )
              )
        return None

    #Prepare the file container for the import
    fiContainer.addSpecfile(specfile, os.path.dirname(filelocation))

    #import featurexml file
    if filelocation.lower().endswith('.featurexml'):
        featureDict = _importFeatureXml(filelocation)
        for featureId, featureEntryDict in viewitems(featureDict):
            rtArea = set()
            for convexHullEntry in featureEntryDict['convexHullDict']['0']:
                rtArea.update([convexHullEntry[0]])

            fi = maspy.core.Fi(featureId, specfile)
            fi.rt = featureEntryDict['rt']
            fi.rtArea = max(rtArea) - min(rtArea)
            fi.rtLow = min(rtArea)
            fi.rtHigh = max(rtArea)
            fi.charge = featureEntryDict['charge']
            fi.mz = featureEntryDict['mz']
            fi.mh = maspy.peptidemethods.calcMhFromMz(featureEntryDict['mz'],
                                                      featureEntryDict['charge'])
            fi.intensity = featureEntryDict['intensity']
            fi.quality = featureEntryDict['overallquality']
            fi.isMatched = False
            fi.isAnnotated = False
            fi.isValid = True

            fiContainer.container[specfile][featureId] = fi

    #import dinosaur tsv file
    elif filelocation.lower().endswith('.features.tsv'):
        featureDict = _importDinosaurTsv(filelocation)
        for featureId, featureEntryDict in viewitems(featureDict):
            fi = maspy.core.Fi(featureId, specfile)
            fi.rt = featureEntryDict['rtApex']
            fi.rtArea = featureEntryDict['rtEnd'] - featureEntryDict['rtStart']
            fi.rtFwhm = featureEntryDict['fwhm']
            fi.rtLow = featureEntryDict['rtStart']
            fi.rtHigh = featureEntryDict['rtEnd']
            fi.charge = featureEntryDict['charge']
            fi.mz = featureEntryDict['mz']
            fi.mh = maspy.peptidemethods.calcMhFromMz(featureEntryDict['mz'],
                                                      featureEntryDict['charge'])
            fi.intensity = featureEntryDict['intensitySum']
            fi.intensityApex = featureEntryDict['intensityApex']
            #Note: not used keys:
            #mostAbundantMz nIsotopes nScans averagineCorr mass massCalib

            fi.isMatched = False
            fi.isAnnotated = False
            fi.isValid = True

            fiContainer.container[specfile][featureId] = fi


def _importFeatureXml(filelocation):
    """Reads a featureXml file.

    :param filelocation: #TODO: docstring
    :returns: {featureKey1: {attribute1:value1, attribute2:value2, ...}, ...}

    See also :func:`importPeptideFeatures`
    """
    with io.open(filelocation, 'r', encoding='utf-8') as openFile:
        readingFeature = False
        readingHull = False
        featureDict = dict()

        for i, line in enumerate(openFile):
            line = line.strip()
            if readingFeature == True:
                if line.find('<convexhull') != -1:
                    readingHull = True
                    hullNr = line.split('<convexhull nr=\"')[1].split('\">')[0]
                    hullList = list()
                elif readingHull == True:
                    if line.find('<pt') != -1:
                        x = float(line.split('x=\"')[1].split('\"')[0])
                        y = float(line.split('y=\"')[1].split('\"')[0])
                        # x = retentiontime, y = m/z
                        #retentionTimeList.append(x)
                        hullList.append([x,y])
                    elif line.find('</convexhull>') != -1:
                        featureDict[featureKey]['convexHullDict'][hullNr] = hullList
                        readingHull = False

                elif line.find('<position dim=\"0\">') != -1:
                    featureDict[featureKey]['dim0'] = float(line.split('<position dim=\"0\">')[1].split('</position>')[0])
                elif line.find('<position dim=\"1\">') != -1:
                    featureDict[featureKey]['dim1'] = float(line.split('<position dim=\"1\">')[1].split('</position>')[0])
                elif line.find('<intensity>') != -1:
                    featureDict[featureKey]['intensity'] = float(line.split('<intensity>')[1].split('</intensity>')[0])
                elif line.find('<overallquality>') != -1:
                    featureDict[featureKey]['overallquality'] = float(line.split('<overallquality>')[1].split('</overallquality>')[0])
                elif line.find('<charge>') != -1:
                    featureDict[featureKey]['charge'] = int( line.split('<charge>')[1].split('</charge>')[0] )

                elif line.find('<userParam') != -1:
                    if line.find('name=\"label\"') != -1:
                        featureDict[featureKey]['label'] = line.split('value=\"')[1].split('\"/>')[0]
                    elif line.find('name=\"score_fit\"') != -1:
                        featureDict[featureKey]['score_fit'] = float(line.split('value=\"')[1].split('\"/>')[0])
                    elif line.find('name=\"score_correlation\"') != -1:
                        featureDict[featureKey]['score_correlation'] = float(line.split('value=\"')[1].split('\"/>')[0])
                    elif line.find('name=\"FWHM\"') != -1:
                        featureDict[featureKey]['FWHM'] = float(line.split('value=\"')[1].split('\"/>')[0])
                    elif line.find('name=\"spectrum_index\"') != -1:
                        featureDict[featureKey]['spectrum_index'] = line.split('value=\"')[1].split('\"/>')[0]
                    elif line.find('name=\"spectrum_native_id\"') != -1:
                        featureDict[featureKey]['spectrum_native_id'] = line.split('value=\"')[1].split('\"/>')[0]

                elif line.find('</feature>') != -1:
                    #mzList = list()
                    #for retentionTime,mz in featureDict[featureKey]['convexHullDict']['0']:
                    #    mzList.append(mz)
                    featureDict[featureKey]['rt'] = featureDict[featureKey]['dim0']#numpy.median(retentionTimeList)
                    featureDict[featureKey]['mz'] = featureDict[featureKey]['dim1']#numpy.median(mzList)

                    readingFeature == False

            if line.find('<feature id') != -1:
                readingFeature = True
                featureKey = line.split('<feature id=\"')[1].split('\">')[0]
                featureDict[featureKey] = dict()
                featureDict[featureKey]['convexHullDict'] = dict()
                #retentionTimeList = list()
    return featureDict


def _importDinosaurTsv(filelocation):
    """Reads a Dinosaur tsv file.

    :returns: {featureKey1: {attribute1:value1, attribute2:value2, ...}, ...}

    See also :func:`importPeptideFeatures`
    """
    with io.open(filelocation, 'r', encoding='utf-8') as openFile:
        #NOTE: this is pretty similar to importing percolator results, maybe unify in a common function
        lines = openFile.readlines()
        headerDict = dict([[y,x] for (x,y) in enumerate(lines[0].strip().split('\t'))])
        featureDict = dict()
        for featureId, line in enumerate(lines[1:]):
            fields = line.strip().split('\t')
            entryDict = dict()
            for headerName, headerPos in viewitems(headerDict):
                entryDict[headerName] = float(fields[headerPos])
                if headerName in ['rtApex', 'rtEnd', 'rtStart', 'fwhm']:
                    #Covnert to seconds
                    entryDict[headerName] *= 60
                elif headerName in ['charge', 'intensitySum', 'nIsotopes', 'nScans', 'intensityApex']:
                    entryDict[headerName] = int(entryDict[headerName])
            featureDict[featureId] = entryDict
    return featureDict
