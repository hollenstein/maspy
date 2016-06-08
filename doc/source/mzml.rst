MS spectra in MasPy - basics
----------------------------

Every vendor software produces mass spectrometer output files in a different
proprietary format. It is a difficult and time consuming task for software
developers to support all of these different formats and the different
versions thereof. Therefore the file format mzML has been developed as the
community standard for represenation of mass spectrometry results. mzML is an
open, XML-based format that not only allows to store recorded mass spectrum
information but also metadata of the instrument configuration, acquisition
settings, software used for data processing and sample descriptions. (see also
`Mass Spectrometer Output File Format mzML
<http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3073315>`_)

.. note::
    Refer to `www.psidev.info <http://www.psidev.info/index.php?q=node/257>`_
    for details on the XML schema definition and mzML file specifications.

Ultimately, it is desirable to use mzML for archiving, sharing, and processing
of mass spectrometry data and thus for all software to support and use the
mzML format.

.. note::
    We recommend using ProteoWizard for conversion of vendor format files to
    mzML. The software can be downloaded from their `website
    <http://proteowizard.sourceforge.net>`_, a detailed protocol how to use
    ProteoWizard can be found `here
    <http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4113728>`_.

    The raw spectral data recorded by an instrument can be either stored as
    profile or centroid data. In centroid mode each signal peak is present
    only as one single pair of a dinstinct m/z value and an intensity. The
    process of converting profile data to centroid data is called peak picking
    and can be applied as a filter while converting vendor format files to
    mzML files using ProteoWizard. The represenation as centroid data is
    easier to work with, saves memory and is sufficient for most applications.
    Therefore we recommend the utilization of centroid data for MasPy.

Modern mass spectrometers can generate tens of thousands spectra per hour
resulting in huge mzML files. Opening and parsing of such large xml files
takes a lot of time. If only some of the spectra have to be accessed at a
certain time, this can be solved partly by using an indexed mzML file.

MsrunContainer
^^^^^^^^^^^^^^

In MasPy we decided to go one step further and split the information that is
present in mzML files into four different groups; run metadata, spectrum
metadata items, spectrum array items and chromatogram items. Each of these
data groups is stored separately in MasPy and has its own file type, thus it
can be accessed, saved and loaded independently of the others. This allows
convenient and very fast reading if only a certain part of the data needs to
be accessed. We call this represenation of an mzML file in MasPy an
:class:`MsrunContainer <maspy.core.MsrunContainer>`. Altough the data is split
into multiple parts, all information originally contained in an mzML file is
still present. This allows the conversion from MsrunContainer to mzML at any
given time (altough the generation of indexed mzML files is not yet
implemented and we rely on ProteoWizard to add an index for now).

See tutorial/docstrings xxx for details on the MsrunContainer file
format.

Fig.: MsrunContainer

* run metadata
* spectrum metadata items
* spectrum array items
* chromatogram items
* spectrum items


Run metadata (rm)
^^^^^^^^^^^^^^^^^

The run metadata element contains all metadata information of the mzML file,
which is not part of the aquired spectra and chromatograms. These mzML
metadata elements are stored as (xml nodes in) a ``lxml.etree.Element``. The
MasPy run metadata element thus holds an exact copy of the respective mzML
elements.

Note: Software which is used to process data of an mzML file should be listed
in the mzML element "softwareList", and the applied data processing steps
should be documented in the "dataProcessingList" element.


Spectrum array item (:class:`Sai <maspy.core.Sai>`), spectrum metadata item (:class:`Smi <maspy.core.Smi>`)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

An mzML spectrum element contains all information of an aquired MS spectrum,
including numerical arrays containing at least recorded m/z and intensity
values of the observed ions but also plenty of metadata data describing for
example details of the aquisition like base peak m/z and intensity, scan start
time, ms level, or precursor information of MS2 scans. In MasPy this
information is split into metadata and spectrum array information and put in
two separate elements(objects better?); spectrum metadata item
(:class:`Smi<maspy.core.Smi>`) and spectraum array item
(:class:`Sai<maspy.core.Sai>`), respectively. ``Smi`` elements are stored in
``MsrunContainer.smic`` (Smi container) and ``Sai`` elements in
``MsrunContainer.saic`` (Sai container). In order to generate an mzML style
spectrum xml element the information of both MasPy elements (``Smi`` and
``Sai``) is required.


Chromatogram item (:class:`maspy.core.Ci`)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A mzML chromatogram element is similar to a spectrum element, containing
metadata and numerical arrays, in which one dimension is typically a time
series. In the current MasPy implemenatation chromatogram elements are not
split but the metadata and chromatogram array information are put in one
single element called chromatogram item (:class:`maspy.core.Ci`), which is
stored in ``MsrunContainer.cic`` (Ci container).


Spectrum item (:class:`maspy.core.Si`)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The mzML file  serves as a data container for active data processing but also
for data sharing and archiving. Thus the spectrum elements contain a lot of
metadata information not needed for most data analysis application. In
addition all information stored in spectrum elements have to be in accordance
with the mzML xml scheme definition and the Controlled Vocabularies (CV's) of
the Proteomic Standard Initiative (`link <http://www.psidev.info/groups
/controlled-vocabularies>`_). Altough in principle this standardization is a
good thing and perfectly reasonable, when actively working with the data this
can be unnecessary and make things quite complicated.

To circumvent this problem MasPy provides a simplier data type for working
with spectrum metdata, called spectrum item (:class:`Si <maspy.core.Si>`). The
``Si`` class has a flat structure, meaning that attributes are not nested
inside other elements but are stored directly as attributes of the class.
``Si`` attributes can be manipulated without restrictions and new attributes
can simply be added. Specific functions can be used to selectively extract
information from ``Smi``. This allows to only import the currently needed
spectrum metadata attributes, thereby making the ``Si`` more memory efficient.
In order to make lasting changes to the mzML file data ``Si`` attributes have
to be translated to the respective ``Smi``. These changes however have to
strictly follow the mzML specifications and syntax. Thus it is recommend to
use existing functions or implement new ones which make changes to ``Smi``
elements in a controlled way.

Each spectrum present in an mzML file is therefore represented threefold in
MasPy. First the ``Smi`` contains a complete representation of all metadata
information present in an mzML spectrum element. Second the ``Sai`` contains
the actual ion information recorded by the mass spectrometer. And third the
``Si``, which can be considered as the spectrum metadata workspace in MasPy,
allowing convinient access to metadata and simple processing of this data
without directly altering the original mzML information.

*MsrunContainer.info -> which specfiles are present, what is to current path
(used for loading or saving) , which data types are currently imported*


Basic code examples
^^^^^^^^^^^^^^^^^^^

**Importing an mzML file**

mzML files can be imported by using the function
:func:`maspy.reader.importMzml()`, the imported specfile is then added to the
``MsrunContainer`` instance passed to the function. ::

    import maspy.core
    import maspy.reader

    mzmlfilepath = 'filedirectory/specfile_name_1.mzML'
    msrunContainer = maspy.core.MsrunContainer()
    maspy.reader.importMzml(mzmlfilepath, msrunContainer)


**Saving an MsrunContainer to the hard disk**

An ``MsrunContainer`` can be saved to the hard disk by calling its
:func:`.save() <maspy.core.MsrunContainer.save>` method. ::

    msrunContainer.save()

By default all files are saved into the folder specified in ``.info``. This can
be altered by changing the ``path`` variable in ``.info`` or temporarely by
passing the "path" parameter to ``.save()``. ::

    msrunContainer.save(path='../an_alternative_location')

In addition, multiple parameters can be set to specify which part of the data
should be written to the hard disk. The keywords "rm", "ci", "smi", "sai" and
"si" can be set to ``True`` or ``False`` and specify which container types are
selected for saving. By default all of them are set to ``False`` which is
however interpreted as selecting all of them. Setting at least one to ``True``
changes this behaviour and only the specified ones are selected. If multiple
specfiles are present in an ``MsrunContainer`` it is possible to only select a
subset for saving by passing the "specfiles" argument to ``.save()``. The value
of "specfiles" can either be the name of one single specfile are a list of
specfile names. In the following example only the spectrum array item container
(saic) and the spectrum metadata item container (smic) of the specfiles
"specfile_name_1" and "specfile_name_3" are saved. ::

    msrunContainer.save(specfiles=["specfile_name_1", "specfile_name_3"],
                        sai=True, smi=True
                        )


**Loading an MsrunContainer from the hard disk**

Before loading an ``MsrunContainer`` from the hard disk, a specfile entry has to
be added to its ``.info`` attribute. This can be done by calling
:func:`.addSpecfile() <maspy.core.MsrunContainer.addSpecfile>` with the name of
the specfile and the path to the filedirectory. Afterwards the files can be
loaded by calling :func:`.load() <maspy.core.MsrunContainer.load>`, which will
import all specfiles present in ``.info`` and update the ``status`` variable of
``.info``. ::

    >>> msrunContainer = maspy.core.MsrunContainer()
    >>> msrunContainer.addSpecfile('specfile_name_1', 'filedirectory')
    >>> msrunContainer.info
    {u'specfile_name_1': {u'path': u'filedirectory',
                          u'status': {u'ci': False,
                                      u'rm': False,
                                      u'sai': False,
                                      u'si': False,
                                      u'smi': False}}}
    >>> msrunContainer.load()
    >>> msrunContainer.info
    {u'specfile_name_1': {u'path': u'filedirectory',
                          u'status': {u'ci': True,
                                      u'rm': True,
                                      u'sai': True,
                                      u'si': True,
                                      u'smi': True}}}

Similar to saving only parts of an ``MsrunContainer`` it is also possible to
only select a subset of specfiles present in ``.info`` and specify which data
types are imported. ::

    >>> msrunContainer = maspy.core.MsrunContainer()
    >>> msrunContainer.addSpecfile('specfile_name_1', 'filedirectory')
    >>> msrunContainer.info
    {u'specfile_name_1': {u'path': u'filedirectory',
                          u'status': {u'ci': False,
                                      u'rm': False,
                                      u'sai': False,
                                      u'si': False,
                                      u'smi': False}}}
    >>> msrunContainer.load(specfiles='specfile_name_1', sai=True, smi=True)
    >>> msrunContainer.info
    {u'specfile_name_1': {u'path': u'filedirectory',
                          u'status': {u'ci': False,
                                      u'rm': False,
                                      u'sai': True,
                                      u'si': False,
                                      u'smi': True}}}


**Deleting data from an MsrunContainer**

If specific data types are not needed anymore, they can be removed to free
memory. This can be done by using :func:`.removeData()
<maspy.core.MsrunContainer.removeData>` and parsing arguments to specify
specfiles and which data types to remove. It is recommended to use this method
to remove data as it automatically updates the ``.info`` attribute of the
``MsrunContainer``. The following command removes the ``Sai`` and ``Smi`` items
of the specfile "specfile_name_1". ::

    >>> msrunContainer.info
    {u'specfile_name_1': {u'path': u'filedirectory',
                          u'status': {u'ci': True,
                                      u'rm': True,
                                      u'sai': True,
                                      u'si': True,
                                      u'smi': True}}}
    >>> msrunContainer.removeData('specfile_name_1', sai=True, smi=True)
    >>> msrunContainer.info
    {u'specfile_name_1': {u'path': u'filedirectory',
                          u'status': {u'ci': True,
                                      u'rm': True,
                                      u'sai': False,
                                      u'si': True,
                                      u'smi': False}}}

A specfile can be completely removed from an ``MsrunContainer`` by calling
:func:`.removeSpecfile() <maspy.core.MsrunContainer.removeSpecfile>`, which
deletes all data from the containers and in addition the entry from the
``.info`` attribute. ::

    msrunContainer.removeSpecfile('specfile_name_1')


**Exporting specfiles from MsrunContainer to mzML files.**

""" Show how to generate a new mzML file, explanations about requirements """ ::

    import maspy.writer
    maspy.writer.writeMzml('specfile_name_1', msrunContainer, 'D:/maspy_test')


**Accessing data from MsrunContainer.**




Introduction:
spectra are the most basic results produced by any mass spectrometry experiment
-> start of the tutorial / introduction

- what is mzML ?
    eg from pyteomics: mzML is an XML-based format for experimental data obtained on MS/MS or LC-MS setups.
    or from http://www.psidev.info/mzml: format for encoding raw spectrometer output
- how to convert vendor format files to mzML?
    - use msConvert from proteoWizard suite
    - also allows peak picking to obtain centroided data

- what exactly is the problem of mzML files?
    - modern mass spectrometers can generate tens of thousands spectra per hour
      resulting in huge mzML files because of spectral information
    - parsing takes a lot of time even tough sometimes only part of the data
      has to be accessed

Maspy implementation:
- maspy features an own internal representation of mzML files
    - contains all information present in the mzML file, split into four groups
    - run metadata, spectrum item metadata, spectrum array item, chromatogram item
    - explain all of them ...
    - explain what spectrum items are, and what their purpose is
- after conversion, each group is saved seperately and can therefore be imported individually
    - run metadata in xml, arrays as binary data, all other data as JSON
    - smaller file size
    - faster data access
- maspy allows writing a new mzML file from the msrunContainer (the internal mzMl represenation)
    - allows to modify mzML within maspy and pass the changed mzML file to external software
    - note: indexing not yet supported



--- ad peak picking, from the internet ---
What is the difference between Profile and Centroid MS data?

MS data collected off an instrument is presented as either profile or
centroid mode. Shown below are two mass spectra illustrating an ion
cluster for profile data and a centroid mass spectrum created from the
profile data.

In profile mode, a peak is represented by a collection of signals over
several scans. The advantage of profile data is it is easier to
classify a signal as a true peak from noise off the instrument.

In centroid mode, the signals are displayed as discrete m/z with zero
line widths. The advantage of centroid data is the file size is
significantly smaller as there is less information describing a
signal.  --- ad peak picking ---
