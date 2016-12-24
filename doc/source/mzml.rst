MS spectra in MasPy
-------------------

The mzML file format
^^^^^^^^^^^^^^^^^^^^

Every vendor software produces mass spectrometer output files in a different
proprietary format. It is a difficult and time consuming task for software
developers to support all of these different formats and format versions.
Therefore the file format mzML has been developed by the Proteomics Standards
Initiative (PSI) as the community standard for representation of mass
spectrometry results. mzML is an open, XML- based format that not only allows to
store recorded mass spectrum information but also metadata of the instrument
configuration, acquisition settings, software used for data processing and
sample descriptions. Ultimately, it is desirable to universally use mzML for
archiving, sharing, and processing of mass spectrometry data and thus for all
software to support and use the mzML format.

.. note::
    Refer to `www.psidev.info <http://www.psidev.info/index.php?q=node/257>`_
    for details on the XML schema definition and mzML file specifications, see
    also the publication `Mass Spectrometer Output File Format mzML
    <http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3073315>`_)

.. note::
    We recommend using ProteoWizard for conversion of vendor format files to
    mzML. The software can be downloaded from their `website
    <http://proteowizard.sourceforge.net>`_, a detailed protocol how to use
    ProteoWizard can be found `here
    <http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4113728>`_.

    The raw spectral data recorded by an instrument can be either stored as
    profile or centroid data. Meassured mass spectra are initially recorded in
    profile mode, where each mass peak is represented by a number of m/z and
    intensity values describing a peak shape. In centroid mode this information
    is reduced to the centroid of the peak shape, storing only one single pair
    of a dinstinct m/z value and an intensity. The process of converting profile
    data to centroid data is called peak picking and can be applied as a filter
    while converting vendor format files to mzML files using ProteoWizard, see
    the ProteoWizard protocol. The representation as centroid data is easier to
    work with, saves memory and is sufficient for most applications. Therefore
    we recommend the utilization of centroid data for MasPy.


MsrunContainer
^^^^^^^^^^^^^^

Modern mass spectrometers can generate tens of thousands of spectra per hour
resulting in huge mzML files. Opening and parsing such large XML files takes a
lot of time. MzML files can contain a byte-offset index which allows directly
reading certain spectra without parsing the whole file. This can increase
performance when only one or a few specific spectra have to be accessed at a
time.

The actual spectral information takes up to largest part of a typical mzML file.
However, sometimes only a certain type of information needs to be accessed, for
example the spectrum metadata. Therefore we split the information that is
contained in mzML files into four data groups; run metadata (``Rm``), spectrum
metadata items (``Smi``), spectrum array items (``Sai``) and chromatogram items
(``Ci``). Each of these data groups is stored separately in MasPy and has its
own file type, thus it can be accessed, saved and loaded independently of the
others. All four data types are stored in the MasPy class :class:`MsrunContainer
<maspy.core.MsrunContainer>`. Altough the data is split into multiple parts, all
information originally contained in an mzML file is still present. This allows
the conversion from MsrunContainer to mzML at any given time. #TODO: Why do we
want to be able to export mzML files? (Preffered data format for archiving and
sharing data and to use as input for other software packages)

See tutorial/docstrings xxx for details on the MsrunContainer file
format. #TODO:

Fig.: MsrunContainer #TODO: make figure

* run metadata
* spectrum metadata items
* spectrum array items
* chromatogram items
* spectrum items


Run metadata (:class:`Rm`)
^^^^^^^^^^^^^^^^^^^^^^^^^^

The run metadata element contains all information of an mzML file, which is not
directly part of the acquired spectra and chromatograms. This covers, amongst
others, a description of the instrument configuration, a list of software used
for data processing and a list of applied data processing steps. In addition it
is possible to add contact information and a description of the analyzed samples
to the mzML file. In MasPy all of these mzML elements are converted to an
``lxml.etree.Element`` and stored in ``MsrunContainer.rmc`` (Rm container).

.. note::
    Software which is used to process data of an mzML file should be listed in
    the mzML element "softwareList", and all applied data processing steps
    should be documented in the "dataProcessingList" element.


Spectrum array item (:class:`Sai <maspy.core.Sai>`), spectrum metadata item (:class:`Smi <maspy.core.Smi>`)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

An mzML spectrum element contains all information of an acquired MS spectrum,
including numerical arrays containing at least recorded m/z and intensity values
of the observed ions but also plenty of metadata describing for example details
of the acquisition like base peak m/z and intensity, scan start time, ms level,
MS2 isolation window or precursor information of MS2 scans. In MasPy this
information is split into a metadata containing part and the spectrum array data
and put into two separate data structures; spectrum metadata item (``Smi``) and
spectrum array item (``Sai``), respectively. ``Smi`` elements are stored in
``MsrunContainer.smic`` (Smi container) and ``Sai`` elements in
``MsrunContainer.saic`` (Sai container). In order to recreate an mzML spectrum
element the information of both MasPy data types (``Smi`` and ``Sai``) is
necessary.


Chromatogram item (:class:`Ci <maspy.core.Ci>`)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

An mzML chromatogram element is similar to a spectrum element, containing
metadata and numerical arrays. Common chromatogram types are ``total ion current
chromatogram``, ``selected ion current chromatogram`` and ``basepeak
chromatogram``. All of them contain time and intensity data points, however,
other chromatogram types can also contain absorption or emission values instead
of intensities. In the current MasPy implementation chromatogram elements are
not split into two data types but the metadata and array information is put into
one single data structure called chromatogram item (``Ci``), which is stored in
``MsrunContainer.cic`` (Ci container).


Spectrum item (:class:`Si <maspy.core.Si>`)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The mzML file  serves as a data container for active data processing but also
for data sharing and archiving. Thus the spectrum elements contain a lot of
metadata information not needed for most data analysis applications. In addition
all information stored in spectrum elements have to be in accordance with the
mzML xml scheme definition and the Controlled Vocabularies (CV's) of the PSI,
`see <http://www.psidev.info/groups /controlled- vocabularies>`_. Altough in
principle this standardization is beneficial and perfectly reasonable, when
actively working with the data it is not always required and can make things
unnecessarily complicated.

To circumvent this problem MasPy provides a simpler data type for working with
spectrum metadata, called spectrum item (:class:`Si <maspy.core.Si>`). The
``Si`` class has a flat structure, meaning that attributes are not nested inside
other elements but are stored directly as attributes of the class. ``Si``
attributes can be manipulated without restrictions and new attributes can simply
be added. Specific functions can be used to selectively extract information from
``Smi``. This allows import only the currently needed spectrum metadata
attributes, like retention time, ms level or MS2 precursor information, thereby
making the ``Si`` more memory efficient. In order to make lasting changes to the
mzML file ``Si`` attributes have to be translated to the respective ``Smi``
elements. These changes however have to strictly follow the mzML specifications
and syntax. Thus it is recommend to use existing functions or implement new ones
that make changes to ``Smi`` elements in a controlled manner.

Each spectrum present in an mzML file is therefore represented threefold in
MasPy. First the ``Smi`` contains a complete representation of all metadata
information present in an mzML spectrum element. However, this data type is not
intended to be used for standard data analysis and will normally only be
accessed to make lasting, documented changes to spectrum metadata and for
generating new mzML files. Second the ``Sai`` contains the actual ion
information recorded by the mass spectrometer. This data type will be used
whenever the ion spectra have to be analyzed or manipulated. In addition it is
also required for generating new mzML files. And third the ``Si``, which can be
considered as the spectrum metadata workspace in MasPy, allowing convenient
access to metadata and simple processing of this data without directly altering
the original mzML information. This data type will be used for most data
processing and analysis steps in MasPy.


MsrunContainer.info
^^^^^^^^^^^^^^^^^^^

*MsrunContainer.info -> which specfiles are present, what is the current path
(used for loading or saving) , which data types are currently imported*


MasPy file formats
^^^^^^^^^^^^^^^^^^

*This section will contain information about how the data contained in an
MsrunContainer is written to the hard drive. (one file type per data type:
mrc_rm, mrc_si, mrc_sai, mrc_smi, mrc_ci)*


Basic code examples
^^^^^^^^^^^^^^^^^^^

Importing an mzML file
""""""""""""""""""""""

mzML files can be imported by using the function
:func:`maspy.reader.importMzml()`, the imported specfile is then added to the
``MsrunContainer`` instance passed to the function. ::

    import maspy.core
    import maspy.reader

    mzmlfilepath = 'filedirectory/specfile_name_1.mzML'
    msrunContainer = maspy.core.MsrunContainer()
    maspy.reader.importMzml(mzmlfilepath, msrunContainer)


Saving an MsrunContainer to the hard disk
"""""""""""""""""""""""""""""""""""""""""

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
of "specfiles" can either be the name of one single specfile or a list of
specfile names. In the following example only the spectrum array item container
(saic) and the spectrum metadata item container (smic) of the specfiles
"specfile_name_1" and "specfile_name_3" are saved. ::

    msrunContainer.save(specfiles=["specfile_name_1", "specfile_name_3"],
                        sai=True, smi=True
                        )


Loading an MsrunContainer from the hard disk
""""""""""""""""""""""""""""""""""""""""""""

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


Deleting data from an MsrunContainer
""""""""""""""""""""""""""""""""""""

If specific data types are not needed anymore, they can be removed to free
memory. This can be done by using :func:`.removeData()
<maspy.core.MsrunContainer.removeData>` and parsing arguments to specify
specfiles and which data types to remove. It is recommended to always use this
method to remove data instead of manually deleting container entries, because
using ``.removeData`` automatically updates the ``.info`` attribute of the
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


Exporting specfiles from MsrunContainer to mzML files
"""""""""""""""""""""""""""""""""""""""""""""""""""""

After working in MasPy it might be desirable to export the MsrunContainer back
into an mzML file which can be used as input for another software or simply for
archiving and sharing mass spectrometry data. An mzML file is generated by using
the function :func:`maspy.writer.writeMzml()` and passing at least the
``specfile`` name that should be exported, an ``MsrunContainer`` and the
``output directory``. In order to write a valid and complete mzML file all data
types except for ``Si`` have to be present in the ``MsrunContainer``. ::

    import maspy.writer
    maspy.writer.writeMzml('specfile_name_1', msrunContainer, '/filedirectory')

.. note::
    Optionally it is possible to supply a list of ``spectrumIds`` and
    ``chromatogramIds`` to only select a subset of spectra and chromatograms
    that should be written to the mzML file. The supplied lists of element ids
    have to be sorted in the order they should be written to the mzML file.


Accessing data from an MsrunContainer
"""""""""""""""""""""""""""""""""""""

#TODO: *examples of .getItem, .getArrays, ... *

