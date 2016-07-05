Peptide LC-MS features in MasPy
-------------------------------

Chromatographic separation coupled directly to the mass spectrometer results in
analytes to appear over time as they emerge from the chromatography column in a
more or less Gaussian peak shape. The elution profile of an analyte can be
recapitulated by using the spectral information present in the MS1 scans. The
simplest way to do this is by extracting the intensities of an ion species with
a given m/z value in consecutive MS1 scans as long as the ion is detectable.
Combining the extracted intensities with the respective MS1 retention times
results in a so called extracted ion chromatogram (EIC or XIC). The intensity
area obtained by integrating the XIC is frequently used as a measure of
abundance in label free quantification (LFQ) and stable isotope labeling (SIL)
workflows.

In MS spectra each peptide consists of an isotope envelope of multiple ion
species with different m/z values. Combining the XICs of the different isotope
states of the same analyte allows the inference of its charge state and
therefore its mass. In addition the availability of more information results in
more accurate intensity area estimates and thus increased accuracy for
quantification. The combined information of XICs and isotope clusters can be
referred to as a peptide LC-MS feature or more commonly simply as a feature.


Representation of LC-MS features in MasPy
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Peptide LC-MS features, but also XICs, are represented in MasPy with the feature
item class :class:`maspy.core.Fi`. Its structure is kept very simple, similar to
:class:`Si <maspy.core.Si>` and :class:`Sii <maspy.core.Sii>`, with only a few
mandatory attributes. Each instance is uniquely identified by the combination of
the ``Fi.id`` and ``Fi.specfile`` attributes. However, the ``id`` attribute is
not associated with any particular MS scan, which was the case for ``Si`` and
``Sii``. Further attributes that should always be supplied when importing
features into MasPy are ``mz``, ``rt``, ``rtLow``, ``rtHigh`` and ``intensity``.
Altough not absolutely mandatory the ``charge`` attribute should also be
supplied whenever possible, since the charge information is used in some
algorithms.

The :class:`FiContainer <maspy.core.FiContainer>` is used to store feature items
of one or multiple specfiles. The container allows saving and loading of
imported results and provides methods for convenient access to the data.

Attribute naming conventions in MasPy and additional attributes that might
be necessary for working with feature items:

#TODO: maybe change mz to obsMz to be consistent between data types

    - ``mz`` the experimentally observed mass to charge ratio (Dalton /
      charge). Normally the m/z value of the monoisotopic peak.

    - ``rt`` the retention time center of the feature.

    - ``rtLow`` the lower retention time boundary of the feature.

    - ``rtHigh`` the upper retention time boundary of the feature.

    - ``intensity`` an estimator for the feature abundance. The preferred value
      is the integrated intensity area, but the feature apex intensity is also
      possible.

    - ``charge`` the charge state of the feature.

    - ``peptide`` the peptide sequence of the :class:`Sii <maspy.core.Sii>`
      that is used for annotating the feature.

    - ``sequence`` the plain amino acid sequence of the
      :class:`Sii <maspy.core.Sii>` that is used for annotating the feature.

    - ``score`` or any other score attribute name of the
      :class:`Sii <maspy.core.Sii>` that is used for annotating the feature.
      It describes the quality of a spectrum identifications.

    - ``obsMz`` the experimentally observed mass to charge ratio of the feature
      (Dalton / charge). Normally the m/z value of the monoisotopic peak.

    - ``obsMh`` the experimentally observed mass to charge ratio of the
      feature, calculated for the mono protonated ion (Dalton / charge).
      Normally the monoisotopic peak.

    - ``obsMass`` the experimentally observed not protonated mass of a feature
      calculated by using the mz and charge values (Dalton / charge).
      Normally the monoisotopic mass.

    - ``excMz`` the exact calculated mass to charge ratio of the peptide
      (Dalton / charge). Normally the monoisotopic ion.

    - ``excMh`` the exact calculated mass to charge ratio of the peptide,
      calculated for the mono protonated state (Dalton / charge). Normally the
      monoisotopic ion.

    - ``excMass`` the exact calculated mass of the not protonated peptide
      (Dalton / charge). Normally the monoisotopic mass.


MasPy internal feature item attributes:

    - ``isValid`` can be used to flag if a Fi has passed a given quality
      threshold.

    - ``isMatched`` can be used to flag if a Fi has been matched to any
      :class:`Si <maspy.core.Si>` or :class:`Sii <maspy.core.Sii>` elements.

    - ``isAnnotated`` can be used to flag if a Fi has been annotated with a
      :class:`Sii <maspy.core.Sii>` element and therefore with an identified
      peptide sequence.

    - ``siIds`` a list of :class:`Si <maspy.core.Si>` elements that have been
      matched to the feature item.

    - ``siiIds`` a list of :class:`Sii <maspy.core.Sii>` elements that have
      been matched to the feature item.



Supported feature detection algorithms
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Currently MasPy supports the import of two feature containing file types; the
openMS feature file format ``.featureXML`` and the ``.feature.tsv`` format
generated by the open source tool Dinosaur. However, adding import routines for
additional file formats should be trivial an can be done on demand.

The FeatureFinderCentroided node from openMS is one of the best established open
source LC-MS feature defining algorithms. It can be used independently of a data
analysis pipeline and other processing steps. It was published in 2013 as part
of a complete openMS pipeline: `An Automated Pipeline for High- Throughput
Label-Free Quantitative Proteomics
<http://pubs.acs.org/doi/abs/10.1021/pr300992u>`_. Since its publication it was
applied in numerous publications and has been reused in at least two additional
open source projects: `DeMix
<http://www.mcponline.org/content/13/11/3211.long>`_ and `DeMix-Q
<http://www.mcponline.org/content/15/4/1467.long>`_.

`Dinosaur: A Refined Open- Source Peptide MS Feature Detector
<http://pubs.acs.org/doi/abs/10.1021/acs.jproteome.6b00016>`_ published in 2016,
is an algorithm based on the graph model concept for feature detection
introduced by MaxQuant in 2008. Dinosaur seems to provide similar or better
results then the FeatureFinderCentroided node of openMS with a substantial
increase in runtime performance. It is available on `Github
<https://github.com/fickludd/dinosaur>`_.


Basic code examples
^^^^^^^^^^^^^^^^^^^

**Importing peptide features**

The function :func:`maspy.reader.importPeptideFeatures()` is used to import LC-
MS features from a file. It automatically recognises the file type by the file
name extension and executes the respective import routine. Therefore the file
extension has to be either ``.featurexml`` (openMS) or ``.feature.tsv``
(Dinosaur) and is not case sensitive. The imported feature items are stored in
the ``FiContainer`` instance passed to the function. ::

    import maspy.core
    import maspy.reader

    fiContainer = maspy.core.FiContainer()
    maspy.reader.importPeptideFeatures(fiContainer, 'filelocation/f.featureXML',
                                       'specfile_name_1')

**Matching spectrum identification items to feature items**

The peptide underlying a LC-MS feature can be determined by using the
information of identified MSn scans. In MasPy this can be achieved by using
:func:`maspy.featuremethods.matchToFeatures()`, which allows matching ``Sii`` to
``Fi`` elements by comparing their m/z, retention time and charge information.
User defined tolerance values for matching should be passed to the function, for
details see the docstring documentation. However, the default settings should be
appropriate for typical high resolution MS1 data as obtained by Thermo Orbitrap
instruments.

#TODO: describe the print output

    >>> import maspy.featuremethods
    >>> maspy.featuremethods.matchToFeatures(fiContainer, siiContainer,
    >>>                                      specfiles='specfile_name_1')
    ------ specfile_name_1 ------
    Annotated features:                      3802 / 20437 = 18.6 %
    Spectra matched to features:             4240 / 4898 = 86.6 %

.. note::

    #TODO: describe which attributes must be present in the Sii items and link
    to the tutorial that describes how to obtain these attributes.
    #charge, m/z, rentention time

**Accessing data stored in a FiContainer**

#TODO: describe .getItem(), .getArrays()
