Peptide features in MasPy - basics
----------------------------------

**Introduction**
    - what are peptide featuers and what is their purpose
        - MS1 mass trace over time as a molecule is eluting from the HPLC
        - very often the isotopes are incorportaed to get a more reliable feature and used to calculate the charge state
        - minimal information: rt boundaries, rt center, m/z monoisotopic peak, charge state, intensity area
    - FeatureItem and FeatureItemContainer
    - Two algorithms used for MasPy: openMS feature finder centroided, Dinosaur (recently published)
      both output formats are supported by MasPy
    - Import of featureXML and dinosaur tsv
    - matching Sii or Si to features
