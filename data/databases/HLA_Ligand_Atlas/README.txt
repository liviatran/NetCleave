The attached files were compiled based on a query against the HLA Ligand Atlas
(hla-ligand-atlas.org). The following files are included:

[LICENSE.txt]
A copy of the CC-BY 4.0 license by which applies to the data provided by the
HLA ligand atlas.

[README.txt]
This document

[AUTHORS.txt]
The list of authors that contributed to the HLA Ligand Atlas

[HLA_aggregated.tsv]
A TSV file with one row per unique peptide sequence comprising the following
columns:
    - peptide_sequence_id
      A unique numeric identifier for peptide sequences
    - peptide_sequence
      The peptide amino acid sequence
    - hla_class
      Indicates whether the peptide was observed in HLA-I, HLA-II or both
      types of experiments (HLA-I+II)
    - donor_alleles
      A comma separated list of all donor alleles of donors in which this
      peptide was observed. The format used is
                [swn]/{gene}*[digits]:[digits]
      where the first character indicates whether the peptides is predicted as
      a [s]trong, [w]eak or [n]on-binder against this allele based on our
      prediction workflow (https://hla-ligand-atlas.org/docs/computational).
    - tissues
      A comma separated list of tissues in which the peptide was observed

Note that from this aggregated data one can not necessarily deduce whether a
certain tissue - HLA allele or tissue - HLA class combination was observed.

[HLA_sample_hits.tsv]
A TSV file with one row for each occurrence of a peptide in a sample
measurement comprising the following columns:
    - peptide_sequence_id
      A unique numeric identifier for peptide sequences, same as above
    - donor
      The donor the sample originated from
    - tissue
      The tissue the sample originated from
    - hla_class
      The HLA class investigated in the measurement

[HLA_donors.tsv]
A TSV file mapping the donors to their HLA alleles

[HLA_protein_map]
A TSV file mapping the peptide sequences (by id, see above) to their source
proteins (by uniprot id)
