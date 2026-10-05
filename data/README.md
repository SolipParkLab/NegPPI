# Input data

This repository does not redistribute the large external datasets used in the analysis.
Users should download the following files from their original sources.

## CPTAC pan-cancer proteome

File:

`Proteome_Broad_Institute_harmonized_v1.tsv`

Source:

CPTAC Pan-Cancer Proteogenomic Data Portal  
https://proteomic.datacommons.cancer.gov/pdc/cptac-pancancer

The file contains harmonized tumor proteome abundance measurements used for
pan-cancer and mutation-associated co-abundance analyses.

## CPTAC somatic mutation data

File:

`PanCan_Union_Maf_Broad_WashU_v1.1.maf`

Source:

CPTAC Pan-Cancer Proteogenomic Data Portal  
https://proteomic.datacommons.cancer.gov/pdc/cptac-pancancer

This file is distributed as part of the Broad/WashU union mutation dataset.

## AlphaMissense predictions

File:

`AlphaMissense_aa_substitutions.tsv`

Source:

AlphaMissense

The file contains amino-acid substitution pathogenicity predictions used to classify
missense variants.

## RefSeq-to-HGNC mapping

File included in this repository:

`data/reference/p2s.rds`

This object contains two columns:

- `refseq_peptide`
- `hgnc_symbol`

The mapping was generated using BioMart and was used to map RefSeq peptide
identifiers in the CPTAC proteome data to HGNC gene symbols.
