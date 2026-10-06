# Input data

This repository does not redistribute the large external datasets used in the
analyses. Users should obtain the required files from their original sources
and place them under `data/raw/` using the directory structure shown below.

The RefSeq-to-HGNC mapping used by the analyses is included under
`data/reference/`.

## Expected directory structure

```text
data/
├── README.md
├── reference/
│   └── p2s.rds
└── raw/
    ├── AlphaMissense_aa_substitutions.tsv
    ├── PanCan_Union_Maf_Broad_WashU_v1.1.maf
    ├── Proteome_Broad_Institute_harmonized_v1.tsv
    ├── Li_2023/
    │   └── mmc7.xlsx
    └── normal_tumor/
        ├── BRCA/
        │   ├── CPTAC2_Breast_Prospective_Collection_BI_Proteome.tmt10.tsv
        │   └── S039_Breast_Cancer_Prospective_Collection_Specimens_r1.xlsx
        ├── ccRCC/
        │   ├── CPTAC3_Clear_Cell_Renal_Cell_Carcinoma_Proteome.tmt10.tsv
        │   └── S044_CPTAC_CCRCC_Discovery_Cohort_Specimens_r1_Sept2018.xlsx
        ├── COAD/
        │   ├── CPTAC2_Colon_Prospective_Collection_PNNL_Proteome.tmt10.tsv
        │   └── CPTAC_Colon_Cancer_Prospective_Collection_Biospecimens_S037.xlsx
        ├── LSCC/
        │   ├── CPTAC3_Lung_Squamous_Cell_Carcinoma_Proteome.tmt11.tsv
        │   └── S063_S058updated_BI_CPTAC3_LSCC_Tumor_Normal_Mapping_r2.xlsx
        ├── LUAD/
        │   ├── CPTAC3_Lung_Adeno_Carcinoma_Proteome.tmt10.tsv
        │   └── S046_S056_BI_CPTAC3_LUAD_Discovery_Cohort_Samples_r2_July2020.xlsx
        ├── PDAC/
        │   ├── CPTAC3_Pancreatic_Ductal_Adenocarcinoma_Proteome.tmt11.tsv
        │   └── S061_CPTAC_PDA_Discovery_Cohort_Specimens_r1_Feb2021.xlsx
        └── UCEC/
            ├── CPTAC3_Uterine_Corpus_Endometrial_Carcinoma_Proteome.tmt10.tsv
            └── S043_CPTAC_UCEC_Discovery_Cohort_Study_Specimens_r1_Sept2018.xlsx
```

Raw input files under `data/raw/` are excluded from version control.

## 1. CPTAC pan-cancer harmonized proteome

File:

`Proteome_Broad_Institute_harmonized_v1.tsv`

Source:

CPTAC Pan-Cancer Proteogenomic Data Portal  
https://proteomic.datacommons.cancer.gov/pdc/cptac-pancancer

This tumor-only harmonized proteome dataset is used for mutation-stratified
co-abundance analyses, including the pan-cancer analysis in
`R/01_mutation_stratified_coabundance.R` and the cancer-type-specific
reanalysis in `R/08_individual_cancer_mutation_rewiring.R`.

Protein abundance values are log2-transformed relative abundance measurements
generated using isobaric labeling with a common pooled reference.

## 2. CPTAC somatic mutation data

File:

`PanCan_Union_Maf_Broad_WashU_v1.1.maf`

Distributed as part of:

`Mutation_Broad_WashU_union_v1.zip`

Source:

CPTAC Pan-Cancer Proteogenomic Data Portal  
https://proteomic.datacommons.cancer.gov/pdc/cptac-pancancer

The same pan-cancer somatic mutation file is used for both the harmonized
tumor analyses and the cancer-type-specific analyses.

For mutation-stratified analyses, samples with any somatic mutation in the
query gene are excluded from the query-specific wild-type group.

The MAF cohort label for clear cell renal cell carcinoma is `CCRCC`, whereas
`ccRCC` is used as the display label in this repository.

## 3. AlphaMissense predictions

File:

`AlphaMissense_aa_substitutions.tsv`

Source:

AlphaMissense (Cheng et al., Science, 2023)

AlphaMissense amino-acid substitution predictions are used to annotate somatic
missense mutations in mutation-stratified analyses.

The original AlphaMissense classifications include benign, ambiguous, and
pathogenic predictions. Downstream analyses in this study compare benign and
pathogenic variants.

## 4. Cancer-type-specific adjacent-normal and tumor proteomes

Normal-tumor co-abundance analyses use the original cohort-specific CPTAC
proteome datasets and specimen annotations rather than the harmonized
tumor-only pan-cancer proteome.

Seven cancer types with available adjacent normal proteomes were analyzed:
BRCA, ccRCC, COAD, LSCC, LUAD, PDAC, and UCEC.

### BRCA

PDC accession: `PDC000120`

Proteome:

`CPTAC2_Breast_Prospective_Collection_BI_Proteome.tmt10.tsv`

Specimen annotation:

`S039_Breast_Cancer_Prospective_Collection_Specimens_r1.xlsx`

Expected samples after specimen annotation and matching to the mutation dataset:

- Adjacent normal: 18
- Tumor: 113

### ccRCC

PDC accession: `PDC000127`

Proteome:

`CPTAC3_Clear_Cell_Renal_Cell_Carcinoma_Proteome.tmt10.tsv`

Specimen annotation:

`S044_CPTAC_CCRCC_Discovery_Cohort_Specimens_r1_Sept2018.xlsx`

Expected samples:

- Adjacent normal: 84
- Tumor: 110

### COAD

PDC accession: `PDC000116`

Proteome:

`CPTAC2_Colon_Prospective_Collection_PNNL_Proteome.tmt10.tsv`

Specimen annotation:

`CPTAC_Colon_Cancer_Prospective_Collection_Biospecimens_S037.xlsx`

Expected samples:

- Adjacent normal: 95
- Tumor: 94

### LSCC

PDC accession: `PDC000234`

Proteome:

`CPTAC3_Lung_Squamous_Cell_Carcinoma_Proteome.tmt11.tsv`

Specimen annotation:

`S063_S058updated_BI_CPTAC3_LSCC_Tumor_Normal_Mapping_r2.xlsx`

Expected samples:

- Adjacent normal: 100
- Tumor: 108

### LUAD

PDC accession: `PDC000153`

Proteome:

`CPTAC3_Lung_Adeno_Carcinoma_Proteome.tmt10.tsv`

Specimen annotation:

`S046_S056_BI_CPTAC3_LUAD_Discovery_Cohort_Samples_r2_July2020.xlsx`

Expected samples:

- Adjacent normal: 99
- Tumor: 106

### PDAC

PDC accession: `PDC000270`

Proteome:

`CPTAC3_Pancreatic_Ductal_Adenocarcinoma_Proteome.tmt11.tsv`

Specimen annotation:

`S061_CPTAC_PDA_Discovery_Cohort_Specimens_r1_Feb2021.xlsx`

Expected samples:

- Adjacent normal: 66
- Tumor: 137

### UCEC

PDC accession: `PDC000125`

Proteome:

`CPTAC3_Uterine_Corpus_Endometrial_Carcinoma_Proteome.tmt10.tsv`

Specimen annotation:

`S043_CPTAC_UCEC_Discovery_Cohort_Study_Specimens_r1_Sept2018.xlsx`

Expected samples:

- Adjacent normal: 30
- Tumor: 95

The cohort-specific proteome files contain several quantitative columns.
The analysis retains columns ending in `Log Ratio` and excludes
`Unshared Log Ratio` measurements.

`R/02_prepare_normal_tumor_proteomes.R` performs cohort-specific specimen
annotation, participant-ID conversion, and matching to the pan-cancer mutation
dataset.

## 5. RefSeq-to-HGNC mapping

File included in this repository:

`data/reference/p2s.rds`

The R object contains two columns:

- `refseq_peptide`
- `hgnc_symbol`

This mapping was generated using BioMart and maps RefSeq peptide identifiers
to HGNC gene symbols. It is used to map proteomic RefSeq peptide identifiers
to gene symbols in the mutation-stratified analyses.

Example:

```text
refseq_peptide   hgnc_symbol
YP_003024026     MT-ND1
YP_003024027     MT-ND2
YP_003024028     MT-CO1
```

## 6. Cancer hallmark annotations

File:

`data/raw/Li_2023/mmc7.xlsx`

The analysis uses hallmark annotations from Supplementary Table S7B
(`Table S7B`) of Li et al. (2023).

The first column is interpreted as the gene identifier, and hallmark membership
is defined by entries marked `Y` in the corresponding annotation columns.

These annotations are used by:

`R/06_hallmark_rewiring.R`

to calculate module-level hallmark loadings and project inter-module rewiring
onto pairs of cancer hallmarks.

## Analysis workflow and input dependencies

The main analysis workflow is:

1. `R/01_mutation_stratified_coabundance.R`  
   Uses the harmonized pan-cancer tumor proteome, somatic mutation calls,
   AlphaMissense predictions, and the RefSeq-to-HGNC mapping to calculate
   mutation-stratified query-target co-abundance associations.

2. `R/02_prepare_normal_tumor_proteomes.R`  
   Prepares cancer-type-specific adjacent-normal and tumor proteome matrices
   from the original CPTAC cohort datasets and specimen annotations.

3. `R/03_normal_tumor_coabundance.R`  
   Calculates query-specific wild-type co-abundance correlations separately in
   adjacent normal and tumor tissues.

4. `R/04_module_architecture.R`  
   Defines normal-tissue co-abundance modules and projects significant
   normal and tumor associations onto the normal-derived module structure.

5. `R/05_module_pair_rewiring.R`  
   Quantifies changes in inter-module co-abundance associations between
   adjacent normal and tumor tissues.

6. `R/06_hallmark_rewiring.R`  
   Uses the Li et al. hallmark annotations to quantify hallmark-pair rewiring
   across the normal-derived module architecture.

7. `R/07_pancancer_mutation_rewiring.R`  
   Classifies pan-cancer mutation-associated strengthening and weakening from
   the mutation-stratified results generated by script 01.

8. `R/08_individual_cancer_mutation_rewiring.R`  
   Recalculates mutation-stratified co-abundance independently within each
   cancer type and generates cancer-specific mutation-associated rewiring
   results.

9. `R/09_mutation_rewiring_by_association_origin.R`  
   Integrates the cancer-specific mutation-rewiring results from script 08
   with the normal-tumor association histories from script 03 to compare
   Pre-existing and Tumor-acquired associations.

10. `R/10_pancancer_query_centric_association_change.R`  
    Performs a pan-cancer query-centric analysis of strengthening versus
    weakening directions. The WT baseline association is required to meet the
    correlation/FDR threshold, whereas no significance threshold is imposed on
    the WT-to-MUT association change itself.
