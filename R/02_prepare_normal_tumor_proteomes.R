# 02_prepare_normal_tumor_proteomes.R
#
# Prepare cohort-specific adjacent-normal and tumor proteome matrices
# for the seven CPTAC cancer types used in normal-tumor analyses.
#
# This script reproduces the cohort-specific sample annotation and
# sample-matching procedures used in the original analysis.
#
# Output:
#   results/normal_tumor_proteomes.rds
#   results/normal_tumor_sample_counts.tsv


# -------------------------------------------------------------------------
# 1. Packages
# -------------------------------------------------------------------------

library(readxl)
library(data.table)
library(stringr)


# -------------------------------------------------------------------------
# 2. Paths
# -------------------------------------------------------------------------

maf_file <- file.path(
  "data",
  "raw",
  "PanCan_Union_Maf_Broad_WashU_v1.1.maf"
)

normal_tumor_dir <- file.path(
  "data",
  "raw",
  "normal_tumor"
)

output_file <- file.path(
  "results",
  "normal_tumor_proteomes.rds"
)

sample_count_file <- file.path(
  "results",
  "normal_tumor_sample_counts.tsv"
)

dir.create(
  "results",
  recursive = TRUE,
  showWarnings = FALSE
)


# -------------------------------------------------------------------------
# 3. Helper functions
# -------------------------------------------------------------------------

# Clean cohort-specific CPTAC proteome files.
#
# Only "Log Ratio" measurements are retained.
# "Unshared Log Ratio" columns are excluded.
# Summary rows (Mean, Median, StdDev) are removed.

clean_proteome <- function(expr) {

  rownames(expr) <- expr$Gene

  expr <- expr[
    ,
    grepl(
      " Log Ratio$",
      colnames(expr)
    ) &
      !grepl(
        "Unshared",
        colnames(expr)
      ),
    drop = FALSE
  ]

  colnames(expr) <- sub(
    " Log Ratio$",
    "",
    colnames(expr)
  )

  expr <- expr[
    !rownames(expr) %in%
      c(
        "Mean",
        "Median",
        "StdDev"
      ),
    ,
    drop = FALSE
  ]

  expr
}


# Match proteome participant IDs to MAF participant IDs.
#
# The original analysis used:
#   1. direct identifier matching
#   2. for unmatched samples, matching after removing one leading
#      alphabetic character from the proteome identifier

match_maf_samples <- function(
    expr_matrix,
    maf_ids
) {

  smpl_res <- unique(
    maf_ids
  )

  smpl_expr <- colnames(
    expr_matrix
  )

  smpl_expr_res_1 <- intersect(
    smpl_expr,
    smpl_res
  )

  smpl_expr_rmn_1 <- setdiff(
    smpl_expr,
    smpl_res
  )

  smpl_res_rmn_1 <- setdiff(
    smpl_res,
    smpl_expr_res_1
  )

  smpl_expr_res_2 <- smpl_expr_rmn_1[
    gsub(
      "^([A-Za-z])",
      "",
      smpl_expr_rmn_1
    ) %in% smpl_res_rmn_1
  ]

  smpl_expr_res <- sort(
    c(
      smpl_expr_res_1,
      smpl_expr_res_2
    )
  )

  expr_matrix[
    ,
    smpl_expr_res,
    drop = FALSE
  ]
}


# Construct one normal or tumor expression matrix after specimen
# annotation and participant-ID conversion.

prepare_expression_subset <- function(
    expr,
    info_subset,
    specimen_column,
    participant_column,
    maf_ids
) {

  specimen_labels <- intersect(
    info_subset[[specimen_column]],
    colnames(expr)
  )

  expr_subset <- expr[
    ,
    specimen_labels,
    drop = FALSE
  ]

  id_map <- info_subset[[participant_column]]

  names(id_map) <- info_subset[[specimen_column]]

  colnames(expr_subset) <- id_map[
    colnames(expr_subset)
  ]

  match_maf_samples(
    expr_subset,
    maf_ids
  )
}

# -------------------------------------------------------------------------
# 4. Load pan-cancer somatic mutation data
# -------------------------------------------------------------------------

message("Loading CPTAC pan-cancer mutation data...")

maf_all <- fread(
  maf_file,
  header = TRUE,
  data.table = TRUE
)


# =========================================================================
# 5. BRCA
# =========================================================================

message("Preparing BRCA...")

maf <- maf_all[
  COHORT == "BRCA"
]

info <- read_excel(
  file.path(
    normal_tumor_dir,
    "BRCA",
    "S039_Breast_Cancer_Prospective_Collection_Specimens_r1.xlsx"
  )
)

expr <- fread(
  file.path(
    normal_tumor_dir,
    "BRCA",
    "CPTAC2_Breast_Prospective_Collection_BI_Proteome.tmt10.tsv"
  ),
  header = TRUE,
  data.table = FALSE
)

expr <- clean_proteome(
  expr
)


# Tumor

info_tumor <- info[
  info$`Sample Type` == "Tumor",
]

BRCA_tumor <- prepare_expression_subset(
  expr = expr,
  info_subset = info_tumor,
  specimen_column = "Specimen Label",
  participant_column =
    "Participant Protocol Identifier : Collection Protocol Registration",
  maf_ids = sub(
    "_T",
    "",
    unique(
      maf$Tumor_Sample_Barcode
    )
  )
)


# Normal

info_normal <- info[
  info$`Sample Type` == "Adjacent_Normal",
]

BRCA_normal <- prepare_expression_subset(
  expr = expr,
  info_subset = info_normal,
  specimen_column = "Specimen Label",
  participant_column =
    "Participant Protocol Identifier : Collection Protocol Registration",
  maf_ids = sub(
    "_N",
    "",
    unique(
      maf$Matched_Norm_Sample_Barcode
    )
  )
)

message(
  "BRCA: ",
  ncol(BRCA_normal),
  " normal, ",
  ncol(BRCA_tumor),
  " tumor"
)


# =========================================================================
# 6. ccRCC
# =========================================================================

message("Preparing ccRCC...")

# NOTE:
# The MAF cohort label is "CCRCC", whereas the manuscript/display
# cancer-type label is "ccRCC".

maf <- maf_all[
  COHORT == "CCRCC"
]

info <- read_excel(
  file.path(
    normal_tumor_dir,
    "ccRCC",
    "S044_CPTAC_CCRCC_Discovery_Cohort_Specimens_r1_Sept2018.xlsx"
  )
)

expr <- fread(
  file.path(
    normal_tumor_dir,
    "ccRCC",
    "CPTAC3_Clear_Cell_Renal_Cell_Carcinoma_Proteome.tmt10.tsv"
  ),
  header = TRUE,
  data.table = FALSE
)

expr <- clean_proteome(
  expr
)


# Tumor

info_tumor <- info[
  info$Group == "Tumor",
]

ccRCC_tumor <- prepare_expression_subset(
  expr = expr,
  info_subset = info_tumor,
  specimen_column = "Aliquot ID",
  participant_column = "ParticipantID",
  maf_ids = sub(
    "_T",
    "",
    unique(
      maf$Tumor_Sample_Barcode
    )
  )
)


# Normal

info_normal <- info[
  info$Group == "Normal",
]

ccRCC_normal <- prepare_expression_subset(
  expr = expr,
  info_subset = info_normal,
  specimen_column = "Aliquot ID",
  participant_column = "ParticipantID",
  maf_ids = sub(
    "_N",
    "",
    unique(
      maf$Matched_Norm_Sample_Barcode
    )
  )
)

message(
  "ccRCC: ",
  ncol(ccRCC_normal),
  " normal, ",
  ncol(ccRCC_tumor),
  " tumor"
)


# =========================================================================
# 7. COAD
# =========================================================================

message("Preparing COAD...")

maf <- maf_all[
  COHORT == "COAD"
]

info <- read_excel(
  file.path(
    normal_tumor_dir,
    "COAD",
    "CPTAC_Colon_Cancer_Prospective_Collection_Biospecimens_S037.xlsx"
  )
)

expr <- fread(
  file.path(
    normal_tumor_dir,
    "COAD",
    "CPTAC2_Colon_Prospective_Collection_PNNL_Proteome.tmt10.tsv"
  ),
  header = TRUE,
  data.table = FALSE
)

expr <- clean_proteome(
  expr
)

info$Group <- str_extract(
  info$`Specimen Title`,
  "^[^ ]+"
)


# Tumor

info_tumor <- info[
  info$Group == "Tumor",
]

COAD_tumor <- prepare_expression_subset(
  expr = expr,
  info_subset = info_tumor,
  specimen_column = "Specimen Label",
  participant_column = "Participant ID",
  maf_ids = sub(
    "_T",
    "",
    unique(
      maf$Tumor_Sample_Barcode
    )
  )
)


# Normal

info_normal <- info[
  info$Group == "Normal",
]

COAD_normal <- prepare_expression_subset(
  expr = expr,
  info_subset = info_normal,
  specimen_column = "Specimen Label",
  participant_column = "Participant ID",
  maf_ids = sub(
    "_N",
    "",
    unique(
      maf$Matched_Norm_Sample_Barcode
    )
  )
)

message(
  "COAD: ",
  ncol(COAD_normal),
  " normal, ",
  ncol(COAD_tumor),
  " tumor"
)


# =========================================================================
# 8. LSCC
# =========================================================================

message("Preparing LSCC...")

maf <- maf_all[
  COHORT == "LSCC"
]

info <- read_excel(
  file.path(
    normal_tumor_dir,
    "LSCC",
    "S063_S058updated_BI_CPTAC3_LSCC_Tumor_Normal_Mapping_r2.xlsx"
  )
)

expr <- fread(
  file.path(
    normal_tumor_dir,
    "LSCC",
    "CPTAC3_Lung_Squamous_Cell_Carcinoma_Proteome.tmt11.tsv"
  ),
  header = TRUE,
  data.table = FALSE
)

expr <- clean_proteome(
  expr
)


# Tumor

info_tumor <- info[
  info$`Tumor or Normal` == "Tumor",
]

LSCC_tumor <- prepare_expression_subset(
  expr = expr,
  info_subset = info_tumor,
  specimen_column = "Specimen_Label",
  participant_column = "Participant_ID",
  maf_ids = sub(
    "_T",
    "",
    unique(
      maf$Tumor_Sample_Barcode
    )
  )
)


# Normal

info_normal <- info[
  info$`Tumor or Normal` == "Normal",
]

LSCC_normal <- prepare_expression_subset(
  expr = expr,
  info_subset = info_normal,
  specimen_column = "Specimen_Label",
  participant_column = "Participant_ID",
  maf_ids = sub(
    "_N",
    "",
    unique(
      maf$Matched_Norm_Sample_Barcode
    )
  )
)

message(
  "LSCC: ",
  ncol(LSCC_normal),
  " normal, ",
  ncol(LSCC_tumor),
  " tumor"
)


# =========================================================================
# 9. LUAD
# =========================================================================

message("Preparing LUAD...")

maf <- maf_all[
  COHORT == "LUAD"
]

info <- read_excel(
  file.path(
    normal_tumor_dir,
    "LUAD",
    "S046_S056_BI_CPTAC3_LUAD_Discovery_Cohort_Samples_r2_July2020.xlsx"
  )
)

expr <- fread(
  file.path(
    normal_tumor_dir,
    "LUAD",
    "CPTAC3_Lung_Adeno_Carcinoma_Proteome.tmt10.tsv"
  ),
  header = TRUE,
  data.table = FALSE
)

expr <- clean_proteome(
  expr
)


# Tumor

info_tumor <- info[
  info$Type == "Tumor",
]

LUAD_tumor <- prepare_expression_subset(
  expr = expr,
  info_subset = info_tumor,
  specimen_column = "Aliquot (Specimen Label)",
  participant_column = "Participant ID (case_id)",
  maf_ids = sub(
    "_T",
    "",
    unique(
      maf$Tumor_Sample_Barcode
    )
  )
)


# Normal

info_normal <- info[
  info$Type == "Normal",
]

LUAD_normal <- prepare_expression_subset(
  expr = expr,
  info_subset = info_normal,
  specimen_column = "Aliquot (Specimen Label)",
  participant_column = "Participant ID (case_id)",
  maf_ids = sub(
    "_N",
    "",
    unique(
      maf$Matched_Norm_Sample_Barcode
    )
  )
)

message(
  "LUAD: ",
  ncol(LUAD_normal),
  " normal, ",
  ncol(LUAD_tumor),
  " tumor"
)


# =========================================================================
# 10. PDAC
# =========================================================================

message("Preparing PDAC...")

maf <- maf_all[
  COHORT == "PDAC"
]

info <- read_excel(
  file.path(
    normal_tumor_dir,
    "PDAC",
    "S061_CPTAC_PDA_Discovery_Cohort_Specimens_r1_Feb2021.xlsx"
  )
)

# Preserve the column selection used in the original analysis.
info <- info[
  ,
  which(
    names(info) == "TMT plex"
  ):
    which(
      names(info) == "Withdrawn"
    )
]

expr <- fread(
  file.path(
    normal_tumor_dir,
    "PDAC",
    "CPTAC3_Pancreatic_Ductal_Adenocarcinoma_Proteome.tmt11.tsv"
  ),
  header = TRUE,
  data.table = FALSE
)

expr <- clean_proteome(
  expr
)


# Tumor

info_tumor <- info[
  grepl(
    "^Tumor",
    info$`Tumor/Normal`
  ),
]

PDAC_tumor <- prepare_expression_subset(
  expr = expr,
  info_subset = info_tumor,
  specimen_column = "Aliquot ID",
  participant_column = "Case ID",
  maf_ids = sub(
    "_T",
    "",
    unique(
      maf$Tumor_Sample_Barcode
    )
  )
)


# Normal

info_normal <- info[
  grepl(
    "^Normal",
    info$`Tumor/Normal`
  ),
]

PDAC_normal <- prepare_expression_subset(
  expr = expr,
  info_subset = info_normal,
  specimen_column = "Aliquot ID",
  participant_column = "Case ID",
  maf_ids = sub(
    "_N",
    "",
    unique(
      maf$Matched_Norm_Sample_Barcode
    )
  )
)

message(
  "PDAC: ",
  ncol(PDAC_normal),
  " normal, ",
  ncol(PDAC_tumor),
  " tumor"
)


# =========================================================================
# 11. UCEC
# =========================================================================

message("Preparing UCEC...")

maf <- maf_all[
  COHORT == "UCEC"
]

info <- read_excel(
  file.path(
    normal_tumor_dir,
    "UCEC",
    "S043_CPTAC_UCEC_Discovery_Cohort_Study_Specimens_r1_Sept2018.xlsx"
  ),
  skip = 4
)

# Preserve the column selection used in the original analysis.
info <- info[
  ,
  which(
    names(info) == "TMT plex"
  ):
    which(
      names(info) == "TCIA Slide ID"
    )
]

expr <- fread(
  file.path(
    normal_tumor_dir,
    "UCEC",
    "CPTAC3_Uterine_Corpus_Endometrial_Carcinoma_Proteome.tmt10.tsv"
  ),
  header = TRUE,
  data.table = FALSE
)

expr <- clean_proteome(
  expr
)


# Tumor

info_tumor <- info[
  info$Group == "Tumor",
]

UCEC_tumor <- prepare_expression_subset(
  expr = expr,
  info_subset = info_tumor,
  specimen_column = "Aliquot ID",
  participant_column = "ParticipantID (Case_ID)",
  maf_ids = sub(
    "_T",
    "",
    unique(
      maf$Tumor_Sample_Barcode
    )
  )
)


# Normal

info_normal <- info[
  info$Group == "Normal",
]

UCEC_normal <- prepare_expression_subset(
  expr = expr,
  info_subset = info_normal,
  specimen_column = "Aliquot ID",
  participant_column = "ParticipantID (Case_ID)",
  maf_ids = sub(
    "_N",
    "",
    unique(
      maf$Matched_Norm_Sample_Barcode
    )
  )
)

message(
  "UCEC: ",
  ncol(UCEC_normal),
  " normal, ",
  ncol(UCEC_tumor),
  " tumor"
)


# -------------------------------------------------------------------------
# 12. Combine prepared matrices
# -------------------------------------------------------------------------

normal_tumor_data <- list(

  BRCA = list(
    normal = BRCA_normal,
    tumor = BRCA_tumor
  ),

  ccRCC = list(
    normal = ccRCC_normal,
    tumor = ccRCC_tumor
  ),

  COAD = list(
    normal = COAD_normal,
    tumor = COAD_tumor
  ),

  LSCC = list(
    normal = LSCC_normal,
    tumor = LSCC_tumor
  ),

  LUAD = list(
    normal = LUAD_normal,
    tumor = LUAD_tumor
  ),

  PDAC = list(
    normal = PDAC_normal,
    tumor = PDAC_tumor
  ),

  UCEC = list(
    normal = UCEC_normal,
    tumor = UCEC_tumor
  )
)


# -------------------------------------------------------------------------
# 13. Check sample numbers
# -------------------------------------------------------------------------

sample_counts <- data.frame(

  cancer = c(
    "BRCA",
    "ccRCC",
    "COAD",
    "LSCC",
    "LUAD",
    "PDAC",
    "UCEC"
  ),

  Normal = c(
    ncol(BRCA_normal),
    ncol(ccRCC_normal),
    ncol(COAD_normal),
    ncol(LSCC_normal),
    ncol(LUAD_normal),
    ncol(PDAC_normal),
    ncol(UCEC_normal)
  ),

  Tumor = c(
    ncol(BRCA_tumor),
    ncol(ccRCC_tumor),
    ncol(COAD_tumor),
    ncol(LSCC_tumor),
    ncol(LUAD_tumor),
    ncol(PDAC_tumor),
    ncol(UCEC_tumor)
  )
)

print(
  sample_counts
)


# Expected sample numbers from the original analysis.

expected_counts <- data.frame(

  cancer = c(
    "BRCA",
    "ccRCC",
    "COAD",
    "LSCC",
    "LUAD",
    "PDAC",
    "UCEC"
  ),

  Normal = c(
    18,
    84,
    95,
    100,
    99,
    66,
    30
  ),

  Tumor = c(
    113,
    110,
    94,
    108,
    106,
    137,
    95
  )
)

if (
  !all(
    sample_counts$Normal ==
      expected_counts$Normal
  ) ||
    !all(
      sample_counts$Tumor ==
        expected_counts$Tumor
    )
) {

  warning(
    "Sample counts differ from those in the original analysis."
  )

} else {

  message(
    "All normal and tumor sample counts match the original analysis."
  )
}


# -------------------------------------------------------------------------
# 14. Save prepared data
# -------------------------------------------------------------------------

saveRDS(
  normal_tumor_data,
  output_file
)

fwrite(
  sample_counts,
  sample_count_file,
  sep = "\t"
)

message(
  "Prepared proteome matrices written to: ",
  output_file
)

message(
  "Sample counts written to: ",
  sample_count_file
)
