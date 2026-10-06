# 09_mutation_rewiring_by_association_origin.R
#
# Compare mutation-associated weakening between Pre-existing and
# Tumor-acquired co-abundance associations.
#
# Mutation-associated rewiring information comes from script 08.
#
# Association origin comes ONLY from the adjacent-normal -> tumor
# association history calculated by script 03.
#
# IMPORTANT:
#
#   mutation-analysis baseline association
#
# and
#
#   tumor association from the normal/tumor analysis
#
# are kept as separate variables and are NOT forced to agree.
#
# Query -> target direction is preserved throughout.
#
# Inputs:
#
#   results/mutation_rewiring/individual/
#     individual_cancer_mutation_rewiring_edges.tsv
#
#   results/
#     exprs.cor.all_genes.partners.wt_only.<cancer>.normal.txt
#     exprs.cor.all_genes.partners.wt_only.<cancer>.tumor.txt
#
# Outputs:
#
#   results/mutation_origin/
#     mutation_edges_with_association_origin.tsv
#     association_origin_analysis_edges.tsv
#     association_origin_summary.tsv
#     association_origin_summary_by_cancer.tsv
#     baseline_tumor_concordance.tsv
#     weakening_model_coefficients.tsv
#     weakening_probabilities.tsv
#     weakening_origin_contrasts.tsv
#     weakening_baseline_contrasts.tsv
#     weakening_mutation_class_contrasts.tsv


library(data.table)
library(dplyr)
library(purrr)
library(emmeans)


# -------------------------------------------------------------------------
# 1. Parameters
# -------------------------------------------------------------------------

r_cut <- 0.5
fdr_cut <- 0.05
interaction_p_cut <- 0.05

cancer_types <- c(
  "brca",
  "ccrcc",
  "coad",
  "lscc",
  "luad",
  "pdac",
  "ucec"
)

cancer_labels <- c(
  brca = "BRCA",
  ccrcc = "ccRCC",
  coad = "COAD",
  lscc = "LSCC",
  luad = "LUAD",
  pdac = "PDAC",
  ucec = "UCEC"
)


# -------------------------------------------------------------------------
# 2. Paths
# -------------------------------------------------------------------------

mutation_file <- file.path(
  "results",
  "mutation_rewiring",
  "individual",
  "individual_cancer_mutation_rewiring_edges.tsv"
)

# Script 03 normal/tumor correlation files.
normal_tumor_dir <- "results"

output_dir <- file.path(
  "results",
  "mutation_origin"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# -------------------------------------------------------------------------
# 3. Load mutation-associated rewiring edges from script 08
# -------------------------------------------------------------------------

mutation_edges <- fread(
  mutation_file,
  header = TRUE,
  data.table = FALSE
)

required_mutation_columns <- c(
  "cancer",
  "query",
  "query_refseq",
  "target",
  "target_refseq",
  "mutation_class",
  "baseline_association",
  "wt_pairwise_n",
  "mut_pairwise_n",
  "wt_cor_r",
  "wt_cor_p_adj",
  "mut_cor_r",
  "mut_cor_p_adj",
  "delta_cor_r",
  "wt_slope",
  "mut_slope",
  "interaction_estimate",
  "interaction_se",
  "interaction_p",
  "interaction_p_adj",
  "rewiring_status"
)

missing_mutation_columns <- setdiff(
  required_mutation_columns,
  colnames(mutation_edges)
)

if (length(missing_mutation_columns) > 0) {
  stop(
    "Missing columns in script 08 output: ",
    paste(
      missing_mutation_columns,
      collapse = ", "
    )
  )
}


# -------------------------------------------------------------------------
# 4. Restrict to cancer types with adjacent normal tissues
# -------------------------------------------------------------------------

mutation_edges <- mutation_edges %>%
  filter(
    cancer %in% unname(cancer_labels)
  )


# -------------------------------------------------------------------------
# 5. Read normal/tumor association history for mutation-analysis edges
#
# IMPORTANT:
# Association origin is determined ONLY from:
#
#   Normal association -> Tumor association
#
# It is NOT determined from the mutation-analysis baseline association.
# -------------------------------------------------------------------------

read_normal_tumor_history <- function(cancer) {

  cancer_label <- unname(
    cancer_labels[cancer]
  )

  message(
    "Processing normal/tumor history: ",
    cancer_label
  )


  # -----------------------------------------------------------------------
  # Only query -> target pairs represented in script 08
  # -----------------------------------------------------------------------

  mutation_keys <- mutation_edges %>%
    filter(
      cancer == cancer_label
    ) %>%
    distinct(
      query,
      target
    )


  # -----------------------------------------------------------------------
  # File paths
  # -----------------------------------------------------------------------

  normal_file <- file.path(
    normal_tumor_dir,
    paste0(
      "exprs.cor.all_genes.partners.wt_only.",
      cancer,
      ".normal.txt"
    )
  )

  tumor_file <- file.path(
    normal_tumor_dir,
    paste0(
      "exprs.cor.all_genes.partners.wt_only.",
      cancer,
      ".tumor.txt"
    )
  )


  if (!file.exists(normal_file)) {
    stop(
      "Normal correlation file not found: ",
      normal_file
    )
  }

  if (!file.exists(tumor_file)) {
    stop(
      "Tumor correlation file not found: ",
      tumor_file
    )
  }


  # -----------------------------------------------------------------------
  # Normal
  # -----------------------------------------------------------------------

  res_normal <- fread(
    normal_file,
    select = c(
      "query",
      "target",
      "wt_cor_r",
      "wt_cor_p_adj"
    ),
    data.table = FALSE
  ) %>%

    semi_join(
      mutation_keys,
      by = c(
        "query",
        "target"
      )
    )


  # -----------------------------------------------------------------------
  # Tumor
  # -----------------------------------------------------------------------

  res_tumor <- fread(
    tumor_file,
    select = c(
      "query",
      "target",
      "wt_cor_r",
      "wt_cor_p_adj"
    ),
    data.table = FALSE
  ) %>%

    semi_join(
      mutation_keys,
      by = c(
        "query",
        "target"
      )
    )


  # -----------------------------------------------------------------------
  # Exact query -> target keys must be unique
  # -----------------------------------------------------------------------

  if (
    any(
      duplicated(
        res_normal[
          c(
            "query",
            "target"
          )
        ]
      )
    )
  ) {
    stop(
      "Duplicate query-target rows in normal file for ",
      cancer_label
    )
  }

  if (
    any(
      duplicated(
        res_tumor[
          c(
            "query",
            "target"
          )
        ]
      )
    )
  ) {
    stop(
      "Duplicate query-target rows in tumor file for ",
      cancer_label
    )
  }


  # -----------------------------------------------------------------------
  # Match the SAME directional query -> target pair
  # -----------------------------------------------------------------------

  inner_join(
    res_normal,
    res_tumor,
    by = c(
      "query",
      "target"
    ),
    suffix = c(
      "_normal",
      "_tumor"
    )
  ) %>%

    mutate(
      cancer =
        cancer_label,


      # -------------------------------------------------------------------
      # Normal association
      # -------------------------------------------------------------------

      normal_association = case_when(

        wt_cor_r_normal >= r_cut &
          wt_cor_p_adj_normal < fdr_cut ~
          "Positive",

        wt_cor_r_normal <= -r_cut &
          wt_cor_p_adj_normal < fdr_cut ~
          "Negative",

        TRUE ~
          "n.s."
      ),


      # -------------------------------------------------------------------
      # Tumor association
      # -------------------------------------------------------------------

      tumor_association = case_when(

        wt_cor_r_tumor >= r_cut &
          wt_cor_p_adj_tumor < fdr_cut ~
          "Positive",

        wt_cor_r_tumor <= -r_cut &
          wt_cor_p_adj_tumor < fdr_cut ~
          "Negative",

        TRUE ~
          "n.s."
      ),


      tumorigenesis_delta_cor_r =
        wt_cor_r_tumor -
        wt_cor_r_normal
    ) %>%

    select(
      cancer,
      query,
      target,

      wt_cor_r_normal,
      wt_cor_p_adj_normal,
      normal_association,

      wt_cor_r_tumor,
      wt_cor_p_adj_tumor,
      tumor_association,

      tumorigenesis_delta_cor_r
    )
}


normal_tumor_history <- map_dfr(
  cancer_types,
  read_normal_tumor_history
)


# -------------------------------------------------------------------------
# 6. Join normal/tumor history to mutation-associated rewiring
#
# Exact directional join:
#
#   cancer + query + target
#
# A -> B and B -> A remain distinct.
# -------------------------------------------------------------------------

edge_annotated <- mutation_edges %>%

  left_join(
    normal_tumor_history,
    by = c(
      "cancer",
      "query",
      "target"
    )
  )


# -------------------------------------------------------------------------
# 7. Define association origin
#
# Pre-existing:
#   Positive -> Positive
#   Negative -> Negative
#
# Tumor-acquired:
#   n.s. -> Positive
#   n.s. -> Negative
#
# Other histories are retained in the annotated output but are not
# included in the Pre-existing vs Tumor-acquired model.
# -------------------------------------------------------------------------

edge_annotated <- edge_annotated %>%

  mutate(
    association_origin = case_when(

      normal_association == "Positive" &
        tumor_association == "Positive" ~
        "Pre-existing",

      normal_association == "Negative" &
        tumor_association == "Negative" ~
        "Pre-existing",

      normal_association == "n.s." &
        tumor_association == "Positive" ~
        "Tumor-acquired",

      normal_association == "n.s." &
        tumor_association == "Negative" ~
        "Tumor-acquired",

      normal_association == "Positive" &
        tumor_association == "Negative" ~
        "Sign-switched",

      normal_association == "Negative" &
        tumor_association == "Positive" ~
        "Sign-switched",

      normal_association %in%
        c(
          "Positive",
          "Negative"
        ) &
        tumor_association == "n.s." ~
        "Lost",

      normal_association == "n.s." &
        tumor_association == "n.s." ~
        "n.s.-to-n.s.",

      TRUE ~
        NA_character_
    ),


    # Mutation-associated significance remains the script 08 test.
    significant_rewiring =
      !is.na(
        interaction_p
      ) &
      interaction_p <
        interaction_p_cut,


    # Diagnostic only.
    # Do NOT use this to redefine association origin.
    baseline_tumor_match =
      !is.na(
        tumor_association
      ) &
      baseline_association ==
        tumor_association
  )


# -------------------------------------------------------------------------
# 8. Diagnostic:
# mutation-analysis baseline vs tumor association
#
# The two variables come from different analyses and are intentionally
# retained separately.
# -------------------------------------------------------------------------

baseline_tumor_concordance <- edge_annotated %>%

  filter(
    association_origin %in%
      c(
        "Pre-existing",
        "Tumor-acquired"
      )
  ) %>%

  count(
    tumor_association,
    baseline_association,
    name = "n"
  ) %>%

  group_by(
    tumor_association
  ) %>%

  mutate(
    fraction =
      n / sum(n)
  ) %>%

  ungroup()


# -------------------------------------------------------------------------
# 9. Analysis dataset
#
# Mutation-associated Baseline association is preserved for defining
# Strengthening/Weakening.
#
# Association origin is independently defined from Normal -> Tumor.
# -------------------------------------------------------------------------

edge_test <- edge_annotated %>%

  filter(
    association_origin %in%
      c(
        "Pre-existing",
        "Tumor-acquired"
      ),

    baseline_association %in%
      c(
        "Positive",
        "Negative"
      ),

    mutation_class %in%
      c(
        "Benign",
        "Pathogenic"
      )
  ) %>%

  mutate(
    association_origin = factor(
      association_origin,
      levels = c(
        "Pre-existing",
        "Tumor-acquired"
      )
    ),

    baseline_association = factor(
      baseline_association,
      levels = c(
        "Negative",
        "Positive"
      )
    ),

    mutation_class = factor(
      mutation_class,
      levels = c(
        "Benign",
        "Pathogenic"
      )
    ),

    cancer = factor(
      cancer,
      levels = unname(
        cancer_labels
      )
    ),

    weakened =
      rewiring_status ==
        "Weakening"
  )


# -------------------------------------------------------------------------
# 10. Descriptive summary
#
# Strengthening is retained here even though the inferential model below
# specifically tests mutation-associated weakening.
# -------------------------------------------------------------------------

association_origin_summary <- edge_test %>%

  group_by(
    mutation_class,
    baseline_association,
    association_origin
  ) %>%

  summarise(
    n =
      n(),

    n_rewired =
      sum(
        significant_rewiring,
        na.rm = TRUE
      ),

    fraction_rewired =
      mean(
        significant_rewiring,
        na.rm = TRUE
      ),

    n_weakened =
      sum(
        rewiring_status ==
          "Weakening",
        na.rm = TRUE
      ),

    fraction_weakened =
      n_weakened / n,

    n_strengthened =
      sum(
        rewiring_status ==
          "Strengthening",
        na.rm = TRUE
      ),

    fraction_strengthened =
      n_strengthened / n,

    .groups = "drop"
  )


# -------------------------------------------------------------------------
# 11. Descriptive summary by cancer type
# -------------------------------------------------------------------------

association_origin_summary_by_cancer <- edge_test %>%

  group_by(
    cancer,
    mutation_class,
    baseline_association,
    association_origin
  ) %>%

  summarise(
    n =
      n(),

    n_rewired =
      sum(
        significant_rewiring,
        na.rm = TRUE
      ),

    fraction_rewired =
      mean(
        significant_rewiring,
        na.rm = TRUE
      ),

    n_weakened =
      sum(
        rewiring_status ==
          "Weakening",
        na.rm = TRUE
      ),

    fraction_weakened =
      n_weakened / n,

    n_strengthened =
      sum(
        rewiring_status ==
          "Strengthening",
        na.rm = TRUE
      ),

    fraction_strengthened =
      n_strengthened / n,

    .groups = "drop"
  )


# -------------------------------------------------------------------------
# 12. Main binomial logistic model
#
# Outcome:
#   mutation-associated weakening
#
# Predictors:
#   Association origin
#   Baseline association
#   Mutation class
#   Cancer type
# -------------------------------------------------------------------------

fit_weakening <- glm(
  weakened ~
    association_origin *
    baseline_association *
    mutation_class +
    cancer,
  data = edge_test,
  family = binomial()
)


# -------------------------------------------------------------------------
# 13. Save model coefficients
# -------------------------------------------------------------------------

coef_table <- as.data.frame(
  summary(
    fit_weakening
  )$coefficients
)

weakening_model_coefficients <- tibble(
  term =
    rownames(
      coef_table
    ),

  estimate =
    coef_table[
      ,
      "Estimate"
    ],

  std_error =
    coef_table[
      ,
      "Std. Error"
    ],

  z_value =
    coef_table[
      ,
      "z value"
    ],

  p_value =
    coef_table[
      ,
      "Pr(>|z|)"
    ]
)


# -------------------------------------------------------------------------
# 14. Estimated marginal probabilities of weakening
#
# Averaged over cancer types, matching the previous emmeans analysis.
# -------------------------------------------------------------------------

emm_weakening_probability <- emmeans(
  fit_weakening,
  ~ association_origin *
    baseline_association *
    mutation_class,
  type = "response"
)

weakening_probabilities <- summary(
  emm_weakening_probability,
  infer = TRUE
) %>%
  as.data.frame()


# -------------------------------------------------------------------------
# 15. Tumor-acquired vs Pre-existing
#
# Main association-origin contrast.
#
# Four contrasts:
#   Negative + Benign
#   Positive + Benign
#   Negative + Pathogenic
#   Positive + Pathogenic
# -------------------------------------------------------------------------

emm_origin <- emmeans(
  fit_weakening,
  ~ association_origin |
    baseline_association *
    mutation_class
)

origin_contrast <- contrast(
  emm_origin,
  method = "revpairwise"
)

weakening_origin_contrasts <- summary(
  origin_contrast,
  infer = TRUE,
  type = "response"
) %>%
  as.data.frame()


# -------------------------------------------------------------------------
# 16. Positive vs Negative baseline association
#
# Useful for testing whether weakening differs by association sign
# within each origin x mutation-class combination.
# -------------------------------------------------------------------------

emm_baseline <- emmeans(
  fit_weakening,
  ~ baseline_association |
    association_origin *
    mutation_class
)

baseline_contrast <- contrast(
  emm_baseline,
  method = "revpairwise"
)

weakening_baseline_contrasts <- summary(
  baseline_contrast,
  infer = TRUE,
  type = "response"
) %>%
  as.data.frame()


# -------------------------------------------------------------------------
# 17. Pathogenic vs Benign
#
# Tested separately within each association-origin x baseline-sign group.
# -------------------------------------------------------------------------

emm_mutation <- emmeans(
  fit_weakening,
  ~ mutation_class |
    association_origin *
    baseline_association
)

mutation_contrast <- contrast(
  emm_mutation,
  method = "revpairwise"
)

weakening_mutation_class_contrasts <- summary(
  mutation_contrast,
  infer = TRUE,
  type = "response"
) %>%
  as.data.frame()


# -------------------------------------------------------------------------
# 18. Convert factor columns back to character before saving edge tables
# -------------------------------------------------------------------------

edge_test_output <- edge_test %>%

  mutate(
    cancer =
      as.character(
        cancer
      ),

    mutation_class =
      as.character(
        mutation_class
      ),

    baseline_association =
      as.character(
        baseline_association
      ),

    association_origin =
      as.character(
        association_origin
      )
  )


# -------------------------------------------------------------------------
# 19. Save
# -------------------------------------------------------------------------

fwrite(
  edge_annotated,
  file = file.path(
    output_dir,
    "mutation_edges_with_association_origin.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

fwrite(
  edge_test_output,
  file = file.path(
    output_dir,
    "association_origin_analysis_edges.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

fwrite(
  association_origin_summary,
  file = file.path(
    output_dir,
    "association_origin_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

fwrite(
  association_origin_summary_by_cancer,
  file = file.path(
    output_dir,
    "association_origin_summary_by_cancer.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

fwrite(
  baseline_tumor_concordance,
  file = file.path(
    output_dir,
    "baseline_tumor_concordance.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

fwrite(
  weakening_model_coefficients,
  file = file.path(
    output_dir,
    "weakening_model_coefficients.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

fwrite(
  weakening_probabilities,
  file = file.path(
    output_dir,
    "weakening_probabilities.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

fwrite(
  weakening_origin_contrasts,
  file = file.path(
    output_dir,
    "weakening_origin_contrasts.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

fwrite(
  weakening_baseline_contrasts,
  file = file.path(
    output_dir,
    "weakening_baseline_contrasts.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

fwrite(
  weakening_mutation_class_contrasts,
  file = file.path(
    output_dir,
    "weakening_mutation_class_contrasts.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# -------------------------------------------------------------------------
# 20. Reproducibility summary
# -------------------------------------------------------------------------

message(
  "Script 08 edges from seven normal/tumor cancer types: ",
  nrow(
    mutation_edges
  )
)

message(
  "Edges with matched normal/tumor history: ",
  sum(
    !is.na(
      edge_annotated$normal_association
    ) &
      !is.na(
        edge_annotated$tumor_association
      )
  )
)

message(
  "Edges used in Pre-existing/Tumor-acquired model: ",
  nrow(
    edge_test
  )
)

print(
  baseline_tumor_concordance,
  n = Inf,
  width = Inf
)

print(
  association_origin_summary,
  n = Inf,
  width = Inf
)

print(
  weakening_origin_contrasts,
  row.names = FALSE
)

message(
  "Association-origin analysis completed."
)
