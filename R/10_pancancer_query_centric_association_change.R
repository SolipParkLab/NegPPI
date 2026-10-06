# 10_pancancer_query_centric_association_change.R
#
# Pan-cancer query-centric analysis of mutation-associated
# co-abundance association changes.
#
# Input:
#   results/mutation_rewiring/
#     pancancer_mutation_rewiring_edges.tsv
#
# IMPORTANT:
#
# WT baseline associations must be significant:
#
#   Positive:
#     WT r > 0.5 and WT FDR < 0.05
#
#   Negative:
#     WT r < -0.5 and WT FDR < 0.05
#
# However, NO significance threshold is applied to the
# WT -> MUT association change.
#
# Association-change direction is determined ONLY from
# the sign of the interaction estimate:
#
#   WT Positive:
#     interaction estimate > 0 -> Strengthening
#     interaction estimate < 0 -> Weakening
#
#   WT Negative:
#     interaction estimate < 0 -> Strengthening
#     interaction estimate > 0 -> Weakening
#
# interaction P values are retained for reference but are
# NOT used to classify association change.
#
# Outputs:
#
#   results/query_centric/
#     pancancer_query_association_change_edges.tsv
#     pancancer_query_association_change.tsv
#     pancancer_query_association_change_glm.tsv
#     pancancer_query_association_change_ecdf_tests.tsv
#     pathogenic_query_sign_bias.tsv


library(data.table)
library(dplyr)


# -------------------------------------------------------------------------
# 1. Parameters
# -------------------------------------------------------------------------

permutation_n <- 10000
permutation_seed <- 1

min_associations_sign_bias <- 3


# -------------------------------------------------------------------------
# 2. Paths
# -------------------------------------------------------------------------

input_file <- file.path(
  "results",
  "mutation_rewiring",
  "pancancer_mutation_rewiring_edges.tsv"
)

output_dir <- file.path(
  "results",
  "query_centric"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# -------------------------------------------------------------------------
# 3. Load pan-cancer association-level results from script 07
# -------------------------------------------------------------------------

edges <- fread(
  input_file,
  header = TRUE,
  data.table = FALSE
)


# -------------------------------------------------------------------------
# 4. Check required columns
# -------------------------------------------------------------------------

required_columns <- c(
  "query",
  "query_refseq",
  "target",
  "target_refseq",
  "mutation_class",
  "baseline_association",
  "wt_cor_r",
  "wt_cor_p_adj",
  "mut_cor_r",
  "mut_cor_p_adj",
  "delta_cor_r",
  "interaction_estimate",
  "interaction_se",
  "interaction_p",
  "interaction_p_adj"
)

missing_columns <- setdiff(
  required_columns,
  colnames(edges)
)

if (length(missing_columns) > 0) {
  stop(
    "Missing required columns: ",
    paste(
      missing_columns,
      collapse = ", "
    )
  )
}


# -------------------------------------------------------------------------
# 5. Retain the four query-analysis groups
#
# Script 07 should already contain only significant WT Positive/Negative
# baseline associations, but this filter is kept explicit.
# -------------------------------------------------------------------------

edges <- edges %>%

  filter(
    mutation_class %in%
      c(
        "Benign",
        "Pathogenic"
      ),

    baseline_association %in%
      c(
        "Positive",
        "Negative"
      )
  )


# -------------------------------------------------------------------------
# 6. Classify association-change direction
#
# IMPORTANT:
#
# interaction_p and interaction_p_adj are NOT used here.
# -------------------------------------------------------------------------

edges <- edges %>%

  mutate(
    association_change = case_when(

      !is.finite(interaction_estimate) ~
        NA_character_,

      baseline_association == "Positive" &
        interaction_estimate > 0 ~
        "Strengthening",

      baseline_association == "Positive" &
        interaction_estimate < 0 ~
        "Weakening",

      baseline_association == "Negative" &
        interaction_estimate < 0 ~
        "Strengthening",

      baseline_association == "Negative" &
        interaction_estimate > 0 ~
        "Weakening",

      interaction_estimate == 0 ~
        "No directional change",

      TRUE ~
        NA_character_
    )
  )


# -------------------------------------------------------------------------
# 7. Save association-level directional classification
# -------------------------------------------------------------------------

edge_output <- edges %>%

  select(
    query,
    query_refseq,
    target,
    target_refseq,

    mutation_class,
    baseline_association,

    wt_cor_r,
    wt_cor_p_adj,

    mut_cor_r,
    mut_cor_p_adj,

    delta_cor_r,

    interaction_estimate,
    interaction_se,
    interaction_p,
    interaction_p_adj,

    association_change
  ) %>%

  arrange(
    mutation_class,
    baseline_association,
    query,
    target
  )


fwrite(
  edge_output,
  file = file.path(
    output_dir,
    "pancancer_query_association_change_edges.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# =========================================================================
# 8. Query-centric summary
#
# Unit:
#
#   query
#   x mutation class
#   x WT baseline association
#
# Fractions are calculated among associations with a directional
# interaction estimate.
# =========================================================================

query_summary <- edges %>%

  group_by(
    query,
    mutation_class,
    baseline_association
  ) %>%

  summarise(
    n_associations =
      n(),

    n_directional_associations =
      sum(
        association_change %in%
          c(
            "Strengthening",
            "Weakening"
          ),
        na.rm = TRUE
      ),

    n_strengthened =
      sum(
        association_change ==
          "Strengthening",
        na.rm = TRUE
      ),

    n_weakened =
      sum(
        association_change ==
          "Weakening",
        na.rm = TRUE
      ),

    n_no_directional_change =
      sum(
        association_change ==
          "No directional change",
        na.rm = TRUE
      ),

    n_not_evaluable =
      sum(
        is.na(
          association_change
        )
      ),

    .groups = "drop"
  ) %>%

  mutate(
    fraction_strengthened =
      ifelse(
        n_directional_associations > 0,
        n_strengthened /
          n_directional_associations,
        NA_real_
      ),

    fraction_weakened =
      ifelse(
        n_directional_associations > 0,
        n_weakened /
          n_directional_associations,
        NA_real_
      )
  )


# -------------------------------------------------------------------------
# 9. Order query-level output
# -------------------------------------------------------------------------

query_summary <- query_summary %>%

  mutate(
    mutation_class = factor(
      mutation_class,
      levels = c(
        "Benign",
        "Pathogenic"
      )
    ),

    baseline_association = factor(
      baseline_association,
      levels = c(
        "Positive",
        "Negative"
      )
    )
  ) %>%

  arrange(
    mutation_class,
    baseline_association,
    query
  ) %>%

  mutate(
    mutation_class =
      as.character(
        mutation_class
      ),

    baseline_association =
      as.character(
        baseline_association
      )
  )


fwrite(
  query_summary,
  file = file.path(
    output_dir,
    "pancancer_query_association_change.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# =========================================================================
# 10. Pathogenic vs Benign:
#     query-centric weakening odds
#
# Separate models are fitted for Positive and Negative WT baseline
# associations.
#
# Outcome:
#
#   Weakening vs Strengthening
#
# Predictor of interest:
#
#   Pathogenic vs Benign
#
# log1p(n_directional_associations) is included to preserve the
# previous query-centric adjustment for the number of evaluable
# partner associations.
# =========================================================================

fit_query_model <- function(baseline_value) {

  model_data <- query_summary %>%

    filter(
      baseline_association ==
        baseline_value,

      n_directional_associations > 0
    ) %>%

    mutate(
      mutation_class = factor(
        mutation_class,
        levels = c(
          "Benign",
          "Pathogenic"
        )
      )
    )


  if (
    n_distinct(
      model_data$mutation_class
    ) < 2
  ) {
    stop(
      "Both mutation classes are not represented for baseline: ",
      baseline_value
    )
  }


  fit <- glm(
    cbind(
      n_weakened,
      n_strengthened
    ) ~
      mutation_class +
      log1p(
        n_directional_associations
      ),
    data = model_data,
    family = binomial()
  )


  coef_table <- summary(
    fit
  )$coefficients


  target_term <- "mutation_classPathogenic"


  if (
    !target_term %in%
      rownames(
        coef_table
      )
  ) {
    stop(
      "Pathogenic coefficient not found for baseline: ",
      baseline_value
    )
  }


  estimate <- coef_table[target_term, "Estimate"]

  std_error <- coef_table[target_term, "Std. Error"]

  z_value <- coef_table[target_term, "z value"]

  p_value <- coef_table[target_term, "Pr(>|z|)"]


  data.frame(
    baseline_association =
      baseline_value,

    contrast =
      "Pathogenic / Benign",

    log_odds_ratio =
      estimate,

    std_error =
      std_error,

    odds_ratio =
      exp(
        estimate
      ),

    lower_95_ci =
      exp(
        estimate -
          1.96 *
          std_error
      ),

    upper_95_ci =
      exp(
        estimate +
          1.96 *
          std_error
      ),

    z_value =
      z_value,

    p_value =
      p_value,

    n_query_rows =
      nrow(
        model_data
      ),

    n_benign_queries =
      sum(
        model_data$mutation_class ==
          "Benign"
      ),

    n_pathogenic_queries =
      sum(
        model_data$mutation_class ==
          "Pathogenic"
      )
  )
}


glm_results <- bind_rows(
  fit_query_model(
    "Positive"
  ),
  fit_query_model(
    "Negative"
  )
)


fwrite(
  glm_results,
  file = file.path(
    output_dir,
    "pancancer_query_association_change_glm.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# =========================================================================
# 11. Query-level ECDF permutation tests
#
# Compare the distribution of fraction_weakened between
# Benign and Pathogenic queries.
#
# The permutation test uses the KS statistic as in the previous
# query-centric analysis.
# =========================================================================

perm_ecdf <- function(
    x,
    y,
    B = 10000,
    seed = 1
) {

  x <- x[
    is.finite(
      x
    )
  ]

  y <- y[
    is.finite(
      y
    )
  ]


  if (
    length(x) < 2 ||
    length(y) < 2
  ) {
    return(
      data.frame(
        statistic = NA_real_,
        p_value = NA_real_,
        n_x = length(x),
        n_y = length(y),
        permutations = B
      )
    )
  }


  d_obs <- suppressWarnings(
    ks.test(
      x,
      y
    )$statistic
  )

  d_obs <- unname(
    d_obs
  )


  pooled <- c(
    x,
    y
  )

  nx <- length(
    x
  )

  ny <- length(
    y
  )


  set.seed(
    seed
  )


  d_perm <- numeric(
    B
  )


  for (i in seq_len(B)) {

    idx <- sample.int(
      nx + ny
    )

    x_perm <- pooled[
      idx[
        seq_len(
          nx
        )
      ]
    ]

    y_perm <- pooled[
      idx[
        (nx + 1):
          (nx + ny)
      ]
    ]

    d_perm[i] <- suppressWarnings(
      unname(
        ks.test(
          x_perm,
          y_perm
        )$statistic
      )
    )
  }


  p_value <- mean(
    d_perm >=
      d_obs
  )


  data.frame(
    statistic =
      d_obs,

    p_value =
      p_value,

    n_x =
      nx,

    n_y =
      ny,

    permutations =
      B
  )
}


run_ecdf_test <- function(
    baseline_value
) {

  x <- query_summary %>%

    filter(
      baseline_association ==
        baseline_value,

      mutation_class ==
        "Benign"
    ) %>%

    pull(
      fraction_weakened
    )


  y <- query_summary %>%

    filter(
      baseline_association ==
        baseline_value,

      mutation_class ==
        "Pathogenic"
    ) %>%

    pull(
      fraction_weakened
    )


  result <- perm_ecdf(
    x = x,
    y = y,
    B = permutation_n,
    seed = permutation_seed
  )


  result %>%

    mutate(
      baseline_association =
        baseline_value,

      comparison =
        "Benign vs Pathogenic",

      n_benign =
        n_x,

      n_pathogenic =
        n_y
    ) %>%

    select(
      baseline_association,
      comparison,
      statistic,
      p_value,
      n_benign,
      n_pathogenic,
      permutations
    )
}


ecdf_results <- bind_rows(
  run_ecdf_test(
    "Positive"
  ),
  run_ecdf_test(
    "Negative"
  )
)


fwrite(
  ecdf_results,
  file = file.path(
    output_dir,
    "pancancer_query_association_change_ecdf_tests.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# =========================================================================
# 12. Within-query Positive-vs-Negative comparison
#
# Pathogenic mutations only.
#
# Restrict to queries with at least the specified number of directional
# associations in BOTH Positive and Negative WT baseline groups.
#
# Positive delta_weakening:
#
#   weakening is more frequent among Negative baseline associations.
#
# Positive log2_weakening_odds_ratio:
#
#   weakening odds are higher among Negative than Positive
#   baseline associations for the same query.
# =========================================================================

pathogenic_query_sign_bias <- query_summary %>%

  filter(
    mutation_class ==
      "Pathogenic"
  ) %>%

  select(
    query,
    baseline_association,
    n_associations,
    n_directional_associations,
    n_strengthened,
    n_weakened,
    fraction_strengthened,
    fraction_weakened
  ) %>%

  tidyr::pivot_wider(
    names_from =
      baseline_association,

    values_from = c(
      n_associations,
      n_directional_associations,
      n_strengthened,
      n_weakened,
      fraction_strengthened,
      fraction_weakened
    ),

    names_glue =
      "{.value}_{baseline_association}"
  ) %>%

  filter(
    n_directional_associations_Positive >=
      min_associations_sign_bias,

    n_directional_associations_Negative >=
      min_associations_sign_bias
  ) %>%

  mutate(
    delta_weakening =
      fraction_weakened_Negative -
      fraction_weakened_Positive,

    log2_weakening_odds_ratio =
      log2(
        (
          (
            n_weakened_Negative +
              0.5
          ) /
            (
              n_strengthened_Negative +
                0.5
            )
        ) /
          (
            (
              n_weakened_Positive +
                0.5
            ) /
              (
                n_strengthened_Positive +
                  0.5
              )
          )
      )
  ) %>%

  arrange(
    desc(
      delta_weakening
    )
  )


fwrite(
  pathogenic_query_sign_bias,
  file = file.path(
    output_dir,
    "pathogenic_query_sign_bias.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# -------------------------------------------------------------------------
# 13. Reproducibility summary
# -------------------------------------------------------------------------

message(
  "Pan-cancer baseline Positive/Negative associations: ",
  nrow(
    edges
  )
)


message(
  "Directional association changes: ",
  sum(
    edges$association_change %in%
      c(
        "Strengthening",
        "Weakening"
      ),
    na.rm = TRUE
  )
)


message(
  "Query-level rows: ",
  nrow(
    query_summary
  )
)


message(
  "Pathogenic queries with at least ",
  min_associations_sign_bias,
  " Positive and ",
  min_associations_sign_bias,
  " Negative directional associations: ",
  nrow(
    pathogenic_query_sign_bias
  )
)


print(
  glm_results,
  row.names = FALSE
)


print(
  ecdf_results,
  row.names = FALSE
)


message(
  "Pan-cancer query-centric association-change analysis completed."
)
