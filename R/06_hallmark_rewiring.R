# 06_hallmark_rewiring.R
#
# Annotate retained normal-derived co-abundance modules with the 21 C3PO
# cancer hallmarks and project module-pair rewiring onto hallmark pairs.
#
# Inputs:
#   results/modules/module_assignments.tsv
#   results/module_rewiring/module_pair_rewiring.tsv
#   data/raw/Li_2023/mmc7.xlsx
#
# Outputs:
#   results/hallmark_rewiring/module_hallmark_loadings.tsv
#   results/hallmark_rewiring/hallmark_pair_rewiring.tsv
#   results/hallmark_rewiring/hallmark_pair_rewiring_significant.tsv


library(data.table)
library(dplyr)
library(tidyr)
library(readxl)
library(purrr)


# -------------------------------------------------------------------------
# 1. Parameters
# -------------------------------------------------------------------------

min_module_size <- 50
n_perm <- 1000
eps <- 1e-12

cancer_types <- c(
  "BRCA",
  "ccRCC",
  "COAD",
  "LSCC",
  "LUAD",
  "PDAC",
  "UCEC"
)


# -------------------------------------------------------------------------
# 2. Paths
# -------------------------------------------------------------------------

module_file <- file.path(
  "results",
  "modules",
  "module_assignments.tsv"
)

module_rewiring_file <- file.path(
  "results",
  "module_rewiring",
  "module_pair_rewiring.tsv"
)

hallmark_file <- file.path(
  "data",
  "raw",
  "Li_2023",
  "mmc7.xlsx"
)

output_dir <- file.path(
  "results",
  "hallmark_rewiring"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

loading_output_file <- file.path(
  output_dir,
  "module_hallmark_loadings.tsv"
)

pair_output_file <- file.path(
  output_dir,
  "hallmark_pair_rewiring.tsv"
)

significant_output_file <- file.path(
  output_dir,
  "hallmark_pair_rewiring_significant.tsv"
)


# -------------------------------------------------------------------------
# 3. Load normal-derived module assignments
# -------------------------------------------------------------------------

module_assignments <- fread(
  module_file,
  header = TRUE
) %>%
  as_tibble() %>%
  mutate(
    cancer = as.character(cancer),
    module = as.integer(module),
    module_size = as.integer(module_size),
    protein = as.character(protein)
  )

retained_membership <- module_assignments %>%
  filter(
    module_size >= min_module_size
  )

module_meta <- retained_membership %>%
  distinct(
    cancer,
    module,
    module_size
  ) %>%
  mutate(
    cancer = factor(
      cancer,
      levels = cancer_types
    )
  ) %>%
  arrange(
    cancer,
    module
  ) %>%
  mutate(
    cancer = as.character(cancer)
  )


# -------------------------------------------------------------------------
# 4. Load C3PO 21 hallmark gene sets
# -------------------------------------------------------------------------

hallmark_tbl <- read_excel(
  hallmark_file,
  sheet = "Table S7B",
  skip = 1
)

colnames(hallmark_tbl)[1] <- "gene"

# Preserve the order of the hallmark columns in Table S7B.
hallmark_order <- colnames(hallmark_tbl)[-1]

hallmark_21_df <- hallmark_tbl %>%
  pivot_longer(
    cols = -gene,
    names_to = "term",
    values_to = "flag"
  ) %>%
  filter(
    flag == "Y"
  ) %>%
  transmute(
    term = as.character(term),
    gene = as.character(gene)
  ) %>%
  distinct()

terms_use <- sort(
  unique(hallmark_21_df$term)
)

if (length(terms_use) != 21) {
  stop(
    "Expected 21 C3PO hallmarks, but found ",
    length(terms_use),
    "."
  )
}

term_pairs <- expand_grid(
  term_1 = terms_use,
  term_2 = terms_use
) %>%
  filter(
    term_1 <= term_2
  )


# -------------------------------------------------------------------------
# 5. Calculate module x hallmark overlap
# -------------------------------------------------------------------------

overlap_tbl <- retained_membership %>%
  select(
    cancer,
    module,
    protein
  ) %>%
  inner_join(
    hallmark_21_df,
    by = c(
      "protein" = "gene"
    )
  ) %>%
  count(
    cancer,
    module,
    term,
    name = "overlap_n"
  )


# -------------------------------------------------------------------------
# 6. Calculate module-specific hallmark loadings
#
# Initial score:
#   overlap / module size
#
# Loading:
#   initial score / sum of initial scores across the 21 hallmarks
# -------------------------------------------------------------------------

module_hallmark_grid <- expand_grid(
  module_meta,
  term = hallmark_order
)

module_hallmark_loadings <- module_hallmark_grid %>%
  left_join(
    overlap_tbl,
    by = c(
      "cancer",
      "module",
      "term"
    )
  ) %>%
  mutate(
    overlap_n = replace_na(
      overlap_n,
      0L
    ),

    hallmark_overlap_fraction =
      overlap_n / module_size
  ) %>%
  group_by(
    cancer,
    module
  ) %>%
  mutate(
    hallmark_loading =
      hallmark_overlap_fraction /
      pmax(
        sum(hallmark_overlap_fraction),
        eps
      )
  ) %>%
  ungroup() %>%
  mutate(
    module_label = paste0(
      "M",
      module
    ),

    term = factor(
      term,
      levels = hallmark_order
    ),

    cancer = factor(
      cancer,
      levels = cancer_types
    )
  ) %>%
  arrange(
    cancer,
    module,
    term
  ) %>%
  transmute(
    cancer = as.character(cancer),
    module = module_label,
    module_size,
    cancer_hallmark = as.character(term),
    overlapping_proteins_n = overlap_n,
    hallmark_overlap_fraction,
    hallmark_loading
  )


# -------------------------------------------------------------------------
# 7. Check loading normalization
# -------------------------------------------------------------------------

loading_check <- module_hallmark_loadings %>%
  group_by(
    cancer,
    module
  ) %>%
  summarise(
    loading_sum = sum(
      hallmark_loading
    ),
    .groups = "drop"
  )

bad_loading <- loading_check %>%
  filter(
    abs(loading_sum - 1) > 1e-10 &
      abs(loading_sum) > 1e-10
  )

if (nrow(bad_loading) > 0) {
  stop(
    "Hallmark-loading normalization failed for ",
    nrow(bad_loading),
    " modules."
  )
}


# -------------------------------------------------------------------------
# 8. Load module-pair rewiring results from script 05
# -------------------------------------------------------------------------

module_pair_rewiring <- fread(
  module_rewiring_file,
  header = TRUE
) %>%
  as_tibble() %>%
  mutate(
    cancer = as.character(cancer),
    module_a = as.character(module_a),
    module_b = as.character(module_b)
  )


# -------------------------------------------------------------------------
# 9. Hallmark-pair projection
# -------------------------------------------------------------------------

compute_pair_rank_stats <- function(
    score_mat,
    edge_df,
    term_pairs
) {

  A <- score_mat[
    edge_df$module_a,
    ,
    drop = FALSE
  ]

  B <- score_mat[
    edge_df$module_b,
    ,
    drop = FALSE
  ]

  map_dfr(
    seq_len(
      nrow(term_pairs)
    ),
    function(k) {

      t1 <- term_pairs$term_1[k]
      t2 <- term_pairs$term_2[k]

      if (t1 == t2) {

        w <-
          A[, t1] *
          B[, t1]

      } else {

        w <-
          A[, t1] *
          B[, t2] +
          A[, t2] *
          B[, t1]
      }

      keep <- w > 0
      wsum <- sum(
        w,
        na.rm = TRUE
      )

      tibble(
        term_1 = t1,
        term_2 = t2,

        n_edges_with_weight =
          sum(keep),

        weight_sum =
          wsum,

        weighted_raw_frac_change =
          ifelse(
            wsum > 0,

            weighted.mean(
              edge_df$fractional_change,
              w,
              na.rm = TRUE
            ),

            NA_real_
          ),

        weighted_rank =
          ifelse(
            wsum > 0,

            weighted.mean(
              edge_df$within_cancer_rewiring_rank,
              w,
              na.rm = TRUE
            ),

            NA_real_
          ),

        weighted_rank_centered =
          ifelse(
            wsum > 0,

            weighted.mean(
              edge_df$centered_within_cancer_rewiring_rank,
              w,
              na.rm = TRUE
            ),

            NA_real_
          ),

        weighted_frac_relatively_weakened =
          ifelse(
            wsum > 0,

            weighted.mean(
              as.numeric(
                edge_df$
                  centered_within_cancer_rewiring_rank < 0
              ),
              w,
              na.rm = TRUE
            ),

            NA_real_
          )
      )
    }
  )
}


# -------------------------------------------------------------------------
# 10. Analyze one cancer type
# -------------------------------------------------------------------------

run_hallmark_rewiring <- function(
    cancer_name
) {

  message(
    "Processing ",
    cancer_name,
    "..."
  )

  edge_df <- module_pair_rewiring %>%
    filter(
      cancer == cancer_name
    )

  loading_df <- module_hallmark_loadings %>%
    filter(
      cancer == cancer_name
    ) %>%
    select(
      module,
      cancer_hallmark,
      hallmark_loading
    )

  score_df <- loading_df %>%
    pivot_wider(
      names_from = cancer_hallmark,
      values_from = hallmark_loading,
      values_fill = 0
    ) %>%
    as.data.frame()

  rownames(score_df) <-
    score_df$module

  score_df$module <- NULL

  score_mat <-
    as.matrix(score_df)


  # -----------------------------------------------------------------------
  # Check module IDs
  # -----------------------------------------------------------------------

  missing_modules <- setdiff(
    unique(
      c(
        edge_df$module_a,
        edge_df$module_b
      )
    ),
    rownames(score_mat)
  )

  if (length(missing_modules) > 0) {
    stop(
      "Missing hallmark loadings for modules in ",
      cancer_name,
      ": ",
      paste(
        missing_modules,
        collapse = ", "
      )
    )
  }


  # -----------------------------------------------------------------------
  # Observed hallmark-pair rewiring scores
  # -----------------------------------------------------------------------

  obs_stats <- compute_pair_rank_stats(
    score_mat = score_mat,
    edge_df = edge_df,
    term_pairs = term_pairs
  )


  # -----------------------------------------------------------------------
  # Permutation test
  #
  # Hallmark-loading profiles are reassigned among modules.
  # Module-pair rewiring scores remain fixed.
  #
  # The original cancer-specific analysis used set.seed(1).
  # Reset it for each cancer to preserve that behavior.
  # -----------------------------------------------------------------------

  set.seed(1)

  perm_stats_list <- vector(
    "list",
    n_perm
  )

  for (i in seq_len(n_perm)) {

    score_perm <- score_mat

    rownames(score_perm) <- sample(
      rownames(score_perm)
    )

    perm_stats_list[[i]] <-
      compute_pair_rank_stats(
        score_mat = score_perm,
        edge_df = edge_df,
        term_pairs = term_pairs
      ) %>%

      select(
        term_1,
        term_2,
        perm_weighted_rank_centered =
          weighted_rank_centered
      ) %>%

      mutate(
        perm_id = i
      )
  }

  perm_stats <- bind_rows(
    perm_stats_list
  )


  # -----------------------------------------------------------------------
  # Permutation distribution
  # -----------------------------------------------------------------------

  perm_summary <- perm_stats %>%
    group_by(
      term_1,
      term_2
    ) %>%
    summarise(
      perm_mean = mean(
        perm_weighted_rank_centered,
        na.rm = TRUE
      ),

      perm_sd = sd(
        perm_weighted_rank_centered,
        na.rm = TRUE
      ),

      .groups = "drop"
    )


  # -----------------------------------------------------------------------
  # One-sided empirical P values
  #
  # Lower tail: weakening
  # Upper tail: strengthening
  # -----------------------------------------------------------------------

  empirical_p <- perm_stats %>%

    left_join(
      obs_stats %>%
        select(
          term_1,
          term_2,
          observed_score =
            weighted_rank_centered
        ),

      by = c(
        "term_1",
        "term_2"
      )
    ) %>%

    group_by(
      term_1,
      term_2
    ) %>%

    summarise(
      p_emp_weaken = mean(
        perm_weighted_rank_centered <=
          observed_score,
        na.rm = TRUE
      ),

      p_emp_strengthen = mean(
        perm_weighted_rank_centered >=
          observed_score,
        na.rm = TRUE
      ),

      .groups = "drop"
    )


  # -----------------------------------------------------------------------
  # Combine results
  # -----------------------------------------------------------------------

  obs_stats %>%

    left_join(
      perm_summary,
      by = c(
        "term_1",
        "term_2"
      )
    ) %>%

    left_join(
      empirical_p,
      by = c(
        "term_1",
        "term_2"
      )
    ) %>%

    mutate(
      cancer = cancer_name,

      z_weaken =
        (
          weighted_rank_centered -
            perm_mean
        ) /
        perm_sd,

      fdr_weaken = p.adjust(
        p_emp_weaken,
        method = "BH"
      ),

      fdr_strengthen = p.adjust(
        p_emp_strengthen,
        method = "BH"
      )
    ) %>%

    select(
      cancer,
      term_1,
      term_2,
      n_edges_with_weight,
      weight_sum,
      weighted_raw_frac_change,
      weighted_rank,
      weighted_rank_centered,
      weighted_frac_relatively_weakened,
      perm_mean,
      perm_sd,
      z_weaken,
      p_emp_weaken,
      p_emp_strengthen,
      fdr_weaken,
      fdr_strengthen
    )
}


# -------------------------------------------------------------------------
# 11. Run all seven cancer types
# -------------------------------------------------------------------------

hallmark_pair_rewiring <- map_dfr(
  cancer_types,
  run_hallmark_rewiring
) %>%
  mutate(
    cancer = factor(
      cancer,
      levels = cancer_types
    )
  ) %>%
  arrange(
    cancer,
    term_1,
    term_2
  ) %>%
  mutate(
    cancer = as.character(cancer)
  )


# -------------------------------------------------------------------------
# 12. Significant cancer-specific hallmark-pair rewiring
#
# Nominal direction-specific one-sided empirical P < 0.05.
# -------------------------------------------------------------------------

hallmark_pair_rewiring_significant <-
  hallmark_pair_rewiring %>%

  mutate(
    rewiring = case_when(

      weighted_rank_centered < 0 &
        p_emp_weaken < 0.05 ~
        "weakening",

      weighted_rank_centered > 0 &
        p_emp_strengthen < 0.05 ~
        "strengthening",

      TRUE ~
        NA_character_
    ),

    empirical_p = case_when(

      rewiring == "weakening" ~
        p_emp_weaken,

      rewiring == "strengthening" ~
        p_emp_strengthen,

      TRUE ~
        NA_real_
    )
  ) %>%

  filter(
    !is.na(rewiring)
  ) %>%

  mutate(
    cancer = factor(
      cancer,
      levels = cancer_types
    ),

    rewiring = factor(
      rewiring,
      levels = c(
        "strengthening",
        "weakening"
      )
    )
  ) %>%

  arrange(
    cancer,
    rewiring,
    weighted_rank_centered
  ) %>%

  transmute(
    cancer = as.character(cancer),
    term_1,
    term_2,
    edges = n_edges_with_weight,
    score = weighted_rank_centered,
    empirical_p,
    rewiring = as.character(rewiring)
  )


# -------------------------------------------------------------------------
# 13. Save results
# -------------------------------------------------------------------------

fwrite(
  module_hallmark_loadings,
  file = loading_output_file,
  sep = "\t",
  quote = FALSE
)

fwrite(
  hallmark_pair_rewiring,
  file = pair_output_file,
  sep = "\t",
  quote = FALSE
)

fwrite(
  hallmark_pair_rewiring_significant,
  file = significant_output_file,
  sep = "\t",
  quote = FALSE
)


# -------------------------------------------------------------------------
# 14. Reproducibility summary
# -------------------------------------------------------------------------

message(
  "Retained modules: ",
  nrow(module_meta)
)

message(
  "Module-hallmark rows: ",
  nrow(module_hallmark_loadings)
)

message(
  "Significant cancer-specific hallmark-pair events: ",
  nrow(hallmark_pair_rewiring_significant)
)

message(
  "Hallmark rewiring analysis completed."
)
