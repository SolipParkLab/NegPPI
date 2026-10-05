# 05_module_pair_rewiring.R
#
# Quantify normal-to-tumor rewiring between modules defined from
# positive co-abundance associations in adjacent normal tissues.
#
# Modules are NOT redefined in tumor tissue.
#
# Input:
#   results/modules/module_assignments.tsv
#   results/modules/association_module_classification.tsv
#
# Output:
#   results/module_rewiring/module_pair_rewiring.tsv


library(data.table)
library(dplyr)
library(tidyr)


# -------------------------------------------------------------------------
# 1. Parameters
# -------------------------------------------------------------------------

min_module_size <- 50

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

association_file <- file.path(
  "results",
  "modules",
  "association_module_classification.tsv"
)

output_dir <- file.path(
  "results",
  "module_rewiring"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

output_file <- file.path(
  output_dir,
  "module_pair_rewiring.tsv"
)


# -------------------------------------------------------------------------
# 3. Load data from script 04
# -------------------------------------------------------------------------

module_assignments <- fread(
  module_file,
  header = TRUE
) %>%
  as_tibble()

association_tbl <- fread(
  association_file,
  header = TRUE
) %>%
  as_tibble()


# -------------------------------------------------------------------------
# 4. Retain normal-derived modules containing at least 50 proteins
# -------------------------------------------------------------------------

retained_modules <- module_assignments %>%
  transmute(
    cancer = as.character(cancer),
    module = as.character(module),
    module_size = as.integer(module_size)
  ) %>%
  distinct() %>%
  filter(
    module_size >= min_module_size
  )


# -------------------------------------------------------------------------
# 5. Retain inter-module associations between retained modules
#
# Positive and negative associations are combined.
# Normal and tumor associations use the SAME normal-derived modules.
# -------------------------------------------------------------------------

inter_module_edges <- association_tbl %>%
  transmute(
    cancer = as.character(cancer),
    tissue = as.character(tissue),
    module1 = as.character(module1),
    module2 = as.character(module2),
    association_class = as.character(association_class)
  ) %>%
  filter(
    association_class == "Inter-module"
  ) %>%

  inner_join(
    retained_modules %>%
      select(
        cancer,
        module1 = module
      ),
    by = c(
      "cancer",
      "module1"
    )
  ) %>%

  inner_join(
    retained_modules %>%
      select(
        cancer,
        module2 = module
      ),
    by = c(
      "cancer",
      "module2"
    )
  ) %>%

  mutate(
    # Keep the original character-based module-pair ordering.
    module_a = pmin(module1, module2),
    module_b = pmax(module1, module2)
  )


# -------------------------------------------------------------------------
# 6. Count significant inter-module associations
#
# Positive and negative associations are counted together.
# -------------------------------------------------------------------------

module_pair_counts <- inter_module_edges %>%
  count(
    cancer,
    tissue,
    module_a,
    module_b,
    name = "n_associations"
  )


normal_counts <- module_pair_counts %>%
  filter(
    tissue == "Normal"
  ) %>%
  transmute(
    cancer,
    module_a,
    module_b,
    normal_associations_n = n_associations
  )


tumor_counts <- module_pair_counts %>%
  filter(
    tissue == "Tumor"
  ) %>%
  transmute(
    cancer,
    module_a,
    module_b,
    tumor_associations_n = n_associations
  )


# -------------------------------------------------------------------------
# 7. Calculate normal-to-tumor module-pair rewiring
#
# Only module pairs with at least one association in normal tissue
# are retained.
#
# Fractional change:
#
#   (Tumor - Normal) / Normal
# -------------------------------------------------------------------------

module_pair_rewiring <- normal_counts %>%

  left_join(
    tumor_counts,
    by = c(
      "cancer",
      "module_a",
      "module_b"
    )
  ) %>%

  mutate(
    tumor_associations_n = replace_na(
      tumor_associations_n,
      0L
    ),

    change_associations_n =
      tumor_associations_n - normal_associations_n,

    fractional_change =
      change_associations_n / normal_associations_n
  ) %>%

  group_by(
    cancer
  ) %>%

  mutate(
    within_cancer_rewiring_rank =
      rank(
        fractional_change,
        ties.method = "average"
      ) / (n() + 1),

    centered_within_cancer_rewiring_rank =
      within_cancer_rewiring_rank - 0.5
  ) %>%

  ungroup()


# -------------------------------------------------------------------------
# 8. Format module labels and order cancer types
# -------------------------------------------------------------------------

module_pair_rewiring <- module_pair_rewiring %>%

  mutate(
    cancer = factor(
      cancer,
      levels = cancer_types
    ),

    module_a = paste0(
      "M",
      module_a
    ),

    module_b = paste0(
      "M",
      module_b
    )
  ) %>%

  arrange(
    cancer,
    within_cancer_rewiring_rank
  ) %>%

  mutate(
    cancer = as.character(cancer)
  ) %>%

  select(
    cancer,
    module_a,
    module_b,
    normal_associations_n,
    tumor_associations_n,
    change_associations_n,
    fractional_change,
    within_cancer_rewiring_rank,
    centered_within_cancer_rewiring_rank
  )


# -------------------------------------------------------------------------
# 9. Reproducibility summary
# -------------------------------------------------------------------------

module_summary <- retained_modules %>%
  count(
    cancer,
    name = "n_retained_modules"
  ) %>%
  mutate(
    possible_module_pairs =
      choose(n_retained_modules, 2)
  )


pair_summary <- module_pair_rewiring %>%
  count(
    cancer,
    name = "n_module_pairs"
  ) %>%

  left_join(
    module_summary,
    by = "cancer"
  ) %>%

  mutate(
    all_possible_pairs_observed =
      n_module_pairs == possible_module_pairs
  )


print(
  pair_summary,
  n = Inf
)

message(
  "Retained modules: ",
  nrow(retained_modules)
)

message(
  "Module pairs with at least one normal association: ",
  nrow(module_pair_rewiring)
)


# -------------------------------------------------------------------------
# 10. Save
# -------------------------------------------------------------------------

fwrite(
  module_pair_rewiring,
  file = output_file,
  sep = "\t",
  quote = FALSE
)

message(
  "Module-pair rewiring results written to: ",
  output_file
)
