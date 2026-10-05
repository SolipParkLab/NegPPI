# 04_module_architecture.R
#
# Define cancer-type-specific protein modules from positive co-abundance
# associations in adjacent normal tissues and classify significant normal
# and tumor associations as intra-module or inter-module.
#
# IMPORTANT:
#   Modules are defined ONLY from the normal positive-association network.
#   Tumor-derived modules are not constructed.
#
# Input:
#   results/normal_tumor/
#     exprs.cor.all_genes.partners.wt_only.<cancer>.normal.txt
#     exprs.cor.all_genes.partners.wt_only.<cancer>.tumor.txt
#
# Output:
#   results/modules/module_assignments.tsv
#   results/modules/association_module_classification.tsv


# -------------------------------------------------------------------------
# 1. Packages
# -------------------------------------------------------------------------

library(data.table)
library(dplyr)
library(igraph)


# -------------------------------------------------------------------------
# 2. Parameters
# -------------------------------------------------------------------------

r_cut <- 0.5
fdr_cut <- 0.05
resolution_use <- 1.5
k_use <- 50

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
# 3. Paths
# -------------------------------------------------------------------------

input_dir <- file.path(
  "results",
  "normal_tumor"
)

output_dir <- file.path(
  "results",
  "modules"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

module_output_file <- file.path(
  output_dir,
  "module_assignments.tsv"
)

association_output_file <- file.path(
  output_dir,
  "association_module_classification.tsv"
)


# -------------------------------------------------------------------------
# 4. Helper: classify associations
# -------------------------------------------------------------------------

classify_state <- function(
    r,
    p_adj
) {

  case_when(

    r >= r_cut &
      p_adj < fdr_cut ~ "Positive",

    r <= -r_cut &
      p_adj < fdr_cut ~ "Negative",

    TRUE ~ "n.s."
  )
}


# -------------------------------------------------------------------------
# 5. Helper: construct significant edge table
# -------------------------------------------------------------------------

make_edge_table <- function(
    merged_df,
    tissue = c(
      "normal",
      "tumor"
    )
) {

  tissue <- match.arg(tissue)

  if (tissue == "normal") {

    state_col <- "normal_state"
    r_col <- "wt_cor_r_normal"
    p_col <- "wt_cor_p_adj_normal"

  } else {

    state_col <- "tumor_state"
    r_col <- "wt_cor_r_tumor"
    p_col <- "wt_cor_p_adj_tumor"
  }

  merged_df %>%
    transmute(
      protein1 = node1,
      protein2 = node2,
      sign = .data[[state_col]],
      r = .data[[r_col]],
      p_adj = .data[[p_col]]
    ) %>%
    filter(
      sign != "n.s."
    )
}

# -------------------------------------------------------------------------
# 6. Helper: assign module class to associations
# -------------------------------------------------------------------------

annotate_module_class <- function(
    edge_df,
    module_tbl
) {

  edge_df %>%

    left_join(

      module_tbl %>%
        rename(
          module1 = module
        ),

      by = c(
        "protein1" = "protein"
      )

    ) %>%

    left_join(

      module_tbl %>%
        rename(
          module2 = module
        ),

      by = c(
        "protein2" = "protein"
      )

    ) %>%

    # Only associations for which both proteins are represented
    # in the normal-derived modular architecture are retained.

    filter(
      !is.na(module1),
      !is.na(module2)
    ) %>%

    mutate(

      association_class = ifelse(
        module1 == module2,
        "Intra-module",
        "Inter-module"
      )
    )
}


# -------------------------------------------------------------------------
# 7. Analysis for one cancer type
# -------------------------------------------------------------------------

run_module_analysis <- function(
    cancer
) {

  message(
    "\nProcessing ",
    cancer,
    "..."
  )


  # -----------------------------------------------------------------------
  # Input files
  # -----------------------------------------------------------------------

  cancer_file_label <- tolower(
    cancer
  )

  normal_file <- file.path(
    input_dir,
    paste0(
      "exprs.cor.all_genes.partners.wt_only.",
      cancer_file_label,
      ".normal.txt"
    )
  )

  tumor_file <- file.path(
    input_dir,
    paste0(
      "exprs.cor.all_genes.partners.wt_only.",
      cancer_file_label,
      ".tumor.txt"
    )
  )


  # -----------------------------------------------------------------------
  # Load normal and tumor correlation results
  # -----------------------------------------------------------------------

  res_normal <- fread(
    normal_file,
    header = TRUE
  )

  res_tumor <- fread(
    tumor_file,
    header = TRUE
  )


  # -----------------------------------------------------------------------
  # Restrict analysis to query-target pairs tested in both tissues
  # -----------------------------------------------------------------------

  merged_df <- inner_join(

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

      normal_state = classify_state(
        wt_cor_r_normal,
        wt_cor_p_adj_normal
      ),

      tumor_state = classify_state(
        wt_cor_r_tumor,
        wt_cor_p_adj_tumor
      ),

      # Convert directional query-target pairs into
      # undirected protein pairs.

      node1 = pmin(
        query,
        target
      ),

      node2 = pmax(
        query,
        target
      )

    ) %>%

    # Preserve the original analysis behavior:
    # only one record is retained for each undirected pair.

    distinct(
      node1,
      node2,
      .keep_all = TRUE
    )


  # -----------------------------------------------------------------------
  # Define normal positive associations
  # -----------------------------------------------------------------------

  normal_positive_edges <- merged_df %>%

    filter(
      normal_state == "Positive"
    ) %>%

    transmute(

      from = node1,

      to = node2,

      r = wt_cor_r_normal,

      p_adj = wt_cor_p_adj_normal,

      weight = wt_cor_r_normal
    )


  if (nrow(normal_positive_edges) == 0) {

    stop(
      "No normal positive associations found for ",
      cancer
    )
  }


  # -----------------------------------------------------------------------
  # Top-k sparsification
  #
  # For each protein, retain up to the 50 strongest positive partners.
  # -----------------------------------------------------------------------

  normal_positive_topk_long <- bind_rows(

    normal_positive_edges %>%

      transmute(
        node = from,
        partner = to,
        r = r,
        p_adj = p_adj,
        weight = weight
      ),

    normal_positive_edges %>%

      transmute(
        node = to,
        partner = from,
        r = r,
        p_adj = p_adj,
        weight = weight
      )

  ) %>%

    group_by(
      node
    ) %>%

    slice_max(
      order_by = weight,
      n = k_use,
      with_ties = FALSE
    ) %>%

    ungroup()


  # Convert back to unique undirected edges.

  normal_positive_topk_edges <-
    normal_positive_topk_long %>%

    transmute(

      from = pmin(
        node,
        partner
      ),

      to = pmax(
        node,
        partner
      ),

      r = r,

      p_adj = p_adj,

      weight = weight

    ) %>%

    distinct(
      from,
      to,
      .keep_all = TRUE
    )


  # -----------------------------------------------------------------------
  # Construct normal positive-association network
  # -----------------------------------------------------------------------

  g_normal <- graph_from_data_frame(

    normal_positive_topk_edges,

    directed = FALSE
  )


  # -----------------------------------------------------------------------
  # Louvain module detection
  # -----------------------------------------------------------------------

  set.seed(1)

  normal_modules <- cluster_louvain(

    g_normal,

    weights = E(g_normal)$weight,

    resolution = resolution_use
  )


  module_tbl <- data.frame(

    cancer = cancer,

    protein = V(g_normal)$name,

    module = as.integer(
      membership(
        normal_modules
      )
    ),

    stringsAsFactors = FALSE
  )


  module_sizes <- module_tbl %>%

    count(
      module,
      name = "module_size"
    )


  module_tbl <- module_tbl %>%

    left_join(
      module_sizes,
      by = "module"
    )


  message(
    "Normal-derived modules: ",
    length(
      unique(
        module_tbl$module
      )
    )
  )

  message(
    "Proteins assigned to modules: ",
    nrow(
      module_tbl
    )
  )


  # -----------------------------------------------------------------------
  # Build significant normal and tumor association tables
  # -----------------------------------------------------------------------

  normal_edges <- make_edge_table(
    merged_df,
    tissue = "normal"
  )

  tumor_edges <- make_edge_table(
    merged_df,
    tissue = "tumor"
  )


  # -----------------------------------------------------------------------
  # Project BOTH tissues onto the SAME normal-derived modules
  # -----------------------------------------------------------------------

  normal_annotated <- annotate_module_class(

    normal_edges,

    module_tbl %>%
      select(
        protein,
        module
      )
  ) %>%

    mutate(
      cancer = cancer,
      tissue = "Normal"
    )


  tumor_annotated <- annotate_module_class(

    tumor_edges,

    module_tbl %>%
      select(
        protein,
        module
      )
  ) %>%

    mutate(
      cancer = cancer,
      tissue = "Tumor"
    )


  association_tbl <- bind_rows(

    normal_annotated,

    tumor_annotated

  ) %>%

    select(

      cancer,

      tissue,

      protein1,

      protein2,

      sign,

      r,

      p_adj,

      module1,

      module2,

      association_class
    )


  # -----------------------------------------------------------------------
  # Return
  # -----------------------------------------------------------------------

  list(

    modules = module_tbl,

    associations = association_tbl
  )
}


# -------------------------------------------------------------------------
# 8. Run all seven cancer types
# -------------------------------------------------------------------------

all_results <- lapply(

  cancer_types,

  run_module_analysis
)

names(
  all_results
) <- cancer_types


# -------------------------------------------------------------------------
# 9. Combine module assignments
# -------------------------------------------------------------------------

module_assignments <- bind_rows(

  lapply(
    all_results,
    function(x) {
      x$modules
    }
  )
)


# -------------------------------------------------------------------------
# 10. Combine association classifications
# -------------------------------------------------------------------------

association_module_classification <- bind_rows(

  lapply(
    all_results,
    function(x) {
      x$associations
    }
  )
)


# -------------------------------------------------------------------------
# 11. Save results
# -------------------------------------------------------------------------

fwrite(

  module_assignments,

  file = module_output_file,

  sep = "\t",

  quote = FALSE
)


fwrite(

  association_module_classification,

  file = association_output_file,

  sep = "\t",

  quote = FALSE
)


message(
  "\nModule assignments written to: ",
  module_output_file
)

message(
  "Association classifications written to: ",
  association_output_file
)

message(
  "\nModule architecture analysis completed."
)
