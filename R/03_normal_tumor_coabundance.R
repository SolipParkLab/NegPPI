# 03_normal_tumor_coabundance.R
#
# Calculate query-specific wild-type co-abundance correlations
# separately in adjacent normal and tumor tissues for seven CPTAC cohorts.
#
# Input:
#   data/raw/PanCan_Union_Maf_Broad_WashU_v1.1.maf
#   results/normal_tumor_proteomes.rds
#
# Output:
#   results/normal_tumor/
#     exprs.cor.all_genes.partners.wt_only.<cancer>.normal.txt
#     exprs.cor.all_genes.partners.wt_only.<cancer>.tumor.txt


# -------------------------------------------------------------------------
# 1. Packages
# -------------------------------------------------------------------------

library(data.table)
library(pbmcapply)


# -------------------------------------------------------------------------
# 2. Paths
# -------------------------------------------------------------------------

maf_file <- file.path(
  "data",
  "raw",
  "PanCan_Union_Maf_Broad_WashU_v1.1.maf"
)

proteome_file <- file.path(
  "results",
  "normal_tumor_proteomes.rds"
)

output_dir <- file.path(
  "results",
  "normal_tumor"
)

temp_root <- file.path(
  "results",
  "tmp_normal_tumor_coabundance"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  temp_root,
  recursive = TRUE,
  showWarnings = FALSE
)


# -------------------------------------------------------------------------
# 3. Load data
# -------------------------------------------------------------------------

message("Loading mutation data...")

maf_all <- fread(
  maf_file,
  header = TRUE,
  data.table = TRUE
)

message("Loading prepared normal/tumor proteomes...")

normal_tumor_data <- readRDS(
  proteome_file
)


# -------------------------------------------------------------------------
# 4. Cancer-type definitions
# -------------------------------------------------------------------------

# Display labels used in output.
cancer_types <- c(
  "BRCA",
  "ccRCC",
  "COAD",
  "LSCC",
  "LUAD",
  "PDAC",
  "UCEC"
)

# Cohort labels used in the pan-cancer MAF.
maf_cohort_labels <- c(
  BRCA  = "BRCA",
  ccRCC = "CCRCC",
  COAD  = "COAD",
  LSCC  = "LSCC",
  LUAD  = "LUAD",
  PDAC  = "PDAC",
  UCEC  = "UCEC"
)


# -------------------------------------------------------------------------
# 5. Helper: identify mutated samples in an expression matrix
# -------------------------------------------------------------------------

# The original analysis used two-step sample matching:
#
# 1. direct match between expression-matrix participant IDs and MAF IDs
# 2. for remaining samples, match after removing one leading
#    alphabetic character from the expression-matrix ID

find_mutated_samples <- function(
    expr_samples,
    mutation_samples
) {

  mutation_samples <- unique(
    mutation_samples
  )

  mutation_samples <- mutation_samples[
    !is.na(mutation_samples) &
      mutation_samples != ""
  ]

  # Direct matches.
  direct_matches <- intersect(
    expr_samples,
    mutation_samples
  )

  # Samples not matched directly.
  expr_remaining <- setdiff(
    expr_samples,
    direct_matches
  )

  mutation_remaining <- setdiff(
    mutation_samples,
    direct_matches
  )

  # Match after removing one leading alphabetic character.
  indirect_matches <- expr_remaining[
    gsub(
      "^([A-Za-z])",
      "",
      expr_remaining
    ) %in% mutation_remaining
  ]

  sort(
    unique(
      c(
        direct_matches,
        indirect_matches
      )
    )
  )
}


# -------------------------------------------------------------------------
# 6. Function for one cancer type and one tissue
# -------------------------------------------------------------------------

run_coabundance_analysis <- function(
    cancer,
    tissue,
    expr_matrix,
    maf
) {

  message(
    "\nRunning ",
    cancer,
    " - ",
    tissue,
    "..."
  )

  # -----------------------------------------------------------------------
  # Define query genes
  # -----------------------------------------------------------------------

  # Preserve the gene universe used in the original analysis:
  # genes represented in the cancer-specific MAF and detected
  # in the corresponding proteome matrix.

  genes <- unique(
    maf$Hugo_Symbol
  )

  genes <- genes[
    !is.na(genes) &
      genes != "" &
      genes %in% rownames(expr_matrix)
  ]

  message(
    "Number of query genes: ",
    length(genes)
  )


  # -----------------------------------------------------------------------
  # Tissue-specific mutation barcode column
  # -----------------------------------------------------------------------

  if (tissue == "tumor") {

    mutation_barcode_column <- "Tumor_Sample_Barcode"
    barcode_suffix <- "_T"

  } else if (tissue == "normal") {

    mutation_barcode_column <- "Matched_Norm_Sample_Barcode"
    barcode_suffix <- "_N"

  } else {

    stop(
      "tissue must be either 'normal' or 'tumor'."
    )
  }


  # -----------------------------------------------------------------------
  # Separate temporary directory for each cancer and tissue
  # -----------------------------------------------------------------------

  result_dir <- file.path(
    temp_root,
    cancer,
    tissue
  )

  dir.create(
    result_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )


  # -----------------------------------------------------------------------
  # Query-level analysis
  # -----------------------------------------------------------------------

  run_gene_analysis_wt <- function(gene) {

    outfile <- file.path(
      result_dir,
      paste0(
        "res_",
        gene,
        ".tsv"
      )
    )

    # Allow interrupted analyses to resume without recalculating
    # completed query genes.
    if (file.exists(outfile)) {
      return(NULL)
    }


    # ---------------------------------------------------------------------
    # Mutation samples for this query gene
    # ---------------------------------------------------------------------

    res_gene <- maf[
      Hugo_Symbol == gene
    ]

    mutation_samples <- unique(
      sub(
        barcode_suffix,
        "",
        res_gene[
          [mutation_barcode_column]
        ]
      )
    )

    expr_samples <- colnames(
      expr_matrix
    )

    mutated_expr_samples <- find_mutated_samples(
      expr_samples = expr_samples,
      mutation_samples = mutation_samples
    )


    # Query-specific WT samples are all expression samples
    # excluding samples corresponding to patients with a somatic
    # mutation in the query gene.

    wt_samples <- setdiff(
      expr_samples,
      mutated_expr_samples
    )

    if (length(wt_samples) < 5) {
      return(NULL)
    }


    # ---------------------------------------------------------------------
    # Target proteins
    # ---------------------------------------------------------------------

    targets <- setdiff(
      genes,
      gene
    )

    local_result <- vector(
      "list",
      length(targets)
    )


    # ---------------------------------------------------------------------
    # Query-target correlations
    # ---------------------------------------------------------------------

    for (i in seq_along(targets)) {

      target <- targets[i]


      # Require at least five non-missing target measurements
      # among query-specific WT samples.

      if (
        sum(
          !is.na(
            expr_matrix[
              target,
              wt_samples
            ]
          )
        ) < 5
      ) {
        next
      }


      pair_expr <- as.data.frame(
        t(
          expr_matrix[
            c(
              gene,
              target
            ),
            wt_samples,
            drop = FALSE
          ]
        )
      )

      colnames(pair_expr) <- c(
        "x",
        "y"
      )


      # Pairwise-complete observations.

      pair_expr <- pair_expr[
        complete.cases(pair_expr),
        ,
        drop = FALSE
      ]

      pairwise_n <- nrow(
        pair_expr
      )

      if (pairwise_n < 5) {
        next
      }


      # Pearson correlation.

      cor_wt <- suppressWarnings(
        cor.test(
          pair_expr$x,
          pair_expr$y
        )
      )


      local_result[[i]] <- data.frame(
        query = gene,
        target = target,
        wt_n = length(wt_samples),
        wt_pairwise_n = pairwise_n,
        wt_cor_r = unname(
          cor_wt$estimate
        ),
        wt_cor_p = cor_wt$p.value
      )
    }


    # Remove targets that did not meet sample-size criteria.

    local_result <- local_result[
      !vapply(
        local_result,
        is.null,
        logical(1)
      )
    ]


    # Save one intermediate file per query gene.

    if (length(local_result) > 0) {

      fwrite(
        rbindlist(
          local_result
        ),
        file = outfile,
        sep = "\t",
        quote = FALSE
      )
    }

    NULL
  }


  # -----------------------------------------------------------------------
  # Run query genes in parallel
  # -----------------------------------------------------------------------

  num_cores <- max(
    1,
    parallel::detectCores() - 1
  )

  pbmclapply(
    genes,
    run_gene_analysis_wt,
    mc.cores = num_cores
  )


  # -----------------------------------------------------------------------
  # Combine query-level files
  # -----------------------------------------------------------------------

  result_files <- list.files(
    result_dir,
    full.names = TRUE,
    pattern = "\\.tsv$"
  )

  if (length(result_files) == 0) {

    warning(
      "No results generated for ",
      cancer,
      " ",
      tissue
    )

    return(NULL)
  }

  final_result <- rbindlist(
    lapply(
      result_files,
      fread
    )
  )

  final_result <- final_result[
    target != query
  ]


  # -----------------------------------------------------------------------
  # FDR correction
  # -----------------------------------------------------------------------

  # Benjamini-Hochberg correction is performed separately
  # within each cancer type and tissue state, as in the
  # original analysis.

  final_result[
    ,
    wt_cor_p_adj :=
      p.adjust(
        wt_cor_p,
        method = "fdr"
      )
  ]


  # -----------------------------------------------------------------------
  # Association classification
  # -----------------------------------------------------------------------

  final_result[
    ,
    association := fifelse(
      wt_cor_r >= 0.5 &
        wt_cor_p_adj < 0.05,
      "Positive",
      fifelse(
        wt_cor_r <= -0.5 &
          wt_cor_p_adj < 0.05,
        "Negative",
        "n.s."
      )
    )
  ]


  # -----------------------------------------------------------------------
  # Save result
  # -----------------------------------------------------------------------

  cancer_file_label <- tolower(
    cancer
  )

  output_file <- file.path(
    output_dir,
    paste0(
      "exprs.cor.all_genes.partners.wt_only.",
      cancer_file_label,
      ".",
      tissue,
      ".txt"
    )
  )

  fwrite(
    final_result,
    file = output_file,
    sep = "\t",
    quote = FALSE
  )

  message(
    "Saved: ",
    output_file
  )

  message(
    "Number of tested query-target pairs: ",
    format(
      nrow(final_result),
      big.mark = ","
    )
  )

  message(
    "Positive associations: ",
    sum(
      final_result$association == "Positive"
    )
  )

  message(
    "Negative associations: ",
    sum(
      final_result$association == "Negative"
    )
  )


  final_result
}


# -------------------------------------------------------------------------
# 7. Run all seven cancer types
# -------------------------------------------------------------------------

all_results <- list()

for (cancer in cancer_types) {

  maf_label <- maf_cohort_labels[
    cancer
  ]

  maf <- maf_all[
    COHORT == maf_label
  ]


  # -----------------------------------------------------------------------
  # Normal
  # -----------------------------------------------------------------------

  normal_result <- run_coabundance_analysis(
    cancer = cancer,
    tissue = "normal",
    expr_matrix =
      normal_tumor_data[
        [cancer]
      ]$normal,
    maf = maf
  )


  # -----------------------------------------------------------------------
  # Tumor
  # -----------------------------------------------------------------------

  tumor_result <- run_coabundance_analysis(
    cancer = cancer,
    tissue = "tumor",
    expr_matrix =
      normal_tumor_data[
        [cancer]
      ]$tumor,
    maf = maf
  )


  all_results[
    [cancer]
  ] <- list(
    normal = normal_result,
    tumor = tumor_result
  )
}


# -------------------------------------------------------------------------
# 8. Save combined R object
# -------------------------------------------------------------------------

combined_output_file <- file.path(
  output_dir,
  "normal_tumor_coabundance_all_cancers.rds"
)

saveRDS(
  all_results,
  combined_output_file
)

message(
  "\nCombined normal/tumor results written to: ",
  combined_output_file
)

message(
  "\nNormal-tumor co-abundance analysis completed."
)
