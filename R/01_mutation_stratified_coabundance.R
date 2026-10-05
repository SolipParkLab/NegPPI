# 01_mutation_stratified_coabundance.R
#
# Construct query-specific wild-type baseline co-abundance correlations
# and evaluate mutation-associated changes in protein co-abundance.
#
# Input:
#   data/raw/AlphaMissense_aa_substitutions.tsv
#   data/raw/PanCan_Union_Maf_Broad_WashU_v1.1.maf
#   data/raw/Proteome_Broad_Institute_harmonized_v1.tsv
#   data/reference/p2s.rds
#
# Output:
#   results/exprs.cor.lm.all_genes.partners.wt_only.am_class.txt


# -------------------------------------------------------------------------
# 1. Packages
# -------------------------------------------------------------------------

library(data.table)
library(dplyr)
library(pbmcapply)


# -------------------------------------------------------------------------
# 2. Input and output paths
# -------------------------------------------------------------------------

alpha_missense_file <- file.path(
  "data",
  "raw",
  "AlphaMissense_aa_substitutions.tsv"
)

mutation_file <- file.path(
  "data",
  "raw",
  "PanCan_Union_Maf_Broad_WashU_v1.1.maf"
)

proteome_file <- file.path(
  "data",
  "raw",
  "Proteome_Broad_Institute_harmonized_v1.tsv"
)

mapping_file <- file.path(
  "data",
  "reference",
  "p2s.rds"
)

result_dir <- file.path(
  "results",
  "tmp_mutation_stratified_coabundance"
)

output_file <- file.path(
  "results",
  "exprs.cor.lm.all_genes.partners.wt_only.am_class.txt"
)

dir.create(
  result_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "results",
  recursive = TRUE,
  showWarnings = FALSE
)


# -------------------------------------------------------------------------
# 3. Load AlphaMissense predictions
# -------------------------------------------------------------------------

message("Loading AlphaMissense predictions...")

am <- fread(
  alpha_missense_file
)


# -------------------------------------------------------------------------
# 4. Load CPTAC somatic mutation data
# -------------------------------------------------------------------------

message("Loading CPTAC somatic mutation data...")

res <- fread(
  mutation_file,
  header = TRUE,
  data.table = TRUE
)

setDT(res)


# -------------------------------------------------------------------------
# 5. Load RefSeq-to-HGNC mapping
# -------------------------------------------------------------------------

message("Loading RefSeq-to-HGNC mapping...")

p2s <- readRDS(
  mapping_file
)


# -------------------------------------------------------------------------
# 6. Load CPTAC harmonized proteome
# -------------------------------------------------------------------------

message("Loading CPTAC harmonized proteome...")

expr <- fread(
  proteome_file,
  header = TRUE,
  data.table = FALSE
)

rownames(expr) <- expr[[1]]
expr <- expr[, -1]

# Harmonize proteome sample identifiers with mutation sample identifiers.
colnames(expr) <- sub(
  "-T$",
  "",
  gsub("\\.", "-", colnames(expr))
)


# -------------------------------------------------------------------------
# 7. Match proteome and mutation samples
# -------------------------------------------------------------------------

smpl_res <- sub(
  "_T",
  "",
  unique(res$Tumor_Sample_Barcode)
)

smpl_expr <- colnames(expr)

# Direct sample matches.
smpl_expr_res_1 <- intersect(
  smpl_expr,
  smpl_res
)

# Remaining proteome samples.
smpl_expr_rmn_1 <- setdiff(
  smpl_expr,
  smpl_res
)

# Remaining mutation samples.
smpl_res_rmn_1 <- setdiff(
  smpl_res,
  smpl_expr_res_1
)

# Some proteome identifiers contain an additional leading alphabetic
# character. Match these after removing the leading character.
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

expr <- expr[, smpl_expr_res]


# -------------------------------------------------------------------------
# 8. Annotate missense mutations using AlphaMissense
# -------------------------------------------------------------------------

message("Annotating missense mutations with AlphaMissense...")

res_ms <- subset(
  res,
  Variant_Classification == "Missense_Mutation"
)

res_ms <- as.data.table(
  res_ms
)

res_ms[
  ,
  aa_change := gsub(
    "^p\\.",
    "",
    Protein_Change
  )
]

setkeyv(
  res_ms,
  c(
    "HGNC_UniProt_ID(supplied_by_UniProt)",
    "aa_change"
  )
)

setkeyv(
  am,
  c(
    "uniprot_id",
    "protein_variant"
  )
)

res_annot <- am[
  res_ms,
  nomatch = 0
]


# -------------------------------------------------------------------------
# 9. Define query and target protein sets
# -------------------------------------------------------------------------

# Query genes:
# genes with missense variants successfully annotated by AlphaMissense.
genes <- unique(
  res_annot$Hugo_Symbol
)

# RefSeq peptide identifiers represented in the proteome dataset.
p_info <- unique(
  sub(
    "\\..{1,2}",
    "",
    rownames(expr)
  )
)

# HGNC symbols represented in the proteome dataset.
s_info <- unique(
  subset(
    p2s,
    refseq_peptide %in% p_info
  )$hgnc_symbol
)


# -------------------------------------------------------------------------
# 10. Query-specific mutation-stratified co-abundance analysis
# -------------------------------------------------------------------------

run_gene_analysis <- function(gene) {

  outfile <- file.path(
    result_dir,
    paste0(
      "res_",
      gene,
      ".tsv"
    )
  )

  # Skip genes that have already been analyzed.
  if (file.exists(outfile)) {
    return(NULL)
  }

  # AlphaMissense-annotated missense variants for the query gene.
  res_gene <- res_annot[
    Hugo_Symbol == gene
  ]

  if (nrow(res_gene) == 0) {
    return(NULL)
  }


  # -----------------------------------------------------------------------
  # Select query-protein RefSeq isoform
  # -----------------------------------------------------------------------

  p2s_intr <- subset(
    p2s,
    hgnc_symbol == gene &
      refseq_peptide != "" &
      !is.na(refseq_peptide)
  )

  ps_intr <- rownames(expr)[
    sub(
      "\\..{1,2}",
      "",
      rownames(expr)
    ) %in% p2s_intr$refseq_peptide
  ]

  if (length(ps_intr) == 0) {
    return(NULL)
  }

  # If multiple RefSeq isoforms map to the same query gene,
  # retain the isoform with the fewest missing measurements.
  if (length(ps_intr) > 1) {

    na_counts_intr <- rowSums(
      is.na(
        expr[
          ps_intr,
          ,
          drop = FALSE
        ]
      )
    )

    ps_intr <- ps_intr[
      which.min(
        na_counts_intr
      )
    ]
  }


  # -----------------------------------------------------------------------
  # Define query-specific wild-type samples
  # -----------------------------------------------------------------------

  # WT samples are tumors without any somatic mutation in the query gene,
  # rather than tumors lacking only the missense mutation class being tested.
  res_gene_bg <- res[
    Hugo_Symbol == gene
  ]

  smpl_mut_bg <- unique(
    sub(
      "_T",
      "",
      res_gene_bg$Tumor_Sample_Barcode
    )
  )

  smpl_expr <- colnames(expr)

  smpl_expr_rmn_bg_1 <- setdiff(
    smpl_expr,
    smpl_mut_bg
  )

  smpl_expr_mut_bg_1 <- intersect(
    smpl_expr,
    smpl_mut_bg
  )

  smpl_mut_rmn_bg_1 <- setdiff(
    smpl_mut_bg,
    smpl_expr_mut_bg_1
  )

  smpl_expr_mut_bg_2 <- smpl_expr_rmn_bg_1[
    gsub(
      "^([A-Za-z])",
      "",
      smpl_expr_rmn_bg_1
    ) %in% smpl_mut_rmn_bg_1
  ]

  smpl_expr_mut_bg <- sort(
    c(
      smpl_expr_mut_bg_1,
      smpl_expr_mut_bg_2
    )
  )

  smpl_expr_wt <- setdiff(
    smpl_expr,
    smpl_expr_mut_bg
  )

  if (length(smpl_expr_wt) < 5) {
    return(NULL)
  }


  # -----------------------------------------------------------------------
  # Analyze each AlphaMissense class separately
  # -----------------------------------------------------------------------

  local_result <- list()

  for (cls in unique(res_gene$am_class)) {

    res_cls <- res_gene[
      am_class == cls
    ]

    smpl_cls <- unique(
      sub(
        "_T",
        "",
        res_cls$Tumor_Sample_Barcode
      )
    )

    smpl_expr_cls_1 <- intersect(
      smpl_expr,
      smpl_cls
    )

    smpl_cls_rmn <- setdiff(
      smpl_cls,
      smpl_expr_cls_1
    )

    smpl_expr_cls_2 <- smpl_expr[
      gsub(
        "^([A-Za-z])",
        "",
        smpl_expr
      ) %in% smpl_cls_rmn
    ]

    smpl_expr_cls <- sort(
      c(
        smpl_expr_cls_1,
        smpl_expr_cls_2
      )
    )

    # Require at least five mutant samples in the mutation class.
    if (length(smpl_expr_cls) < 5) {
      next
    }


    # ---------------------------------------------------------------------
    # Test all proteome-mapped target proteins
    # ---------------------------------------------------------------------

    for (s_ppi in s_info) {

      p2s_ppi <- subset(
        p2s,
        hgnc_symbol == s_ppi &
          refseq_peptide != "" &
          !is.na(refseq_peptide)
      )

      ps_ppi <- rownames(expr)[
        sub(
          "\\..{1,2}",
          "",
          rownames(expr)
        ) %in% p2s_ppi$refseq_peptide
      ]

      if (length(ps_ppi) == 0) {
        next
      }

      # If multiple RefSeq isoforms map to a target gene,
      # retain the isoform with the fewest missing measurements
      # in the current mutant sample set.
      na_counts <- rowSums(
        is.na(
          expr[
            ps_ppi,
            smpl_expr_cls,
            drop = FALSE
          ]
        )
      )

      ps_ppi <- ps_ppi[
        which.min(
          na_counts
        )
      ]

      if (
        sum(
          !is.na(
            expr[
              ps_ppi,
              smpl_expr_cls
            ]
          )
        ) < 5
      ) {
        next
      }


      # -------------------------------------------------------------------
      # Construct mutant and wild-type expression data
      # -------------------------------------------------------------------

      expr_mut <- as.data.frame(
        t(
          expr[
            c(
              ps_intr,
              ps_ppi
            ),
            smpl_expr_cls
          ]
        )
      )

      expr_wt <- as.data.frame(
        t(
          expr[
            c(
              ps_intr,
              ps_ppi
            ),
            smpl_expr_wt
          ]
        )
      )

      colnames(expr_mut) <- c(
        "x",
        "y"
      )

      colnames(expr_wt) <- c(
        "x",
        "y"
      )

      # Require at least five pairwise-complete samples
      # in both WT and MUT groups.
      if (
        sum(
          complete.cases(
            expr_mut
          )
        ) < 5 ||
          sum(
            complete.cases(
              expr_wt
            )
          ) < 5
      ) {
        next
      }


      # -------------------------------------------------------------------
      # Pearson correlations
      # -------------------------------------------------------------------

      cor_mut <- suppressWarnings(
        cor.test(
          expr_mut$x,
          expr_mut$y,
          na.action = na.omit
        )
      )

      cor_wt <- suppressWarnings(
        cor.test(
          expr_wt$x,
          expr_wt$y,
          na.action = na.omit
        )
      )

      cor_df <- data.frame(
        wt_cor_r = cor_wt$estimate,
        wt_cor_p = cor_wt$p.value,
        mut_cor_r = cor_mut$estimate,
        mut_cor_p = cor_mut$p.value
      )


      # -------------------------------------------------------------------
      # Linear regression interaction model
      #
      # y ~ x * group
      #
      # group = 0: query-specific WT
      # group = 1: mutation-class-specific MUT
      # -------------------------------------------------------------------

      expr_mut$group <- 1
      expr_wt$group <- 0

      expr_cmb <- rbind(
        expr_mut,
        expr_wt
      )

      model <- lm(
        y ~ x * group,
        data = expr_cmb,
        na.action = na.omit
      )

      coef_mat <- summary(
        model
      )$coefficients

      coef_vec <- as.vector(
        t(
          coef_mat
        )
      )

      coef_df <- as.data.frame(
        t(
          coef_vec
        )
      )

      colnames(coef_df) <- gsub(
        ":",
        "_",
        gsub(
          "\\.",
          "",
          gsub(
            "[() >]",
            "",
            paste(
              rep(
                rownames(coef_mat),
                each = ncol(coef_mat)
              ),
              colnames(coef_mat),
              sep = "_"
            )
          )
        )
      )


      # -------------------------------------------------------------------
      # Store result
      # -------------------------------------------------------------------

      result <- data.frame(
        query = gene,
        query_refseq = paste(
          ps_intr,
          collapse = ";"
        ),
        mutation = cls,
        target = s_ppi,
        target_refseq = ps_ppi,
        wt_n = nrow(expr_wt),
        mut_n = nrow(expr_mut),
        cor_df,
        coef_df
      )

      local_result[[
        paste(
          gene,
          cls,
          s_ppi,
          sep = "_"
        )
      ]] <- result
    }
  }


  # -----------------------------------------------------------------------
  # Save one intermediate result per query gene
  # -----------------------------------------------------------------------

  if (length(local_result) > 0) {

    fwrite(
      do.call(
        rbind,
        local_result
      ),
      file = outfile,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }

  return(NULL)
}


# -------------------------------------------------------------------------
# 11. Run analysis in parallel
# -------------------------------------------------------------------------

num_cores <- max(
  1,
  parallel::detectCores() - 1
)

message(
  "Running ",
  length(genes),
  " query genes using ",
  num_cores,
  " cores..."
)

pbmclapply(
  genes,
  run_gene_analysis,
  mc.cores = num_cores
)


# -------------------------------------------------------------------------
# 12. Combine query-level result files
# -------------------------------------------------------------------------

message("Combining query-level results...")

result_files <- list.files(
  result_dir,
  full.names = TRUE,
  pattern = "\\.tsv$"
)

final_result <- rbindlist(
  lapply(
    result_files,
    fread
  )
)

# Remove self-correlations.
final_result <- subset(
  final_result,
  target != query
)


# -------------------------------------------------------------------------
# 13. Multiple-testing correction
# -------------------------------------------------------------------------

final_result$wt_cor_p_adj <- p.adjust(
  final_result$wt_cor_p,
  method = "fdr"
)

final_result$mut_cor_p_adj <- p.adjust(
  final_result$mut_cor_p,
  method = "fdr"
)

int_p_col <- grep(
  "x_group_Pr",
  colnames(final_result),
  value = TRUE
)

if (length(int_p_col) > 0) {

  final_result$int_p_adj <- p.adjust(
    final_result[[int_p_col[1]]],
    method = "fdr"
  )
}


# -------------------------------------------------------------------------
# 14. Save combined result
# -------------------------------------------------------------------------

write.table(
  final_result,
  file = output_file,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# -------------------------------------------------------------------------
# 15. Summary
# -------------------------------------------------------------------------

message(
  "Number of result rows: ",
  format(
    nrow(final_result),
    big.mark = ","
  )
)

message(
  "Number of query proteins: ",
  length(
    unique(
      final_result$query
    )
  )
)

message(
  "Output written to: ",
  output_file
)
