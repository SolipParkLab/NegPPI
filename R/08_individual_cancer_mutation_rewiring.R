# 08_individual_cancer_mutation_rewiring.R
#
# Recalculate mutation-stratified co-abundance independently within each
# CPTAC cancer type, classify mutation-associated rewiring, and calculate
# WT/MUT pairwise-complete sample counts for baseline positive/negative edges.
#
# IMPORTANT:
# Cancer-specific results are recalculated from mutation and proteome data.
# They are NOT obtained by splitting the pan-cancer results from script 01.
#
# Main output used by script 09:
#   results/mutation_rewiring/individual/
#     individual_cancer_mutation_rewiring_edges.tsv


library(data.table)
library(dplyr)
library(purrr)
library(pbmcapply)


# -------------------------------------------------------------------------
# 1. Parameters
# -------------------------------------------------------------------------

r_cut <- 0.5
fdr_cut <- 0.05
interaction_p_cut <- 0.05

num_cores <- max(
  1,
  parallel::detectCores() - 1
)

cancer_types <- c(
  "brca",
  "ccrcc",
  "coad",
  "gbm",
  "hnscc",
  "lscc",
  "luad",
  "ov",
  "pdac",
  "ucec"
)

cancer_labels <- c(
  brca = "BRCA",
  ccrcc = "ccRCC",
  coad = "COAD",
  gbm = "GBM",
  hnscc = "HNSCC",
  lscc = "LSCC",
  luad = "LUAD",
  ov = "OV",
  pdac = "PDAC",
  ucec = "UCEC"
)

maf_cohorts <- c(
  brca = "BRCA",
  ccrcc = "CCRCC",
  coad = "COAD",
  gbm = "GBM",
  hnscc = "HNSCC",
  lscc = "LSCC",
  luad = "LUAD",
  ov = "OV",
  pdac = "PDAC",
  ucec = "UCEC"
)


# -------------------------------------------------------------------------
# 2. Paths
# -------------------------------------------------------------------------

alphamissense_file <- file.path(
  "data",
  "raw",
  "AlphaMissense_aa_substitutions.tsv"
)

maf_file <- file.path(
  "data",
  "raw",
  "PanCan_Union_Maf_Broad_WashU_v1.1.maf"
)

proteome_file <- file.path(
  "data",
  "raw",
  "Proteome_Broad_Institute_harmonized_v1.tsv"
)

p2s_file <- file.path(
  "data",
  "reference",
  "p2s.rds"
)

output_dir <- file.path(
  "results",
  "mutation_rewiring",
  "individual"
)

raw_output_dir <- file.path(
  output_dir,
  "raw"
)

tmp_root <- file.path(
  "results",
  "tmp_individual_mutation_rewiring"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  raw_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  tmp_root,
  recursive = TRUE,
  showWarnings = FALSE
)


# -------------------------------------------------------------------------
# 3. Load input data
# -------------------------------------------------------------------------

am <- fread(
  alphamissense_file
)

maf <- fread(
  maf_file,
  header = TRUE,
  data.table = FALSE
)

p2s <- readRDS(
  p2s_file
)

expr <- fread(
  proteome_file,
  header = TRUE,
  data.table = FALSE
)

rownames(expr) <- expr[[1]]

expr <- expr[
  ,
  -1,
  drop = FALSE
]

colnames(expr) <- sub(
  "-T$",
  "",
  gsub(
    "\\.",
    "-",
    colnames(expr)
  )
)

expr_samples <- colnames(
  expr
)


# -------------------------------------------------------------------------
# 4. Match mutation IDs to proteome sample IDs
# -------------------------------------------------------------------------

match_mutation_samples <- function(
    expression_ids,
    mutation_ids,
    leading_pool = NULL
) {

  mutation_ids <- unique(
    mutation_ids[
      !is.na(mutation_ids) &
        mutation_ids != ""
    ]
  )

  exact <- intersect(
    expression_ids,
    mutation_ids
  )

  mutation_remaining <- setdiff(
    mutation_ids,
    exact
  )

  if (is.null(leading_pool)) {

    leading_pool <- setdiff(
      expression_ids,
      mutation_ids
    )
  }

  leading <- leading_pool[
    gsub(
      "^([A-Za-z])",
      "",
      leading_pool
    ) %in%
      mutation_remaining
  ]

  sort(
    unique(
      c(
        exact,
        leading
      )
    )
  )
}


# =========================================================================
# 5. Recalculate WT/MUT co-abundance separately in each cancer type
# =========================================================================

for (cancer in cancer_types) {

  cancer_label <- unname(
    cancer_labels[[cancer]]
  )

  maf_cohort <- unname(
    maf_cohorts[[cancer]]
  )

  message(
    "Processing cancer-specific analysis: ",
    cancer_label
  )


  # -----------------------------------------------------------------------
  # 5.1 Cancer-specific mutation data
  # -----------------------------------------------------------------------

  res_cohort <- as.data.table(
    maf[
      maf$COHORT == maf_cohort,
      ,
      drop = FALSE
    ]
  )

  res_ms_cohort <- res_cohort[
    Variant_Classification == "Missense_Mutation"
  ]

  res_ms_cohort[
    ,
    aa_change := gsub(
      "^p\\.",
      "",
      Protein_Change
    )
  ]

  res_annot_cohort <- am[
    res_ms_cohort,
    on = .(
      uniprot_id =
        `HGNC_UniProt_ID(supplied_by_UniProt)`,
      protein_variant =
        aa_change
    ),
    nomatch = 0
  ]


  # -----------------------------------------------------------------------
  # 5.2 Cancer-specific proteome
  # -----------------------------------------------------------------------

  smpl_res_cohort <- unique(
    sub(
      "_T",
      "",
      res_cohort$Tumor_Sample_Barcode
    )
  )

  smpl_expr_cohort <- match_mutation_samples(
    expression_ids = expr_samples,
    mutation_ids = smpl_res_cohort
  )

  expr_cohort <- expr[
    ,
    smpl_expr_cohort,
    drop = FALSE
  ]


  # -----------------------------------------------------------------------
  # 5.3 Query and target universe
  # -----------------------------------------------------------------------

  genes_cohort <- unique(
    res_annot_cohort$Hugo_Symbol
  )

  genes_cohort <- genes_cohort[
    !is.na(genes_cohort) &
      genes_cohort != ""
  ]

  p_info <- unique(
    sub(
      "\\..{1,2}",
      "",
      rownames(expr_cohort)
    )
  )

  s_info <- unique(
    subset(
      p2s,
      refseq_peptide %in% p_info
    )$hgnc_symbol
  )

  s_info <- s_info[
    !is.na(s_info) &
      s_info != ""
  ]


  # -----------------------------------------------------------------------
  # 5.4 Resumable temporary directory
  # -----------------------------------------------------------------------

  tmp_dir <- file.path(
    tmp_root,
    cancer
  )

  dir.create(
    tmp_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )


  # -----------------------------------------------------------------------
  # 5.5 Analyze one query gene
  # -----------------------------------------------------------------------

  run_gene_analysis <- function(gene) {

    outfile <- file.path(
      tmp_dir,
      paste0(
        "res_",
        gene,
        ".tsv"
      )
    )

    if (file.exists(outfile)) {
      return(NULL)
    }

    res_gene <- res_annot_cohort[
      Hugo_Symbol == gene
    ]

    if (nrow(res_gene) == 0) {
      return(NULL)
    }


    # ---------------------------------------------------------------------
    # Query RefSeq
    # ---------------------------------------------------------------------

    p2s_intr <- subset(
      p2s,
      hgnc_symbol == gene &
        refseq_peptide != "" &
        !is.na(refseq_peptide)
    )

    ps_intr <- rownames(expr_cohort)[
      sub(
        "\\..{1,2}",
        "",
        rownames(expr_cohort)
      ) %in%
        p2s_intr$refseq_peptide
    ]

    if (length(ps_intr) == 0) {
      return(NULL)
    }

    if (length(ps_intr) > 1) {

      na_counts_intr <- rowSums(
        is.na(
          expr_cohort[
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


    # ---------------------------------------------------------------------
    # AlphaMissense-annotated missense samples
    # ---------------------------------------------------------------------

    smpl_expr <- colnames(
      expr_cohort
    )

    smpl_ms <- unique(
      sub(
        "_T",
        "",
        res_gene$Tumor_Sample_Barcode
      )
    )

    if (length(smpl_ms) < 5) {
      return(NULL)
    }

    # Preserve the original leading-letter matching pool.
    smpl_expr_rmn_annot <- setdiff(
      smpl_expr,
      smpl_ms
    )

    smpl_expr_ms <- match_mutation_samples(
      expression_ids = smpl_expr,
      mutation_ids = smpl_ms,
      leading_pool = smpl_expr_rmn_annot
    )

    if (length(smpl_expr_ms) < 5) {
      return(NULL)
    }


    # ---------------------------------------------------------------------
    # WT samples
    #
    # WT = no somatic mutation of any class in the query gene
    # ---------------------------------------------------------------------

    res_gene_bg <- res_cohort[
      Hugo_Symbol == gene
    ]

    smpl_mut_bg <- unique(
      sub(
        "_T",
        "",
        res_gene_bg$Tumor_Sample_Barcode
      )
    )

    smpl_expr_mut_bg <- match_mutation_samples(
      expression_ids = smpl_expr,
      mutation_ids = smpl_mut_bg
    )

    smpl_expr_wt <- setdiff(
      smpl_expr,
      smpl_expr_mut_bg
    )

    if (length(smpl_expr_wt) < 5) {
      return(NULL)
    }


    # ---------------------------------------------------------------------
    # Mutation classes
    # ---------------------------------------------------------------------

    local_result <- list()

    for (cls in unique(res_gene$am_class)) {

      res_cls <- res_gene[
        am_class == cls
      ]

      smpl_ms_flt <- unique(
        sub(
          "_T",
          "",
          res_cls$Tumor_Sample_Barcode
        )
      )

      if (length(smpl_ms_flt) < 5) {
        next
      }

      smpl_expr_ms_flt <- match_mutation_samples(
        expression_ids = smpl_expr,
        mutation_ids = smpl_ms_flt,
        leading_pool = smpl_expr_rmn_annot
      )

      if (length(smpl_expr_ms_flt) < 5) {
        next
      }


      # -------------------------------------------------------------------
      # Targets
      # -------------------------------------------------------------------

      for (s_ppi in s_info) {

        p2s_ppi <- subset(
          p2s,
          hgnc_symbol == s_ppi &
            refseq_peptide != "" &
            !is.na(refseq_peptide)
        )

        ps_ppi <- rownames(expr_cohort)[
          sub(
            "\\..{1,2}",
            "",
            rownames(expr_cohort)
          ) %in%
            p2s_ppi$refseq_peptide
        ]

        if (length(ps_ppi) == 0) {
          next
        }

        # Target isoform selection is performed in the CURRENT
        # mutation-class sample set.
        na_counts <- rowSums(
          is.na(
            expr_cohort[
              ps_ppi,
              smpl_expr_ms_flt,
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
              expr_cohort[
                ps_ppi,
                smpl_expr_ms_flt
              ]
            )
          ) < 5
        ) {
          next
        }


        # -----------------------------------------------------------------
        # WT and MUT abundance data
        # -----------------------------------------------------------------

        expr_mut <- as.data.frame(
          t(
            expr_cohort[
              c(
                ps_intr,
                ps_ppi
              ),
              smpl_expr_ms_flt,
              drop = FALSE
            ]
          )
        )

        expr_wt <- as.data.frame(
          t(
            expr_cohort[
              c(
                ps_intr,
                ps_ppi
              ),
              smpl_expr_wt,
              drop = FALSE
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


        # -----------------------------------------------------------------
        # Correlations
        # -----------------------------------------------------------------

        expr_mut$group <- 1
        expr_wt$group <- 0

        expr_cmb <- rbind(
          expr_mut,
          expr_wt
        )

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
          wt_cor_r =
            unname(
              cor_wt$estimate
            ),

          wt_cor_p =
            cor_wt$p.value,

          mut_cor_r =
            unname(
              cor_mut$estimate
            ),

          mut_cor_p =
            cor_mut$p.value
        )


        # -----------------------------------------------------------------
        # Interaction model
        #
        # y ~ x * group
        # -----------------------------------------------------------------

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


        # -----------------------------------------------------------------
        # Store result
        # -----------------------------------------------------------------

        result <- data.frame(
          query = gene,

          query_refseq =
            paste(
              ps_intr,
              collapse = ";"
            ),

          mutation = cls,

          target = s_ppi,
          target_refseq = ps_ppi,

          wt_n =
            nrow(
              expr_wt
            ),

          mut_n =
            nrow(
              expr_mut
            ),

          cor_df,
          coef_df
        )

        local_result[[paste(
          gene,
          cls,
          s_ppi,
          sep = "_"
        )]] <- result
      }
    }


    # ---------------------------------------------------------------------
    # Save per-query result
    # ---------------------------------------------------------------------

    if (length(local_result) > 0) {

      fwrite(
        rbindlist(
          local_result,
          use.names = TRUE,
          fill = TRUE
        ),
        file = outfile,
        sep = "\t",
        quote = FALSE,
        na = "NA"
      )
    }

    return(NULL)
  }


  # -----------------------------------------------------------------------
  # 5.6 Run query genes
  # -----------------------------------------------------------------------

  pbmclapply(
    genes_cohort,
    run_gene_analysis,
    mc.cores = num_cores
  )


  # -----------------------------------------------------------------------
  # 5.7 Combine results
  # -----------------------------------------------------------------------

  result_files <- list.files(
    tmp_dir,
    full.names = TRUE,
    pattern = "\\.tsv$"
  )

  if (length(result_files) == 0) {

    stop(
      "No result files generated for ",
      cancer_label
    )
  }

  final_result <- rbindlist(
    lapply(
      result_files,
      fread
    ),
    use.names = TRUE,
    fill = TRUE
  )

  final_result <- final_result[
    target != query
  ]


  # -----------------------------------------------------------------------
  # 5.8 Cancer-specific FDR
  # -----------------------------------------------------------------------

  final_result[
    ,
    mut_cor_p_adj :=
      p.adjust(
        mut_cor_p,
        method = "fdr"
      )
  ]

  final_result[
    ,
    wt_cor_p_adj :=
      p.adjust(
        wt_cor_p,
        method = "fdr"
      )
  ]

  int_p_col <- grep(
    "^x_group_Pr",
    colnames(final_result),
    value = TRUE
  )

  if (length(int_p_col) != 1) {

    stop(
      "Could not uniquely identify interaction P-value column for ",
      cancer_label
    )
  }

  final_result[
    ,
    int_p_adj :=
      p.adjust(
        get(
          int_p_col
        ),
        method = "fdr"
      )
  ]


  # -----------------------------------------------------------------------
  # 5.9 Save raw cancer-specific result
  # -----------------------------------------------------------------------

  raw_file <- file.path(
    raw_output_dir,
    paste0(
      "exprs.cor.lm.all_genes.partners.wt_only.",
      cancer,
      ".am_class.txt"
    )
  )

  fwrite(
    final_result,
    file = raw_file,
    sep = "\t",
    quote = FALSE,
    na = "NA"
  )

  message(
    "Written: ",
    raw_file
  )
}


# =========================================================================
# 6. Classify cancer-specific mutation-associated rewiring
# =========================================================================

read_subtype_result <- function(cancer) {

  fname <- file.path(
    raw_output_dir,
    paste0(
      "exprs.cor.lm.all_genes.partners.wt_only.",
      cancer,
      ".am_class.txt"
    )
  )

  x <- fread(
    fname,
    header = TRUE,
    data.table = FALSE
  )

  x %>%

    filter(
      mutation %in%
        c(
          "benign",
          "pathogenic"
        )
    ) %>%

    mutate(
      cancer_type =
        unname(
          cancer_labels[[cancer]]
        ),

      baseline_association = case_when(

        wt_cor_r >= r_cut &
          wt_cor_p_adj < fdr_cut ~
          "Positive",

        wt_cor_r <= -r_cut &
          wt_cor_p_adj < fdr_cut ~
          "Negative",

        TRUE ~
          "n.s."
      ),

      delta_cor_r =
        mut_cor_r -
        wt_cor_r,

      wt_slope =
        x_Estimate,

      mut_slope =
        x_Estimate +
        x_group_Estimate,

      rewiring_status = case_when(

        is.na(
          x_group_Pr.t.
        ) |
          x_group_Pr.t. >=
          interaction_p_cut ~
          "Not significant",

        baseline_association ==
          "Positive" &
          x_group_Estimate > 0 ~
          "Strengthening",

        baseline_association ==
          "Positive" &
          x_group_Estimate < 0 ~
          "Weakening",

        baseline_association ==
          "Negative" &
          x_group_Estimate < 0 ~
          "Strengthening",

        baseline_association ==
          "Negative" &
          x_group_Estimate > 0 ~
          "Weakening",

        TRUE ~
          "Not significant"
      )
    ) %>%

    distinct(
      query_refseq,
      target_refseq,
      mutation,
      .keep_all = TRUE
    )
}


subtype_pairs <- map_dfr(
  cancer_types,
  read_subtype_result
)


# =========================================================================
# 7. Individual-cancer rewiring-frequency summary
#
# This intentionally uses the RAW cancer-specific files before distinct(),
# matching the previous heatmap analysis.
# =========================================================================

summarize_cancer_rewiring <- function(cancer) {

  fname <- file.path(
    raw_output_dir,
    paste0(
      "exprs.cor.lm.all_genes.partners.wt_only.",
      cancer,
      ".am_class.txt"
    )
  )

  x <- fread(
    fname,
    header = TRUE,
    data.table = FALSE
  ) %>%

    filter(
      mutation %in%
        c(
          "benign",
          "pathogenic"
        )
    ) %>%

    mutate(
      baseline = case_when(

        wt_cor_r >= r_cut &
          wt_cor_p_adj < fdr_cut ~
          "pos",

        wt_cor_r <= -r_cut &
          wt_cor_p_adj < fdr_cut ~
          "neg",

        TRUE ~
          "none"
      ),

      significant_rewiring =
        !is.na(
          int_p_adj
        ) &
        !is.na(
          x_group_Pr.t.
        ) &
        x_group_Pr.t. <
          interaction_p_cut
    )


  bg <- x %>%
    count(
      baseline,
      mutation,
      name = "bg"
    )


  sig <- x %>%

    filter(
      significant_rewiring
    ) %>%

    count(
      baseline,
      mutation,
      name = "sig"
    )


  n_query <- x %>%

    filter(
      significant_rewiring
    ) %>%

    group_by(
      baseline,
      mutation
    ) %>%

    summarise(
      n_query =
        n_distinct(
          query
        ),

      .groups = "drop"
    )


  bg %>%

    left_join(
      sig,
      by = c(
        "baseline",
        "mutation"
      )
    ) %>%

    left_join(
      n_query,
      by = c(
        "baseline",
        "mutation"
      )
    ) %>%

    mutate(
      sig =
        coalesce(
          sig,
          0L
        ),

      n_query =
        coalesce(
          n_query,
          0L
        ),

      prop =
        sig / bg,

      prop_min5_queries =
        ifelse(
          n_query < 5,
          NA_real_,
          prop
        ),

      cancer =
        unname(
          cancer_labels[[cancer]]
        ),

      pair =
        paste(
          baseline,
          mutation,
          sep = "_"
        )
    )
}


individual_fraction <- map_dfr(
  cancer_types,
  summarize_cancer_rewiring
)


pair_levels <- c(
  "pos_benign",
  "neg_benign",
  "none_benign",
  "pos_pathogenic",
  "neg_pathogenic",
  "none_pathogenic"
)


individual_fraction <- individual_fraction %>%

  mutate(
    pair = factor(
      pair,
      levels = pair_levels
    ),

    cancer = factor(
      cancer,
      levels = unname(
        cancer_labels
      )
    )
  ) %>%

  arrange(
    cancer,
    pair
  ) %>%

  mutate(
    cancer =
      as.character(
        cancer
      ),

    pair =
      as.character(
        pair
      )
  )


cross_cancer_summary <- individual_fraction %>%

  group_by(
    pair
  ) %>%

  summarise(
    mean_prop =
      mean(
        prop_min5_queries,
        na.rm = TRUE
      ),

    median_prop =
      median(
        prop_min5_queries,
        na.rm = TRUE
      ),

    sd_prop =
      sd(
        prop_min5_queries,
        na.rm = TRUE
      ),

    n_cancer =
      sum(
        !is.na(
          prop_min5_queries
        )
      ),

    .groups = "drop"
  ) %>%

  mutate(
    pair = factor(
      pair,
      levels = pair_levels
    )
  ) %>%

  arrange(
    pair
  ) %>%

  mutate(
    pair =
      as.character(
        pair
      )
  )


# =========================================================================
# 8. WT/MUT pairwise-complete sample counts
#
# Only baseline Positive/Negative edges are counted here because these
# are the edges that will be passed to script 09.
# =========================================================================

count_pairwise_n <- function(cancer) {

  cancer_label <- unname(
    cancer_labels[[cancer]]
  )

  maf_cohort <- unname(
    maf_cohorts[[cancer]]
  )

  message(
    "Counting pairwise-complete samples: ",
    cancer_label
  )


  # -----------------------------------------------------------------------
  # 8.1 Cancer-specific MAF and AlphaMissense annotation
  # -----------------------------------------------------------------------

  maf_cancer <- maf %>%

    filter(
      COHORT ==
        maf_cohort
    )


  maf_ms <- maf_cancer %>%

    filter(
      Variant_Classification ==
        "Missense_Mutation"
    ) %>%

    mutate(
      aa_change =
        gsub(
          "^p\\.",
          "",
          Protein_Change
        )
    ) %>%

    as.data.table()


  maf_annot <- am[
    maf_ms,
    on = .(
      uniprot_id =
        `HGNC_UniProt_ID(supplied_by_UniProt)`,
      protein_variant =
        aa_change
    ),
    nomatch = 0
  ]


  # -----------------------------------------------------------------------
  # 8.2 Cancer-specific proteome
  # -----------------------------------------------------------------------

  smpl_cancer <- unique(
    sub(
      "_T",
      "",
      maf_cancer$Tumor_Sample_Barcode
    )
  )


  smpl_expr_cancer <- match_mutation_samples(
    expression_ids = expr_samples,
    mutation_ids = smpl_cancer
  )


  expr_cancer <- expr[
    ,
    smpl_expr_cancer,
    drop = FALSE
  ]


  # -----------------------------------------------------------------------
  # 8.3 Positive/Negative baseline pairs retained in this cancer
  # -----------------------------------------------------------------------

  pairs <- subtype_pairs %>%

    filter(
      cancer_type ==
        cancer_label,

      baseline_association %in%
        c(
          "Positive",
          "Negative"
        )
    )


  pair_groups <- pairs %>%

    group_split(
      query,
      mutation
    )


  out <- map_dfr(
    pair_groups,
    function(df) {

      gene <- df$query[1]
      cls <- df$mutation[1]


      # -------------------------------------------------------------------
      # WT samples
      # -------------------------------------------------------------------

      mut_any <- maf_cancer %>%

        filter(
          Hugo_Symbol ==
            gene
        ) %>%

        pull(
          Tumor_Sample_Barcode
        ) %>%

        sub(
          "_T",
          "",
          .
        ) %>%

        unique()


      mut_expr_all <- match_mutation_samples(
        expression_ids =
          smpl_expr_cancer,

        mutation_ids =
          mut_any
      )


      wt_samples <- setdiff(
        smpl_expr_cancer,
        mut_expr_all
      )


      # -------------------------------------------------------------------
      # MUT samples for this AlphaMissense class
      # -------------------------------------------------------------------

      mut_cls <- maf_annot[
        Hugo_Symbol == gene &
          am_class == cls
      ]


      mut_samples_raw <- unique(
        sub(
          "_T",
          "",
          mut_cls$Tumor_Sample_Barcode
        )
      )


      mut_samples <- match_mutation_samples(
        expression_ids =
          smpl_expr_cancer,

        mutation_ids =
          mut_samples_raw
      )


      # -------------------------------------------------------------------
      # Pairwise-complete counts
      # -------------------------------------------------------------------

      qref <- df$query_refseq[1]


      if (!qref %in% rownames(expr_cancer)) {

        return(
          df %>%

            transmute(
              query_refseq,
              target_refseq,
              mutation,

              wt_qt_n =
                NA_integer_,

              mut_qt_n =
                NA_integer_
            )
        )
      }


      targets <- df$target_refseq


      valid_target <-
        targets %in%
        rownames(expr_cancer)


      wt_qt_n <- rep(
        NA_integer_,
        length(
          targets
        )
      )


      mut_qt_n <- rep(
        NA_integer_,
        length(
          targets
        )
      )


      if (any(valid_target)) {

        target_valid <-
          targets[
            valid_target
          ]


        # WT
        q_wt <- !is.na(
          as.numeric(
            expr_cancer[
              qref,
              wt_samples,
              drop = TRUE
            ]
          )
        )


        t_wt <- !is.na(
          as.matrix(
            expr_cancer[
              target_valid,
              wt_samples,
              drop = FALSE
            ]
          )
        )


        wt_qt_n[
          valid_target
        ] <- rowSums(
          sweep(
            t_wt,
            2,
            q_wt,
            `&`
          )
        )


        # MUT
        q_mut <- !is.na(
          as.numeric(
            expr_cancer[
              qref,
              mut_samples,
              drop = TRUE
            ]
          )
        )


        t_mut <- !is.na(
          as.matrix(
            expr_cancer[
              target_valid,
              mut_samples,
              drop = FALSE
            ]
          )
        )


        mut_qt_n[
          valid_target
        ] <- rowSums(
          sweep(
            t_mut,
            2,
            q_mut,
            `&`
          )
        )
      }


      df %>%

        transmute(
          query_refseq,
          target_refseq,
          mutation,
          wt_qt_n,
          mut_qt_n
        )
    }
  )


  stopifnot(
    nrow(out) ==
      nrow(pairs)
  )


  pairs %>%

    left_join(
      out,
      by = c(
        "query_refseq",
        "target_refseq",
        "mutation"
      )
    )
}


subtype_pairs_n <- map_dfr(
  cancer_types,
  count_pairwise_n
)


# =========================================================================
# 9. Prepare edge-level input for script 09
# =========================================================================

individual_edges <- subtype_pairs_n %>%

  transmute(
    cancer =
      cancer_type,

    query,
    query_refseq,

    target,
    target_refseq,

    mutation_class =
      recode(
        mutation,
        benign =
          "Benign",
        pathogenic =
          "Pathogenic"
      ),

    baseline_association,

    wt_samples_n =
      wt_n,

    mut_samples_n =
      mut_n,

    wt_pairwise_n =
      wt_qt_n,

    mut_pairwise_n =
      mut_qt_n,

    wt_cor_r,
    wt_cor_p_adj,

    mut_cor_r,
    mut_cor_p_adj,

    delta_cor_r,

    wt_slope,
    mut_slope,

    interaction_estimate =
      x_group_Estimate,

    interaction_se =
      x_group_StdError,

    interaction_p =
      x_group_Pr.t.,

    interaction_p_adj =
      int_p_adj,

    rewiring_status
  )


# -------------------------------------------------------------------------
# 10. Save
# -------------------------------------------------------------------------

fwrite(
  individual_edges,
  file = file.path(
    output_dir,
    "individual_cancer_mutation_rewiring_edges.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


fwrite(
  individual_fraction,
  file = file.path(
    output_dir,
    "individual_cancer_rewiring_fraction.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


fwrite(
  cross_cancer_summary,
  file = file.path(
    output_dir,
    "individual_cancer_rewiring_cross_cancer_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# -------------------------------------------------------------------------
# 11. Reproducibility summary
# -------------------------------------------------------------------------

message(
  "Positive/Negative cancer-specific edges for script 09: ",
  nrow(
    individual_edges
  )
)


print(
  cross_cancer_summary,
  n = Inf,
  width = Inf
)


message(
  "Individual-cancer mutation-associated rewiring analysis completed."
)
