#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript run_coloc_stream.R <trait> <tissue>\n", call. = FALSE)
}
trait  <- args[1]
tissue <- args[2]

suppressPackageStartupMessages({
  library(coloc)
  library(dplyr)
  library(data.table)
})

base_dir <- paste0("/faststorage/project/farmgtex/zhudi/human_GTEx_GWAS_coloc/",
                   trait, "/", tissue)
dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)

eqtl_path      <- paste0("/faststorage/project/farmgtex/zhudi/human_GTEx/GTEx_v10_eqtl_summary/",
                         tissue, "/egene_SNP_pair_hg19_filtered_X23.txt")
gwas_path      <- Sys.glob(paste0("/faststorage/project/farmgtex/zhudi/human_GWAS_ldsc/",
                                  trait, ".tsv"))
condition_path <- paste0(base_dir, "/save_egene.txt")

summary_out_path <- paste0(base_dir, "/coloc_result_summary")
pip_out_path     <- paste0(base_dir, "/all_PIP0.8")
log_path         <- paste0(base_dir, "/coloc_log.txt")


eqtl <- fread(eqtl_path,      header = TRUE)
gwas <- fread(gwas_path,      header = F)

conditions <- fread(condition_path, header = FALSE)
gwas[[1]] <- as.character(gwas[[1]])
conditions[[2]] <- as.character(conditions[[2]])

summary_written <- FALSE
pip_written     <- FALSE

total_genes     <- nrow(conditions)
processed_count <- 0
success_count   <- 0
skip_count      <- 0
error_count     <- 0

for (i in seq_len(nrow(conditions))) {
  processed_count <- processed_count + 1
  
  if (processed_count %% 10 == 0 || processed_count == 1) {
    cat(sprintf("process: %d/%d (%.1f%%)\n", processed_count, total_genes, 
                100 * processed_count / total_genes))
  }

  gene  <- conditions[[1]][i]
  chr   <- conditions[[2]][i]
  start <- conditions[[3]][i]
  end   <- conditions[[4]][i]

  eqtl_use <- eqtl[eqtl[[1]] == gene, ]

  if (nrow(eqtl_use) == 0) {
    skip_count <- skip_count + 1
    next
  }

  eqtl_input <- data.table(
    rs_id        = eqtl_use[[4]],
    phenotype_id = eqtl_use[[1]],
    variant_id   = eqtl_use[[4]],
    af           = eqtl_use[[6]],
    pval_nominal = eqtl_use[[9]],
    slope        = eqtl_use[[10]],
    slope_se     = eqtl_use[[11]]
  )

  chr <- as.character(chr)
  
  gwas_use <- gwas[
    gwas[[2]] >= start &
      gwas[[2]] <= end &
      gwas[[1]] == chr,
  ]
  
  if (nrow(gwas_use) == 0) {
    skip_count <- skip_count + 1
    next
  }

  gwas_input <- data.table(
    rs_id        = gwas_use[[3]],
    variant_id   = gwas_use[[3]],
    beta         = gwas_use[[4]],
    sebeta       = gwas_use[[5]],
    pval_nominal = gwas_use[[6]]
  )

  input <- merge(
    eqtl_input, gwas_input,
    by = "rs_id",
    all = FALSE,
    suffixes = c("_eqtl", "_gwas")
  )

  if (nrow(input) == 0) {
    skip_count <- skip_count + 1
    next
  }

  input$pval_nominal_eqtl <- as.numeric(input$pval_nominal_eqtl)
  input$pval_nominal_gwas <- as.numeric(input$pval_nominal_gwas)
  input$af                <- as.numeric(input$af)

  keep_idx <- !is.na(input$pval_nominal_eqtl) &
              !is.na(input$pval_nominal_gwas) &
              !is.na(input$af)
  input <- input[keep_idx, ]
  if (nrow(input) == 0) {
    skip_count <- skip_count + 1
    next
  }

  if (nrow(input) < 2) {
    skip_count <- skip_count + 1
    next
  }
  
  valid_maf <- input$af > 0 & input$af < 1
  if (sum(valid_maf) < 2) {
    skip_count <- skip_count + 1
    next
  }
  input <- input[valid_maf, ]
  
  valid_pval_eqtl <- input$pval_nominal_eqtl > 0 & input$pval_nominal_eqtl <= 1
  valid_pval_gwas <- input$pval_nominal_gwas > 0 & input$pval_nominal_gwas <= 1
  if (sum(valid_pval_eqtl & valid_pval_gwas) < 2) {
    skip_count <- skip_count + 1
    next
  }
  input <- input[valid_pval_eqtl & valid_pval_gwas, ]
  
  if (nrow(input) < 2) {
    skip_count <- skip_count + 1
    next
  }

if (any(duplicated(input$variant_id_gwas)) || any(duplicated(input$variant_id_eqtl))) {
  input <- input[!duplicated(input$variant_id_gwas), ]
}

if (nrow(input) < 2) {
  skip_count <- skip_count + 1
  next
}
  result <- tryCatch(
    {
      coloc.abf(
        dataset1 = list(
          pvalues = input$pval_nominal_gwas,
          snp     = input$variant_id_gwas,
          type    = "quant",
          N       = 50000
        ),
        dataset2 = list(
          pvalues = input$pval_nominal_eqtl,
          snp     = input$variant_id_eqtl,
          type    = "quant",
          N       = 120
        ),
        MAF = input$af
      )
    },
    error = function(e) {
      error_msg <- conditionMessage(e)
      cat(sprintf("worning: gene %s failed: %s\n", gene, error_msg))
      return(list(error = error_msg))
    }
  )

  if (!is.null(result$error) || inherits(result, "try-error")) {
    error_count <- error_count + 1
    next
  }
  if (is.null(result$results) || nrow(result$results) == 0) {
    skip_count <- skip_count + 1
    next
  }
  
  success_count <- success_count + 1

  top_result <- result$results %>%
    slice_max(SNP.PP.H4, n = 1, with_ties = FALSE)

  top_snp_id   <- top_result$snp
  top_snp_pph4 <- top_result$SNP.PP.H4

  summary_row <- data.table(
    trait        = trait,
    tissue       = tissue,
    gene         = gene,
    chr          = chr,
    start        = start,
    end          = end,
    top_snp      = top_snp_id,
    top_snp_pph4 = top_snp_pph4,
    summary5     = result$summary[5],
    summary6     = result$summary[6]
  )

  fwrite(
    summary_row,
    file       = summary_out_path,
    sep        = "\t",
    quote      = FALSE,
    append     = summary_written,
    col.names  = !summary_written
  )
  summary_written <- TRUE
  high_pip <- result$results %>%
    dplyr::filter(SNP.PP.H4 > 0.8)

  if (nrow(high_pip) > 0) {
    high_pip <- high_pip %>%
      mutate(
        trait  = trait,
        tissue = tissue,
        gene   = gene,
        chr    = chr,
        start  = start,
        end    = end
      )

    high_pip_dt <- as.data.table(high_pip)

    fwrite(
      high_pip_dt,
      file       = pip_out_path,
      sep        = "\t",
      quote      = FALSE,
      append     = pip_written,
      col.names  = !pip_written
    )
    pip_written <- TRUE
  }
}

cat("\n", summary_msg, "\n", sep = "")

writeLines(summary_msg, con = log_path)

