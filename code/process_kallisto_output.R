library(data.table)
library(tximport)

# Creates EXPR_gene_tpm.tsv, EXPR_gene_counts.tsv, EXPR_tx_tpm.tsv, EXPR_tx_counts.tsv 
process_kallisto_output <- function(work_dir, kallisto_zip, tx2gene){
  dir.create(file.path(work_dir, 'rnaseq'))
  unzip(file.path(work_dir, kallisto_zip), exdir=file.path(work_dir, 'rnaseq'))
  
  samples <- list.dirs(file.path(work_dir, 'rnaseq'), full.names=FALSE, recursive = FALSE)
  files <- file.path(work_dir, 'rnaseq', samples, "abundance.h5")
  names(files) <- samples
  
  expr_tx <- tximport(files, type = "kallisto", txOut = TRUE)
  expr_gene <- tximport(files, type = "kallisto", tx2gene = tx2gene)
  
  expr_list <- list()
  expr_list[['expr_gene_tpm']] <- expr_gene$abundance
  expr_list[['expr_gene_counts']] <- expr_gene$counts
  expr_list[['expr_isoform_tpm']] <- expr_tx$abundance
  expr_list[['expr_isoform_counts']] <- expr_tx$counts
  
  saveRDS(expr_list, file.path(work_dir, 'expr_list.rds'))
  
  unlink(file.path(work_dir, 'rnaseq'), recursive = TRUE)
}
