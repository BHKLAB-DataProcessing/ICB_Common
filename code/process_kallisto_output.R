library(data.table)
library(tximport)

# Creates EXPR_gene_tpm.tsv, EXPR_gene_counts.tsv, EXPR_tx_tpm.tsv, EXPR_tx_counts.tsv 
process_kallisto_output <- function(work_dir, tx2gene){
  samples <- list.dirs(file.path(work_dir, 'rnaseq'), full.names=FALSE, recursive = FALSE)
  files <- file.path(work_dir, 'rnaseq', samples, "abundance.h5")
  names(files) <- samples
  
  expr_tx <- tximport(files, type = "kallisto", txOut = TRUE)
  expr_gene <- tximport(files, type = "kallisto", tx2gene = tx2gene)
  
  expr_list <- list()
  expr_list[['expr_gene_tpm']] <- log2(expr_gene$abundance + 0.001)
  expr_list[['expr_gene_counts']] <- log2(expr_gene$counts + 1)
  expr_list[['expr_isoform_tpm']] <- log2(expr_tx$abundance + 0.001)
  expr_list[['expr_isoform_counts']] <- log2(expr_tx$counts + 1)
  
  saveRDS(expr_list, file.path(work_dir, 'expr_list.rds'))
  
  unlink(file.path(work_dir, 'rnaseq'), recursive = TRUE)
}
