library(data.table)

process_abundance <- function(sample_dir, samples){
  transcripts <- c()
  for(sample in samples){
    abundance <- read.table(file.path(sample_dir, sample, 'abundance.tsv'), header=TRUE, sep='\t')
    transcripts <- c(transcripts, abundance$target_id)
    transcripts <- unique(transcripts)
  }
  
  expr <- data.frame(matrix(nrow=length(transcripts), ncol=length(samples)))
  rownames(expr) <- transcripts
  colnames(expr) <- samples
  expr <- expr[order(rownames(expr)), ]
  
  for(sample in samples){
    print(sample)
    abundance <- read.table(file.path(sample_dir, sample, 'abundance.tsv'), header=TRUE, sep='\t')
    rownames(abundance) <- abundance$target_id
    missing <- rownames(expr)[!rownames(expr) %in% rownames(abundance)]
    if(length(missing) > 0){
      added <- data.frame(matrix(nrow=length(missing), ncol=length(colnames(abundance))))
      rownames(added) <- missing
      colnames(added) <- colnames(abundance)
      abundance <- rbind(abundance, added)
    }
    abundance <- abundance[order(rownames(abundance)), ]
    expr[[sample]] <- abundance$tpm
  }
  return(expr)
}
