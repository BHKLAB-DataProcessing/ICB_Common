
format_transcript2gene <- function(tx2gene, expr){
  expr <- expr[rownames(expr) %in% tx2gene$transcripts, ]
  expr <- expr[order(rownames(expr)), ]
  tx2gene <- tx2gene[order(tx2gene$transcripts), ]
  
  expr_gene <- data.frame(matrix(ncol = length(colnames(expr)), nrow = length(unique(tx2gene$genes))))
  colnames(expr_gene) <- colnames(expr)
  rownames(expr_gene) <- unique(tx2gene$genes)
  
  for(gene in rownames(expr_gene)){
    transcripts <- tx2gene[tx2gene$genes == gene, ]$transcripts
    if(length(transcripts) > 1){
      df <- expr[transcripts, ]
      expr_gene[rownames(expr_gene) == gene, ] <- unlist(lapply(colnames(df), function(col){
        return(sum(df[, col]))
      }))
    }else{
      expr_gene[rownames(expr_gene) == gene, ] <- expr[transcripts, ]
    }
  }
  return(expr_gene)
}