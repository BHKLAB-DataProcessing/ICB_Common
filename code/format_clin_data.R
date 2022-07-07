format_clin_data <- function(clin_original, id_colname, selected_cols, clin){
  rownames(clin) <- clin$patient
  appended <- clin_original[clin_original[[id_colname]] %in% rownames(clin), ]
  rownames(appended) <- appended[[id_colname]]
  appended <- appended[, colnames(appended)[!colnames(appended) %in% selected_cols], drop=FALSE]
  appended <- appended[order(rownames(appended)), , drop=FALSE]
  clin <- clin[order(rownames(clin)), ]
  clin <- cbind(clin, appended)
  return(clin)
}