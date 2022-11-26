library(data.table)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
work_dir <- args[1]
study <- args[2]
snv_bool <- args[3]
cna_bool <- args[4]
expr_with_counts_isoforms <- args[5]
annotation_file <- args[6]

Create_EXP_SummarizedExperiment = function( study, assay_name, case , clin, expr, feat_snv, feat_cna, feat_cin, is_isoform=FALSE ){
  
  study_with_gene_id <- c("Braun", "Gide", "Hugo", "INSPIRE", "Jung", "Kim", "Miao.1", "Nathanson", "Riaz", "Snyder", "Van_Allen", "VanDenEnde")
  
  rownames(clin) = clin$patient
  added_df <- as.data.frame(matrix(NA, nrow = nrow(clin), ncol = length(added_cols)))
  colnames(added_df) <- added_cols
  clin <- data.frame(cbind(
    clin[, clin_cols],
    added_df,
    clin[, !colnames(clin) %in% clin_cols]
  ))
  
  if(snv_bool){ 
    clin[ rownames(feat_snv) , colnames(feat_snv) ] = feat_snv
  }
  
  if(cna_bool){ 
    rownames(feat_cin) = feat_cin$patient
    clin[ rownames(feat_cin) , "CIN" ] = feat_cin[ , "CIN" ]
    rownames(feat_cna) = feat_cna$patient
    clin[ rownames(feat_cna) , colnames(feat_cna)[-1] ] = feat_cna[ , c( "CNA_tot" , "AMP" , "DEL" ) ]
    
  }
  
  patient = intersect( colnames(expr) , rownames(clin) )
  clin = clin[ patient , ]
  expr = expr[ , patient ]
  
  expr_se <- format_se(
    assay=expr, 
    coldata=clin, 
    assay_type='expr', 
    convert_gene_name=!study %in% study_with_gene_id, 
    is_isoform=is_isoform, 
    gene_annotation_file=annotation_file
  )
  saveRDS(expr_se, file.path(work_dir, paste0(assay_name, '.rds')))
}

case <- readRDS(file.path(work_dir, 'cased_sequenced.rds'))
clin <- readRDS( file.path(work_dir, 'CLIN.rds'))
feat_snv <- NA
feat_cna <- NA
feat_cin <- NA

if(snv_bool){
  feat_snv <- readRDS(file.path(work_dir, 'feat_snv.rds'))
}

if(cna_bool){
  feat_cna <- readRDS(file.path(work_dir, 'feat_cna.rds'))
  feat_cin <- readRDS(file.path(work_dir, 'feat_cin.rds'))
}

if(expr_with_counts_isoforms){
  expr_files <- c('EXPR_gene_tpm.csv', 'EXPR_gene_counts.csv', 'EXPR_isoform_tpm.csv', 'EXPR_isoform_counts.csv')
  for(expr_file in expr_files){
    expr = read.csv( file.path(input_dir, expr_file) , sep=";" , stringsAsFactors=FALSE )
    colnames(expr) <- str_replace_all(colnames(expr), '[-\\.]', '_')
    
    Create_EXP_SummarizedExperiment(
      study=study,
      assay_name=str_replace(expr_file, '.csv', ''),
      case=case, 
      clin=clin, 
      expr=expr, 
      feat_snv=feat_snv, 
      feat_cna=feat_cna, 
      feat_cin=feat_cin, 
      is_isoform=str_detect(expr_file, 'isoform')
    )
  }
}else{
  expr = read.csv( file.path(work_dir, 'EXPR.csv') , sep=";" , stringsAsFactors=FALSE )
  colnames(expr) <- str_replace_all(colnames(expr), '[-\\.]', '_')
  
  Create_EXP_SummarizedExperiment(
    study=study,
    assay_name='EXPR',
    case=case, 
    clin=clin, 
    expr=expr, 
    feat_snv=feat_snv, 
    feat_cna=feat_cna, 
    feat_cin=feat_cin
  )
}