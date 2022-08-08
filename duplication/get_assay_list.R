library(MultiAssayExperiment)

input_dir <- "~/Documents/ICBCuration/results_from_source"
output_dir <- "~/Documents/ICBCuration/duplication/results"
assay_name <- 'expr'

## EXPR
studies <- c(
  "Braun", "Gide", "Hugo", "Hwang", "Jerby_Arnon", "Jung", "Kim", "Liu", "Mariathasan", "Miao1", 
  "Nathanson", "Padron", "Puch", "Riaz", "Roh", "Shiuan", "Snyder", "Van_Allen", "VanDenEnde"
)
studies_with_multiple_expr <- c(
  "Gide", "Hugo", "Jung", "Kim", "Riaz"
)

# SNV
studies <- c(
  "Braun", "Jung", "Liu", "Miao1", "Nathanson", 
  "Riaz", "Roh", "Van_Allen"
)

esets <- lapply(studies, function(study){
  print(study)
  dat <- readRDS(file.path(input_dir, paste0('ICB_', study, '.rds')))
  exp <- experiments(dat)
  
  if(study %in% studies_with_multiple_expr){
    assay_name <- 'expr_gene_tpm'
  }
  
  dat_exp_bhk <- assay(exp[[assay_name]])
  dat_exp_bhk <- dat_exp_bhk[order(rownames(dat_exp_bhk)), ]
  
  dat_clinic_bhk <- data.frame(colData(dat))
  dat_clinic_bhk  <- dat_clinic_bhk [order(rownames(dat_clinic_bhk )), ]
  
  if(assay_name == 'snv'){
    dat_gene_bhk <- data.frame(rowData(exp[[assay_name]]))
    dat_gene_bhk$mod_gene_name <- sapply(1:nrow(dat_exp_bhk), function(j){
      dat_gene_bhk[dat_gene_bhk$gene_name == rownames(dat_exp_bhk)[j], "gene_name"]
    })
    rownames(dat_exp_bhk) <- dat_gene_bhk$mod_gene_name
  }

  dat_exp_bhk <- dat_exp_bhk[order(rownames(dat_exp_bhk)), ]
  dat_exp_bhk <- dat_exp_bhk[, order(colnames(dat_exp_bhk))]
  
  int <- intersect(colnames(dat_exp_bhk), rownames(dat_clinic_bhk))
  dat_clinic_bhk <- dat_clinic_bhk[rownames(dat_clinic_bhk) %in% int, ]
  
  
  eset <- ExpressionSet(assayData = as.matrix(dat_exp_bhk),
                        phenoData = AnnotatedDataFrame(dat_clinic_bhk))
  return(eset)
})

names(esets) <- studies
saveRDS(esets, file.path(output_dir, paste0(assay_name, '_list.rds')))


