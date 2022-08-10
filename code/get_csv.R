library(stringr)
library(data.table)

input_dir <- '~/Documents/ICBCuration/results_from_source'
output_dir <- '~/Documents/ICBCuration/results_csv' 

assay_names <- c(
  'expr',
  'expr_gene_tpm',
  'expr_gene_counts',
  'expr_isoform_tpm',
  'expr_isoform_counts',
  'snv',
  'cna'
)

objects <- list.files(input_dir)

for(file in objects){
  obj_name <- str_extract(file, '.*(?=.rds)')
  obj <- readRDS(file.path(input_dir, file))
  dir.create(file.path(output_dir, obj_name))
  coldata <- data.frame(colData(obj))
  write.table(
    coldata, 
    file=file.path(output_dir, obj_name, paste0(obj_name, '_coldata', '.tsv')), 
    sep='\t',
    col.names=TRUE, 
    row.names=TRUE
  )
  
  for(assay_name in assay_names){
    experiment <- experiments(obj)[[assay_name]]
    if(!is.null(experiment)){
      assay <- data.frame(assay(experiment))
      assay_genes <- NA
      if(str_detect(assay_name, 'snv')){
        assay_genes <- data.frame(elementMetadata(experiment))
      }else{
        assay_genes <- data.frame(rowRanges(experiment))
      }
      write.table(
        assay, 
        file=file.path(output_dir, obj_name, paste0(obj_name, '_', assay_name, '.tsv')),  
        sep='\t',
        col.names=TRUE, 
        row.names=TRUE
      )
      write.table(
        assay_genes, 
        file=file.path(output_dir, obj_name, paste0(obj_name, '_', assay_name, '_genes.tsv')),  
        sep='\t',
        col.names=TRUE, 
        row.names=TRUE
      )
    }
  }

  zip(
    zipfile=file.path(output_dir, paste0(obj_name, '.zip')), 
    files=list.files(file.path(output_dir, obj_name), full.names = TRUE),
    flags = '-r9Xj'
  )
  
  unlink(file.path(output_dir, obj_name), recursive=TRUE)
}
