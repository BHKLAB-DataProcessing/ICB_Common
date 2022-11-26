library(Biobase)
library(SummarizedExperiment)
library(MultiAssayExperiment)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
work_dir <- args[1]
output_dir <- args[2]
study <- args[3]
assays <- args[4]

se_list <- list()
for(assay in unlist(str_split(assays, ':'))){
  se_list[[tolower(assay)]] <- readRDS(file.path(work_dir, paste0(assay, '.rds')))
}

cols <-list()
for(assay_name in names(se_list)){
  cols[[assay_name]]<- data.frame(colData(se_list[[assay_name]]))
}

# Format and merge coldata
allcols <- lapply(cols, function(col){
  return(rownames(col))
})
allcols <- unique(unlist(allcols))

coldata <- NA
for(col in cols){
  if(!is.data.frame(coldata)){
    coldata <- col
  }
  missing <- allcols[!allcols %in% rownames(coldata)]
  filtered <- col[rownames(col) %in% missing, ]
  if(nrow(filtered) > 0){
    coldata <- rbind(coldata, filtered)
  }
}
coldata <- coldata[order(rownames(coldata)), ]

multiassay <- MultiAssayExperiment(experiments=se_list, colData=coldata)
saveRDS(multiassay, paste0(output_dir, study, 'rds'))