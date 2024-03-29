library(Biobase)
library(SummarizedExperiment)
library(MultiAssayExperiment)
library(stringr)

source_location <- "https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main"
# source_location <- "~/Documents/GitHub/Pachyderm/PredictIO_ICB/ICB_Common"

get_MultiAssayExp <- function(study, input_dir, expr_with_counts_isoforms=FALSE){ 
  source(paste(source_location, "code", "CNA_Functions.R", sep = "/"))
  source(paste(source_location, "code", "SNV_Functions.R", sep = "/"))
  source(paste(source_location, "code", "Create_SummarizedExp.R", sep = "/"))
  data = read.csv( paste(source_location, "data", "DATASET_LOAD_INFO.csv", sep="/") , sep=";" , stringsAsFactors=FALSE )
  data <- data[data$study == study, ]
  
  se_list <- Create_SummarizedExperiments( 
    input_dir=input_dir,
    study= data$study, 
    expr_bool= data$expr_bool, 
    snv_bool= data$snv_bool, 
    cna_bool= data$cna_bool, 
    cin_bool= data$cin_bool, 
    coverage= data$coverage, 
    indel_bool= data$indel_bool,
    expr_with_counts_isoforms=expr_with_counts_isoforms
  )
  
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

  return(MultiAssayExperiment(experiments=se_list, colData=coldata))
}
