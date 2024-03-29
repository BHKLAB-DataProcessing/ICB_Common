library(data.table)
library(SummarizedExperiment)
library(stringr)

source_location <- "https://raw.githubusercontent.com/BHKLAB-DataProcessing/ICB_Common/main"
source(paste(source_location, "code", "format_se.R", sep = "/"))

args <- commandArgs(trailingOnly = TRUE)
work_dir <- args[1]
snv_bool <- args[2]
annotation_file <- args[3]

cna <- read.csv( file.path(work_dir, 'CNA_gene.csv') , sep=";" , stringsAsFactors=FALSE )
clin <- readRDS( file.path(work_dir, 'CLIN.rds'))
cases <- readRDS( file.path(work_dir, 'cased_sequenced.rds'))

feat_cin <- readRDS(file.path(work_dir, 'feat_cin.rds')) 
feat_cna <- readRDS(file.path(work_dir, 'feat_cna.rds'))  

rownames(clin) = clin$patient
added_df <- as.data.frame(matrix(NA, nrow = nrow(clin), ncol = length(added_cols)))
colnames(added_df) <- added_cols
clin <- data.frame(cbind(
  clin[, clin_cols],
  added_df,
  clin[, !colnames(clin) %in% clin_cols]
))

if(snv_bool){ 
  feat_snv <- readRDS(file.path(work_dir, 'feat_snv.rds'))
  clin[ rownames(feat_snv) , colnames(feat_snv) ] = feat_snv
}

rownames(feat_cin) = feat_cin$patient
clin[ rownames(feat_cin) , "CIN" ] = feat_cin[ , "CIN" ]

rownames(feat_cna) = feat_cna$patient
clin[ rownames(feat_cna) , colnames(feat_cna)[-1] ] = feat_cna[ , c( "CNA_tot" , "AMP" , "DEL" ) ]
clin = clin[ colnames(cna) , ]

cna_se <- format_se(
  assay=cna, 
  coldata=clin, 
  assay_type='cna', 
  gene_annotation_file=annotation_file
)

saveRDS(cna_se, file.path(work_dir, 'CNA.rds'))