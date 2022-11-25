library(data.table)
library(stringr)

source_location <- "https://raw.githubusercontent.com/BHKLAB-DataProcessing/ICB_Common/main"
source(paste(source_location, "code", "format_se.R", sep = "/"))

args <- commandArgs(trailingOnly = TRUE)
work_dir <- args[1]
cna_bool <- args[2]
annotation_file <- args[3]

case = read.csv( file.path(work_dir, 'cased_sequenced.csv') , sep=";" , stringsAsFactors=FALSE )
snv = read.csv( file.path(work_dir, 'SNV.csv') , sep=";" , stringsAsFactors=FALSE )
clin = read.csv( file.path(work_dir, 'CLIN.csv') , sep=";" , stringsAsFactors=FALSE )
feat_nsv <- readRDS(file.path(work_dir, 'feat_snv.rds'))

snv$Sample <- str_replace_all(snv$Sample, '[-\\.]', '_')

snv = snv[ snv$Effect %in% c("In_Frame_Del" , "In_Frame_Ins" , "Start_Codon_Ins" , "Frame_Shift_Del" ,
                             "Frame_Shift_Ins" , "Missense_Mutation" , "Nonsense_Mutation" , "Nonstop_Mutation" ,
                             "Splice_Site" , "Stop_Codon_Del" , "De_novo_Start_OutOfFrame" , "Start_Codon_SNP") ,
]

genes = sort( unique( snv$Gene ) )
patient = case[ case$snv %in% 1 , ]$patient

mat_snv = matrix( NA , nrow=length(genes) , ncol=length(patient) )
colnames(mat_snv) = patient
rownames(mat_snv) = genes

sample = patient

for(i in 1:length(sample)){
  s = snv[ snv$Sample %in% sample[i], ]
  if( nrow(s) ){
    for(j in 1:nrow(s)){
      if( !is.na( s$Gene[j]) & nchar(s$Gene[j]) > 0){ 
        if( !is.na( mat_snv[ s$Gene[j] , sample[i] ] ) ){
          mat_snv[ s$Gene[j] , sample[i] ]  = paste( mat_snv[ s$Gene[j] , sample[i] ] , paste( s$Ref[j] , s$Alt[j] , sep=">" ) , sep=";" )
        } else{
          mat_snv[ s$Gene[j] , sample[i] ]  = paste( s$Ref[j] , s$Alt[j] , sep=">" )					
        }
      }
    }
  }
}

rownames(clin) = clin$patient
added_df <- as.data.frame(matrix(NA, nrow = nrow(clin), ncol = length(added_cols)))
colnames(added_df) <- added_cols
clin <- data.frame(cbind(
  clin[, clin_cols],
  added_df,
  clin[, !colnames(clin) %in% clin_cols]
))

clin[ rownames(feat_snv) , colnames(feat_snv) ] = feat_snv

if(cna_bool){ 
  feat_cin <- readRDS(file.path(work_dir, 'feat_cin.rds'))
  feat_cna <- readRDS(file.path(work_dir, 'feat_cna.rds'))
  rownames(feat_cin) = feat_cin$patient
  clin[ rownames(feat_cin) , "CIN" ] = feat_cin[ , "CIN" ]
  rownames(feat_cna) = feat_cna$patient
  clin[ rownames(feat_cna) , colnames(feat_cna)[-1] ] = feat_cna[ , c( "CNA_tot" , "AMP" , "DEL" ) ]
}

clin = clin[ patient , ]

snv_se <- format_se(
  assay=mat_snv, 
  coldata=clin, 
  assay_type='snv', 
  gene_annotation_file=annotation_file
)

saveRDS(snv_se, file.path(work_dir, 'SNV.rds'))
