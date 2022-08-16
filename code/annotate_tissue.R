library(data.table)
library(stringr)

clin_cols <- c(
  "patient" , "sex" , "age" , "primary" , "histo" , "tissueid", "stage" , 
  "response.other.info" , "recist" , "response" , "drug_type" , 
  "dna" , "rna" , "t.pfs" , "pfs" , "t.os" , "os", 
  "survival_unit", "survival_type"
)

annotate_tissue <- function(clin, study, annotation_tissue, check_histo=FALSE){
  
  clin$survival_unit <- "month"
  clin$survival_type <- "PFS"
  
  study_annotation <- annotation_tissue[annotation_tissue$study == study, ]
  
  tissue_ids <- c()
  if(check_histo){
    tissue_ids <- unlist(lapply(clin$patient, function(patient){
      patient_cancer_type <- clin[clin$patient == patient, c('primary', 'histo')]
      cancer <- patient_cancer_type$histo
      if(!is.na(cancer) & length(study_annotation[study_annotation$cancer_type == cancer, c('tissueid')]) > 0){
        return(
          study_annotation[study_annotation$cancer_type == cancer, c('tissueid')]
        )
      }else{
        return(
          study_annotation[study_annotation$cancer_type == patient_cancer_type$primary, c('tissueid')]
        )
      }
    }))
  }else{
    tissue_ids <- unlist(lapply(clin$primary, function(cancer){
      return(
        study_annotation[study_annotation$cancer_type == cancer, c('tissueid')]
      )
    }))
  }

  clin$tissueid <- tissue_ids
  
  clin <- clin[, c(clin_cols, colnames(clin)[!colnames(clin) %in% clin_cols])]
  clin[, c('t.pfs', 't.os')] <- sapply(clin[, c('t.pfs', 't.os')], as.numeric)
  
  return(clin)
}

# mapping <- read.csv(file.path(annot_dir, 'curationTissue_mapping - cancer_types.tsv'), sep='\t')
# mapping <- mapping[, c('study', 'cancer_type', 'unique_tissueid')]
# annotation_tissue <- data.frame(matrix(nrow=0, ncol=3))
# colnames(annotation_tissue) <- c('study', 'cancer_type', 'unique_tissueid')
# 
# for(study in mapping$study){
#   study_annot <- mapping[mapping$study == study, ]
#   if(str_detect(study_annot$cancer_type, ';')){
#     cancer_types <- unlist(str_split(study_annot$cancer_type, ';'))
#     tissue_ids <- unlist(str_split(study_annot$unique_tissueid, ';'))
#     study_annot <- data.frame(
#       study = study,
#       cancer_type = cancer_types,
#       unique_tissueid = tissue_ids
#     )
#   }
#   study_annot$unique_tissueid[study_annot$unique_tissueid == '<blank>' | study_annot$unique_tissueid == "NA"] <- ''
#   annotation_tissue <- rbind(annotation_tissue, study_annot)
# }
# 
# write.table(annotation_tissue, file=file.path(annot_dir, 'curation_tissue.csv'), sep=',', row.names = FALSE, col.names = TRUE)