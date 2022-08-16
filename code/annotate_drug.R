library(data.table)
library(stringr)

annotate_drug <- function(study, drug_names, annotation_drug){
  annotation_drug <- annotation_drug[annotation_drug$study == study, ]
  drug_ids <- unlist(lapply(drug_names, function(drug){
    if(!is.na(drug)){
      drugid <- annotation_drug[annotation_drug$dataset.drug.name == drug, c('treatmentid')]
      if(!is.na(drugid)){
        return(drugid)
      }else{
        return('')
      }
    }else{
      return('')
    }
  }))
  return(drug_ids)
}

# mapping <- read.csv(file.path(annot_dir, 'curationDrug_mapping - Sheet1.tsv'), sep='\t')
# mapping <- mapping[, c('study', 'dataset.drug.name', 'unique_drugid')]
# mapping <- mapping[mapping$dataset.drug.name != '', ]
# mapping$study <- str_replace(mapping$study, 'ICB_', '')
# 
# annotation_drug <- data.frame(matrix(nrow=0, ncol=3))
# colnames(annotation_drug) <- c('study', 'dataset.drug.name', 'unique_drugid')
# 
# for(study in mapping$study){
#   study_annot <- mapping[mapping$study == study, ]
#   if(str_detect(study_annot$dataset.drug.name, ';')){
#     drugs <- unlist(str_split(study_annot$dataset.drug.name, ';'))
#     drug_ids <- unlist(str_split(study_annot$unique_drugid, ';'))
#     study_annot <- data.frame(
#       study = study,
#       dataset.drug.name = drugs,
#       unique_drugid = drug_ids
#     )
#   }
#   annotation_drug <- rbind(annotation_drug, study_annot)
# }
# 
# write.table(annotation_drug, file=file.path(annot_dir, 'curation_drug.csv'), sep=',', row.names = FALSE, col.names = TRUE)