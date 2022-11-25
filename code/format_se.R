library(data.table)
library(SummarizedExperiment)
library(GenomicRanges)
library(stringr)

clin_cols <- c(
  "patient" , "sex" , "age" , "primary" , "histo" , "tissueid", "treatmentid", "stage" , 
  "response.other.info" , "recist" , "response" , "drug_type" , 
  "dna" , "rna" , "t.pfs" , "pfs" , "t.os" , "os", 
  "survival_unit", "survival_type"
)

added_cols <- c(
  "TMB_raw" , "nsTMB_raw" , "indel_TMB_raw" , "indel_nsTMB_raw" , 
  "TMB_perMb" , "nsTMB_perMb" , "indel_TMB_perMb" , "indel_nsTMB_perMb" ,
  "CIN" , "CNA_tot" , "AMP" , "DEL" 
)

renamed_cols <- list(
  drug_type = "treatment",
  primary = "cancer_type",
  t.pfs = "survival_time_pfs",
  pfs = "event_occurred_pfs",
  t.os = "survival_time_os",
  os = "event_occurred_os",
  patient = "patientid"
)

format_se <- function(assay, coldata, assay_type, convert_gene_name=TRUE, is_isoform=FALSE){
  # colnames(assay) <- str_replace_all(colnames(assay), '[-\\.]', '_')
  # rownames(coldata) <- str_replace_all(rownames(coldata), '[-\\.]', '_')
  # coldata$patient <- rownames(coldata)
  # assay <- expr
  # coldata <- clin
  # assay_type <- 'expr'
  
  for(renamed_col in names(renamed_cols)){
    colnames(coldata)[colnames(coldata) == renamed_col] <- renamed_cols[[renamed_col]]
  }
  
  # coldata$survival_unit <- "month"
  # coldata$survival_type <- "PFS"
  
  if(assay_type == "snv"){
    # Subset features_gene by gene_name and remove duplicates, if any.
    assay_genes <- features_gene[features_gene$gene_name %in% rownames(assay), ]
    duplicates <- lapply(unique(assay_genes$gene_name), function(gene){
      filtered <- assay_genes[assay_genes$gene_name == gene, ]
      if(nrow(filtered) > 1){
        return(rownames(filtered)[2:length(rownames(filtered))])
      }else{
        return(NA)
      }
    })
    duplicates <- unlist(duplicates[!is.na(duplicates)])
    assay_genes <- assay_genes[!rownames(assay_genes) %in% duplicates, ]
    rownames(assay_genes) <- assay_genes$gene_name
    
    # Add missing genes to assay_genes dataframe.
    missing <- rownames(assay)[!rownames(assay) %in% assay_genes$gene_name]
    if(length(missing) > 0){
      added <- data.frame(matrix(data=NA, ncol=ncol(assay_genes), nrow=(length(missing))))
      colnames(added) <- colnames(assay_genes)
      rownames(added) <- missing
      added$gene_name <- missing
      assay_genes <- rbind(assay_genes, added)
    }
    
    # fill in the range values
    assay_genes <- as.data.table(assay_genes)
    assay_genes[
      is.na(start),
      c("start", "end", "length", "strand") := list(-1, -1, 0, "*")
    ]
    
    assay_genes <- data.frame(assay_genes)
    rownames(assay_genes) <- assay_genes$gene_name
    assay_genes <- assay_genes[order(rownames(assay_genes)), ]
    assay <- assay[order(rownames(assay)), ]
    return(SummarizedExperiment(assays=list("snv"=assay), colData=coldata, rowData=assay_genes))
    
  }else{
    
    gene_ids <- c()
    features_df <- features_gene
    
    if(is_isoform){
      features_df <- features_transcript
    }
    
    if(convert_gene_name){
      # replace assay gene names with gene id
      assay <- assay[rownames(assay) %in% features_df$gene_name, ]
      gene_ids <- unlist(lapply(rownames(assay), function(assay_row){
        vals <- rownames(features_df[features_df$gene_name == assay_row, ])
        if(length(vals) > 1){
          return(vals[1])
        }else{
          return(vals)
        }
      }))
      rownames(assay) <- gene_ids 
    }else{
      assay <- assay[rownames(assay) %in% rownames(features_df), ]
      gene_ids <- rownames(assay)
    }
    
    assay <- assay[order(rownames(assay)), ]
    
    # build the GRanges object ussed as rowRanges (rowData)
    assay_genes <- as.data.table(features_df[rownames(features_df) %in% gene_ids, ])
    assay_genes <- assay_genes[order(rownames(assay_genes)), ]
    assay_genes <- assay_genes[, gene_id_no_ver := gsub("\\..*$", "", gene_id)]
    assay_genes[
      is.na(start),
      c("start", "end", "length", "strand") := list(-1, -1, 0, "*")
    ]
    assay_genes <- assay_genes[order(assay_genes$gene_id), ]
    row_ranges <- makeGRangesFromDataFrame(
      assay_genes,
      keep.extra.columns=TRUE  # retain metadata
    )
    
    names(row_ranges) <- row_ranges$rownames
    
    assay_list <- list()
    assay_list[[assay_type]] <- assay
    return(SummarizedExperiment(assays=assay_list, colData=coldata, rowRanges=row_ranges))
    
  }
}