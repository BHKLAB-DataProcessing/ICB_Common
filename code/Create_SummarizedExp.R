library(data.table)
library(Biobase)
library(SummarizedExperiment)
library(GenomicRanges)
library(biomaRt)

clin_cols <- c(
  "patient" , "sex" , "age" , "primary" , "histo" , "stage" , 
  "response.other.info" , "recist" , "response" , "drug_type" , 
  "dna" , "rna" , "t.pfs" , "pfs" , "t.os" , "os"
)

added_cols <- c(
  "TMB_raw" , "nsTMB_raw" , "indel_TMB_raw" , "indel_nsTMB_raw" , 
  "TMB_perMb" , "nsTMB_perMb" , "indel_TMB_perMb" , "indel_nsTMB_perMb" ,
  "CIN" , "CNA_tot" , "AMP" , "DEL" 
)

renamed_cols <- list(
  drug_type = "treatment",
  primary = "cancer_type",
  t.pfs = "survival_time",
  pfs = "event_occured",
  patient = "unique_patient_ID"
)

format_se <- function(assay, coldata, assay_type, convert_gene_name=FALSE){
  # colnames(assay) <- str_replace_all(colnames(assay), '[-\\.]', '_')
  # rownames(coldata) <- str_replace_all(rownames(coldata), '[-\\.]', '_')
  # coldata$patient <- rownames(coldata)
  # assay <- expr
  # coldata <- clin
  # assay_type <- 'expr'
  
  
  
  for(renamed_col in names(renamed_cols)){
    colnames(coldata)[colnames(coldata) == renamed_col] <- renamed_cols[[renamed_col]]
  }
  
  coldata$survival_unit <- "month"
  coldata$survival_type <- "PFS"
  
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
    
    if(convert_gene_name){
      # replace assay gene names with gene id
      assay <- assay[rownames(assay) %in% features_gene$gene_name, ]
      gene_ids <- unlist(lapply(rownames(assay), function(assay_row){
        vals <- rownames(features_gene[features_gene$gene_name == assay_row, ])
        if(length(vals) > 1){
          return(vals[1])
        }else{
          return(vals)
        }
      }))
      rownames(assay) <- gene_ids 
    }
    
    # build the GRanges object ussed as rowRanges (rowData)
    assay_genes <- as.data.table(features_gene[rownames(features_gene) %in% gene_ids, ])
    assay_genes <- assay_genes[order(rownames(assay_genes)), ]
    assay_genes <- assay_genes[, gene_id_no_ver := gsub("\\..*$", "", gene_id)]
    assay_genes[
      is.na(start),
      c("start", "end", "length", "strand") := list(-1, -1, 0, "*")
    ]
    row_ranges <- makeGRangesFromDataFrame(
      assay_genes,
      keep.extra.columns=TRUE  # retain metadata
    )
    names(row_ranges) <- row_ranges$rownames
    
    assay <- assay[order(rownames(assay)), ]
    assay_list <- list()
    assay_list[[assay_type]] <- assay
    return(SummarizedExperiment(assays=assay_list, colData=coldata, rowRanges=row_ranges))
    
  }
}


Create_CNA_SummarizedExperiment = function( case, clin, cna, feat_snv , feat_cna , feat_cin , snv_bool , cna_bool ){
  # case = read.csv( case_file , sep=";" , stringsAsFactors=FALSE )
  # cna = read.csv( cna_file , sep=";" , stringsAsFactors=FALSE )
  # clin = read.csv( clin_file , sep=";" , stringsAsFactors=FALSE )
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
  rownames(feat_cin) = feat_cin$patient
  clin[ rownames(feat_cin) , "CIN" ] = feat_cin[ , "CIN" ]
  
  rownames(feat_cna) = feat_cna$patient
  clin[ rownames(feat_cna) , colnames(feat_cna)[-1] ] = feat_cna[ , c( "CNA_tot" , "AMP" , "DEL" ) ]
  clin = clin[ colnames(cna) , ]
  
  return(format_se(assay=cna, coldata=clin, assay_type='cna'))
}

Create_EXP_SummarizedExperiment = function( study, case , clin, expr, feat_snv, feat_cna, feat_cin, snv_bool, cna_bool ){
  # case = read.csv( case_file , sep=";" , stringsAsFactors=FALSE )
  # expr = read.csv( expr_file , sep=";" , stringsAsFactors=FALSE )
  # clin = read.csv( clin_file , sep=";" , stringsAsFactors=FALSE )
  
  study_with_gene_id <- c("Miao.1")
  
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
  
  return(format_se(assay=expr, coldata=clin, assay_type='expr', convert_gene_name=study %in% study_with_gene_id))
}

Create_SNV_SummarizedExperiment = function( case, clin, snv, feat_snv , feat_cna , feat_cin , cna_bool, snv_bool ){
  # case = read.csv( case_file , sep=";" , stringsAsFactors=FALSE )
  # snv = read.csv( snv_file , sep=";" , stringsAsFactors=FALSE )
  # clin = read.csv( clin_file , sep=";" , stringsAsFactors=FALSE )
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
    rownames(feat_cin) = feat_cin$patient
    clin[ rownames(feat_cin) , "CIN" ] = feat_cin[ , "CIN" ]
    rownames(feat_cna) = feat_cna$patient
    clin[ rownames(feat_cna) , colnames(feat_cna)[-1] ] = feat_cna[ , c( "CNA_tot" , "AMP" , "DEL" ) ]
  }
  
  clin = clin[ patient , ]

  return(format_se(assay=mat_snv, coldata=clin, assay_type='snv'))
}


Create_SummarizedExperiments = function( input_dir, study , expr_bool , snv_bool , cna_bool, cin_bool , coverage , indel_bool ){
  # study= data$study
  # expr_bool= data$expr_boo
  # snv_bool= data$snv_boo
  # cna_bool= data$cna_boo
  # cin_bool= data$cin_bool
  # coverage= data$coverage
  # indel_bool= data$indel_bool

  # Path to processed data 
	case_file = paste( input_dir , "cased_sequenced.csv" , sep="/" )
	clin_file = paste( input_dir , "CLIN.csv" , sep="/" )
	expr_file = paste( input_dir , "EXPR.csv" , sep="/" )
	snv_file = paste( input_dir , "SNV.csv" , sep="/" )
	cna_file = paste( input_dir , "CNA_gene.csv" , sep="/" )

	feat_snv <- NA
	feat_cna <- NA
	feat_cin <- NA
  se_list <- list()
	
  if( cna_bool ){ 	
	  feat_cin <- Get_CIN( case_file=case_file, cna_file=cna_file, cin_bool=cin_bool ) 
		feat_cna <- Get_CNA_feature( case_file=case_file, cna_file=cna_file, coverage=coverage, cna_bool=cna_bool ) 
	}

	if( snv_bool ){ 
		feat_snv <- Get_SNV_feature( case=case_file, file=snv_file, coverage=coverage , indel_bool=indel_bool , mutsig_bool=mutsig_bool ) 
		rownames(feat_snv) <- str_replace_all(rownames(feat_snv), '[-\\.]', '_')
	}
  
  case <- read.csv( case_file , sep=";" , stringsAsFactors=FALSE )
  case$patient <- str_replace_all(case$patient, '[-\\.]', '_')
  clin <- read.csv( clin_file , sep=";" , stringsAsFactors=FALSE )
  clin$patient <- str_replace_all(clin$patient, '[-\\.]', '_')

	if( cna_bool ){
	  cna = read.csv( cna_file , sep=";" , stringsAsFactors=FALSE )
	  colnames(cna) <- str_replace_all(colnames(cna), '[-\\.]', '_')
	  se_list[["cna"]] <- Create_CNA_SummarizedExperiment( case=case , clin=clin , cna=cna , 
									feat_snv=feat_snv , feat_cna=feat_cna , feat_cin=feat_cin , cna_bool=cna_bool , snv_bool=snv_bool )
	}

	if( expr_bool ){
	  expr = read.csv( expr_file , sep=";" , stringsAsFactors=FALSE )
	  colnames(expr) <- str_replace_all(colnames(expr), '[-\\.]', '_')
	  se_list[["expr"]] <- Create_EXP_SummarizedExperiment(study=study, case=case, clin=clin, expr=expr, 
									feat_snv=feat_snv , feat_cna=feat_cna, feat_cin=feat_cin , cna_bool=cna_bool , snv_bool=snv_bool )
	}

	if( snv_bool ){
	  snv = read.csv( snv_file , sep=";" , stringsAsFactors=FALSE )
	  snv$Sample <- str_replace_all(snv$Sample, '[-\\.]', '_')
	  se_list[["snv"]] <- Create_SNV_SummarizedExperiment( case=case , clin=clin , snv=snv , 
									feat_snv=feat_snv , feat_cna=feat_cna , feat_cin=feat_cin , cna_bool=cna_bool , snv_bool=snv_bool )
	}
  return(se_list)
}
