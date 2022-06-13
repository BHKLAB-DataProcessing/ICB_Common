library(data.table)
library(Biobase)
library(SummarizedExperiment)
library(biomaRt)

se_colnames <- c( 
  "patient" , "sex" , "age" , "primary" , "histo" , "stage" , 
  "response.other.info" , "recist" , "response" , "drug_type" , 
  "dna" , "rna" , "t.pfs" , "pfs" , "t.os" , "os" ,
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

format_se <- function(assay, coldata){
  colnames(assay) <- str_replace_all(colnames(assay), '[-\\.]', '_')
  rownames(coldata) <- str_replace_all(rownames(coldata), '[-\\.]', '_')
  coldata$patient <- rownames(coldata)
  
  colnames(coldata) <- unlist(lapply(se_colnames, function(col){
    if(is.null(renamed_cols[[col]])){
      return(col)
    }
    return(renamed_cols[[col]])
  }))
  
  coldata$survival_unit <- "month"
  coldata$survival_type <- "PFS"
  
  return(SummarizedExperiment(assays=list("expr"=assay), colData=coldata))
}


Create_CNA_SummarizedExperiment = function( cna_output , case_file , clin_file , cna_file , feat_snv , feat_cna , feat_cin , snv_bool , cna_bool ){
  case = read.csv( case_file , sep=";" , stringsAsFactors=FALSE )
  cna = read.csv( cna_file , sep=";" , stringsAsFactors=FALSE )
  
  clin = read.csv( clin_file , sep=";" , stringsAsFactors=FALSE )
  rownames(clin) = clin$patient
  clin = cbind( clin , matrix( NA , nrow=nrow(clin) , ncol=12 ) )
  colnames(clin) <- se_colnames
  
  if(snv_bool){ 
    clin[ rownames(feat_snv) , colnames(feat_snv) ] = feat_snv
  }
  rownames(feat_cin) = feat_cin$patient
  clin[ rownames(feat_cin) , "CIN" ] = feat_cin[ , "CIN" ]
  
  rownames(feat_cna) = feat_cna$patient
  clin[ rownames(feat_cna) , colnames(feat_cna)[-1] ] = feat_cna[ , c( "CNA_tot" , "AMP" , "DEL" ) ]
  clin = clin[ colnames(cna) , ]
  
  return(format_se(assay=cna, coldata=clin))
}

Create_EXP_SummarizedExperiment = function( case_file , clin_file , expr_file , feat_snv , feat_cna , feat_cin , snv_bool , cna_bool ){
  case = read.csv( case_file , sep=";" , stringsAsFactors=FALSE )
  expr = read.csv( expr_file , sep=";" , stringsAsFactors=FALSE )

  clin = read.csv( clin_file , sep=";" , stringsAsFactors=FALSE )
  rownames(clin) = clin$patient
  
  clin = cbind( clin , matrix( NA , nrow=nrow(clin) , ncol=12 ) )
  colnames(clin) <- se_colnames
  
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
  
  return(format_se(assay=expr, coldata=clin))
}

Create_SNV_SummarizedExperiment = function( snv_output , case_file , clin_file , snv_file , feat_snv , feat_cna , feat_cin , cna_bool, snv_bool ){
  case = read.csv( case_file , sep=";" , stringsAsFactors=FALSE )
  snv = read.csv( snv_file , sep=";" , stringsAsFactors=FALSE )
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
        if( !is.na( s$Gene[j]) ){ 
          if( !is.na( mat_snv[ s$Gene[j] , sample[i] ] ) ){
            mat_snv[ s$Gene[j] , sample[i] ]  = paste( mat_snv[ s$Gene[j] , sample[i] ] , paste( s$Ref[j] , s$Alt[j] , sep=">" ) , sep=";" )
          } else{
            mat_snv[ s$Gene[j] , sample[i] ]  = paste( s$Ref[j] , s$Alt[j] , sep=">" )					
          }
        }
      }
    }
  }
  
  clin = read.csv( clin_file , sep=";" , stringsAsFactors=FALSE )
  rownames(clin) = clin$patient
  clin = cbind( clin , matrix( NA , nrow=nrow(clin) , ncol=12 ) )
  colnames(clin) <- se_colnames
  clin[ rownames(feat_snv) , colnames(feat_snv) ] = feat_snv
  
  if(cna_bool){ 
    rownames(feat_cin) = feat_cin$patient
    clin[ rownames(feat_cin) , "CIN" ] = feat_cin[ , "CIN" ]
    rownames(feat_cna) = feat_cna$patient
    clin[ rownames(feat_cna) , colnames(feat_cna)[-1] ] = feat_cna[ , c( "CNA_tot" , "AMP" , "DEL" ) ]
  }
  
  clin = clin[ patient , ]

  return(format_se(assay=mat_snv, coldata=clin))
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
	}

	if( cna_bool ){
	  se_list[["cna"]] <- Create_CNA_SummarizedExperiment( case_file=case_file , clin_file=clin_file , cna_file=cna_file , 
									feat_snv=feat_snv , feat_cna=feat_cna , feat_cin=feat_cin , cna_bool=cna_bool , snv_bool=snv_bool )
	}

	if( expr_bool ){
	  se_list[["expr"]] <- Create_EXP_SummarizedExperiment( case_file=case_file , clin_file=clin_file , expr_file=expr_file , 
									feat_snv=feat_snv , feat_cna=feat_cna, feat_cin=feat_cin , cna_bool=cna_bool , snv_bool=snv_bool )
	}

	if( snv_bool ){
	  se_list[["snv"]] <- Create_SNV_SummarizedExperiment( case_file=case_file , clin_file=clin_file , snv_file=snv_file , 
									feat_snv=feat_snv , feat_cna=feat_cna , feat_cin=feat_cin , cna_bool=cna_bool , snv_bool=snv_bool )
	}
  return(se_list)
}
