library(data.table)
library(stringr)

Get_TMB_raw = function( file , cases ){
  snv = read.csv( file , sep=";" , stringsAsFactors=FALSE )
  patient = cases[ cases$snv %in% 1, ]$patient
  
  tmb = rep( 0 , length(patient))
  names(tmb) = as.character( patient )
  
  snv_patient = sort(unique(snv$Sample))
  for(i in 1:length(snv_patient)){
    s = snv[ snv$Sample %in% snv_patient[i] , ]
    tmb[ as.character( snv_patient[i] ) ] = nrow(s)
  }
  tmb
}

Get_nsTMB_raw = function( file , cases ){
  snv = read.csv(  file , sep=";" , stringsAsFactors=FALSE )
  patient = cases[ cases$snv %in% 1, ]$patient
  
  tmb = rep( 0 , length(patient))
  names(tmb) = as.character( patient )
  
  snv_patient = sort(unique(snv$Sample))
  for(i in 1:length(snv_patient)){
    s = snv[ snv$Sample %in% snv_patient[i] & 
               snv$Effect %in% c("In_Frame_Del" , "In_Frame_Ins" , "Start_Codon_Ins" , "Frame_Shift_Del" ,
                                 "Frame_Shift_Ins" , "Missense_Mutation" , "Nonsense_Mutation" , "Nonstop_Mutation" ,
                                 "Splice_Site" , "Stop_Codon_Del" , "De_novo_Start_OutOfFrame" , "Start_Codon_SNP") ,
    ]
    tmb[ as.character( snv_patient[i] ) ] = nrow(s)
  }
  tmb
}

Get_indel_TMB_raw = function( file , cases , indel_bool ){
  snv = read.csv( file , sep=";" , stringsAsFactors=FALSE )
  patient = cases[ cases$snv %in% 1, ]$patient
  snv_patient = sort(unique(snv$Sample))
  if(indel_bool){
    
    tmb = rep( 0 , length(patient))
    names(tmb) = as.character( patient )
    
    # snv_patient = sort(unique(snv$Sample))
    for(i in 1:length(snv_patient)){
      s = snv[ snv$Sample %in% snv_patient[i] & snv$MutType %in% "INDEL" , ]
      tmb[ as.character( snv_patient[i] ) ] = nrow(s)
    }
    
  } else{
    # tmb = rep( NA , length(patient))
    # names(tmb) = as.character( patient )
    tmb = rep( NA , length(snv_patient))
    names(tmb) = as.character( snv_patient )
  }
  tmb
}

Get_indel_nsTMB_raw = function( file , cases  , indel_bool ){
  snv = read.csv( file , sep=";" , stringsAsFactors=FALSE )
  patient = cases[ cases$snv %in% 1, ]$patient
  snv_patient = sort(unique(snv$Sample))
  if(indel_bool){
    
    tmb = rep( 0 , length(patient))
    names(tmb) = as.character( patient )
    
    # snv_patient = sort(unique(snv$Sample))
    for(i in 1:length(snv_patient)){
      s = snv[ 	snv$Sample %in% snv_patient[i] & 
                  snv$MutType %in% "INDEL" & 
                  snv$Effect %in% c("In_Frame_Del","In_Frame_Ins","Start_Codon_Ins","Frame_Shift_Del",
                                    "Frame_Shift_Ins","Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation",
                                    "Splice_Site","Stop_Codon_Del") , ]
      tmb[ as.character( snv_patient[i] ) ] = nrow(s) 
    }
  } else{
    # tmb = rep( NA , length(patient))
    # names(tmb) = as.character( patient )
    tmb = rep( NA , length(snv_patient))
    names(tmb) = as.character( snv_patient )
  }
  tmb
}

Get_TMB_perMb = function( file , cases , coverage ){
  snv = read.csv( file , sep=";" , stringsAsFactors=FALSE )
  patient = cases[ cases$snv %in% 1, ]$patient
  
  tmb = rep( 0 , length(patient))
  names(tmb) = as.character( patient )
  
  snv_patient = sort(unique(snv$Sample))
  for(i in 1:length(snv_patient)){
    s = snv[ snv$Sample %in% snv_patient[i] , ]
    tmb[ as.character( snv_patient[i] ) ] = nrow(s) / coverage
  }
  tmb	
}

Get_nsTMB_perMb = function( file , cases , coverage ){
  snv = read.csv( file , sep=";" , stringsAsFactors=FALSE )
  patient = cases[ cases$snv %in% 1, ]$patient
  
  tmb = rep( 0 , length(patient))
  names(tmb) = as.character( patient )
  
  snv_patient = sort(unique(snv$Sample))
  for(i in 1:length(snv_patient)){
    s = snv[ snv$Sample %in% snv_patient[i] & 
               snv$Effect %in% c("In_Frame_Del","In_Frame_Ins","Start_Codon_Ins","Frame_Shift_Del",
                                 "Frame_Shift_Ins","Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation",
                                 "Splice_Site","Stop_Codon_Del") ,
    ]
    tmb[ as.character( snv_patient[i] ) ] = nrow(s) / coverage
  }
  tmb
}

Get_indel_TMB_perMb = function( file , cases , coverage , indel_bool ){
  snv = read.csv( file , sep=";" , stringsAsFactors=FALSE )
  patient = cases[ cases$snv %in% 1, ]$patient
  snv_patient = sort(unique(snv$Sample))
  if(indel_bool){
    
    tmb = rep( 0 , length(patient))
    names(tmb) = as.character( patient )
    
    # snv_patient = sort(unique(snv$Sample))
    for(i in 1:length(snv_patient)){
      s = snv[ snv$Sample %in% snv_patient[i] & snv$MutType %in% "INDEL" , ]
      tmb[ as.character( snv_patient[i] ) ] = nrow(s) / coverage
    }
  } else{
    # tmb = rep( NA , length(patient))
    # names(tmb) = as.character( patient )
    tmb = rep( NA , length(snv_patient))
    names(tmb) = as.character( snv_patient )
  }
  tmb
}

Get_indel_nsTMB_perMb = function( file , cases , coverage , indel_bool ){
  snv = read.csv( file , sep=";" , stringsAsFactors=FALSE )
  patient = cases[ cases$snv %in% 1, ]$patient
  snv_patient = sort(unique(snv$Sample))
  if(indel_bool){
    tmb = rep( 0 , length(patient))
    names(tmb) = as.character( patient )
    
    # snv_patient = sort(unique(snv$Sample))
    for(i in 1:length(snv_patient)){
      s = snv[ 	snv$Sample %in% snv_patient[i] & 
                  snv$MutType %in% "INDEL" & 
                  snv$Effect %in% c("In_Frame_Del","In_Frame_Ins","Start_Codon_Ins","Frame_Shift_Del",
                                    "Frame_Shift_Ins","Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation",
                                    "Splice_Site","Stop_Codon_Del") , ]
      tmb[ as.character( snv_patient[i] ) ] = nrow(s) / coverage
    }
  } else{
    # tmb = rep( NA , length(patient))
    # names(tmb) = as.character( patient )
    tmb = rep( NA , length(snv_patient))
    names(tmb) = as.character( snv_patient )
  }
  tmb
}

args <- commandArgs(trailingOnly = TRUE)
work_dir <- args[1]
coverage <- args[2]
indel_bool <- args[3]

coverage <- as.numeric(coverage)
file <- file.path(work_dir, 'SNV.csv')
cases <- readRDS(file.path(work_dir, 'cased_sequenced.rds'))

feat_snv = cbind( 
  Get_TMB_raw( file=file , cases=cases ) ,
  Get_nsTMB_raw( file=file , cases=cases ) ,
  Get_indel_TMB_raw( file=file , cases=cases , indel_bool=indel_bool ) ,
  Get_indel_nsTMB_raw( file=file , cases=cases , indel_bool=indel_bool ) ,
  Get_TMB_perMb( file=file , cases=cases , coverage=coverage ) ,
  Get_nsTMB_perMb( file=file , cases=cases , coverage=coverage) ,
  Get_indel_TMB_perMb( file=file , cases=cases , coverage=coverage , indel_bool=indel_bool ) ,
  Get_indel_nsTMB_perMb( file=file , cases=cases , coverage=coverage , indel_bool=indel_bool )
)
colnames(feat_snv) = c( "TMB_raw" , "nsTMB_raw" , "indel_TMB_raw" , "indel_nsTMB_raw" , 
                 "TMB_perMb" , "nsTMB_perMb" , "indel_TMB_perMb" , "indel_nsTMB_perMb" )

rownames(feat_snv) <- str_replace_all(rownames(feat_snv), '[-\\.]', '_')
saveRDS(feat_snv, file.path(work_dir, 'feat_snv.rds'))