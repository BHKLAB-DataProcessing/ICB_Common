library(data.table)

Get_CIN = function( case_file, cna_file , cin_bool ){
  cases = read.csv( case_file , sep=";" , stringsAsFactors=FALSE )
  patient = cases[ cases$cna %in% 1, ]$patient
  
  if(cin_bool){
    cna = read.csv( cna_file , sep=";" , stringsAsFactors=FALSE )
    cin = rep( 0 , length(patient))
    names(cin) = as.character( patient )
    
    cna_patient = sort(unique(cna$sample))
    for(i in 1:length(cna_patient)){
      c = cna[ cna$sample %in% cna_patient[i] , ]
      cin[ cna_patient[i] ] <- sum( c[ abs( round( c$seg.mean , 2 ) ) >= 0.2 , ]$loc.end - c[ abs( round( c$seg.mean , 2 ) ) >= 0.2 , ]$loc.start ) / sum( c$loc.end - c$loc.start )
    }
  } else{
    cin = rep( NA , length(patient))
    names(cin) = as.character( patient )
  }
  
  cin  = as.data.frame( cbind( names(cin) , cin ) )
  colnames(cin) = c( "patient" , "CIN" )
  
  return(cin)
  
}

Get_CNA_feature = function( case_file, cna_file, coverage , cna_bool ){
  cases = read.csv( case_file , sep=";" , stringsAsFactors=FALSE )
  patient = cases[ cases$cna %in% 1, ]$patient
  
  if(cna_bool){
    cna = read.csv( cna_file , sep=";" , stringsAsFactors=FALSE )
    
    # CNA total
    cna_tot = cna
    cna_tot[ cna_tot != 0 ] = 1
    cna_tot = colSums( cna_tot ) / coverage	
    
    # AMP total
    amp = cna
    amp[ amp < 0 ] = 0
    amp[ amp > 0 ] = 1
    amp = colSums( amp ) / coverage
    
    # DEL total
    del = cna
    del[ del > 0 ] = 0
    del[ del < 0 ] = 1
    del = colSums( del ) / coverage
    
    # Output
    out = cbind( names(cna_tot) , cna_tot , amp , del )
    colnames(out) = c( "patient" , "CNA_tot" , "AMP" , "DEL" )
    
    return(as.data.frame(out))
  }
}