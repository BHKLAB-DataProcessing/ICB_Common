Get_Response = function( data ){

	response = rep( NA , nrow(data) )
	names( response) = rownames(data)

	for( i in 1:nrow(data)){

		other.response = as.character( data$response.other.info[i] )
		recist = as.character( data$recist[i] )
		t.pfs = as.numeric( as.character( data$t.pfs[i] ))
		pfs = as.numeric( as.character( data$pfs[i] ))

		if( is.na( recist ) ){
			if( !is.na( other.response ) ){
					response[i] = other.response
			} else{
				if( !is.na( t.pfs ) ){
					if( t.pfs <= 6 & pfs %in% 1  ){
						response[i] = "NR"
					} else{
						if( t.pfs >= 6 & pfs %in% 0 ){
							response[i] =  "R"
						}
					}
				}				
			}
		} else{
			if( recist %in% c( "CR" , "PR" ) ){
				response[i] = "R"
			}
			if( recist %in%  "SD" ){
				if( !is.na( t.pfs ) ){
					if( t.pfs <= 6 & pfs %in% 1 ){
						response[i] = "NR"
					} else{
						if( t.pfs >= 6 & pfs %in% 0 ){
							response[i] = "R"
						} 				
					}
				}
			}
			if( recist %in%  "PD" ){
				response[i] = "NR"
			}			
		}
	}
	response
}