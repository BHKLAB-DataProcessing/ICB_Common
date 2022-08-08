library(vcd)

work_dir <- "~/Documents/ICBCuration/duplication/results"

snv <- readRDS(file.path(work_dir, 'snv_list.rds'))

study = names(snv)
nb_before = NULL
for(i in 1:length(study) ){
	nb_before = c( nb_before , ncol( exprs(snv[[ study[i] ]] ) ) ) 
}
names(nb_before) = study

duplicated_patients <- list()

for( i in 1: ( length(study)-1 ) ){
  
  print(paste('Filtering ', study[i]))
  removed_patients = c()
	for( j in (i+1) : length(study) ) {
    print(paste(study[i], 'Samples:', length(rownames(pData(snv[[ study[i] ]])))))
	  print(paste('--', study[j]))
	  
		snv1 = exprs(snv[[ study[i] ]])
		snv2 = exprs(snv[[ study[j] ]])

		clin1 = pData(snv[[ study[i] ]])
		clin2 = pData(snv[[ study[j] ]])

		tumor1 = sort( unique( clin1$cancer_type ))

		for( l  in 1:length(tumor1)){
			s1 = snv1[ , clin1$unique_patient_ID[ clin1$cancer_type %in% tumor1[l] ] ]
			c1 = clin1[ clin1$unique_patient_ID[ clin1$cancer_type %in% tumor1[l] ] , ]
			s2 = snv2[ , clin2$unique_patient_ID[ clin2$cancer_type %in% tumor1[l] ] ] 
			if( is.null( dim(s2) )){ 
				s2 = as.matrix( s2 ) 
				colnames(s2) = clin2$unique_patient_ID[ clin2$cancer_type %in% tumor1[l] ]
			}
			c2 = clin2[ clin2$unique_patient_ID[ clin2$cancer_type %in% tumor1[l] ] , ]
      
			print(paste('Tumor:', tumor1[l], 'Subset Num:', length(colnames(s1))))
			
			if( ncol(s2) & !is.null( dim(s1) ) ){

				geneID = intersect( rownames(s1) , rownames(s2) )
				patient1 = colnames(s1)
				patient2 = colnames(s2)
				
				kappa=NULL

				for( m in 1:length(patient1)){

					kappa_p1 = NULL
					for( n in 1:length(patient2)){

						sex = ifelse( !is.na( c1[ patient1[m] , ]$sex ) & !is.na( c2[ patient2[n] , ]$sex ) , c1[ patient1[m] , ]$sex == c2[ patient2[n] , ]$sex , TRUE  )
						age = ifelse( !is.na( c1[ patient1[m] , ]$age ) & !is.na( c2[ patient2[n] , ]$age ) , abs( c1[ patient1[m] , ]$age - c2[ patient2[n] , ]$age ) <= 5 , TRUE  )
						
						if( sex & age ){
							x = s1[ geneID , patient1[m] ]
							y = s2[ geneID , patient2[n] ]

							level = sort(unique(c( x , y )))
							x = factor( x , levels= level )
							y = factor( y , levels= level ) 

							kappa_p1 = c( kappa_p1 , ifelse( sum( table( x , y ) ) >= 10 , Kappa( table( x , y ) )$Unweighted[1] , 0 ) )
						} else{
							kappa_p1 = c( kappa_p1 , 0 )
						}
					}

					kappa = c( kappa , max( kappa_p1 , na.rm=TRUE ) )
				}

				names(kappa) = patient1

				remove.patient = names( kappa[ kappa >= .9 ] )

				if( length( remove.patient ) ){
					print( paste( "Removed " , length( remove.patient) , " (" , study[i] , "|" , tumor1[l] , ")/" , ncol(s1) , " (" , study[j] , ")" , sep="" ))
		  		flush.console()
          removed_patients <- c(removed_patients, remove.patient)
				} 
			} 
		}

		if( length(removed_patients) > 0 ){
		  duplicated_patients[[study[i]]] <- removed_patients
		}	
	}
}

saveRDS(duplicated_patients, file.path(work_dir, 'duplicated_patients_snv.rds'))

nb_after = NULL
for(i in 1:length(study) ){
  duplicated_num <- 0
  if(study[i] %in% names(duplicated_patients)){
    duplicated_num <- length(duplicated_patients[[study[i]]])
  }
	nb_after = c( nb_after , nb_before[[ study[i] ]] - duplicated_num ) 
}
names(nb_after) = study

nb = cbind(nb_before , nb_after )
colnames(nb) = c( "before" , "after" )
rownames(nb) = study

write.table( nb , file=file.path(work_dir, "duplicated_filtered_snv.csv") , quote=FALSE , sep=";" , col.names=TRUE , row.names=TRUE )
