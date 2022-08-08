library(vcd)

work_dir <- "~/Documents/ICBCuration/duplication/results"
expr <- readRDS(file.path(work_dir, 'expr_list.rds'))

study = names(expr)

study <- c(
  "Jung", "Kim", "Liu", "Mariathasan", "Miao1", "Nathanson", "Padron", "Puch", "Riaz", "Roh", 
  "Shiuan", "Snyder", "Van_Allen", "VanDenEnde"
)

nb_before = NULL
for(i in 1:length(study) ){

	nb_before = c( nb_before , ncol( exprs(expr[[ study[i] ]] ) ) ) 
}
names(nb_before) = study

duplicated_patients <- list()

for( i in 1: ( length(study)-1 ) ){
  
  print(paste('Filtering ', study[i]))
  removed_patients = c()
  
	for( j in (i+1) : length(study) ) {
	  print(paste(study[i], 'Samples:', length(rownames(pData(expr[[ study[i] ]])))))
	  print(paste('--', study[j]))

		expr1 = exprs(expr[[ study[i] ]])
		expr2 = exprs(expr[[ study[j] ]])

		clin1 = pData(expr[[ study[i] ]])
		clin2 = pData(expr[[ study[j] ]])

		tumor1 = sort( unique( clin1$cancer_type ))

		expr_out = clin_out = NULL
		for( l  in 1:length(tumor1)){

			e1 = expr1[ , clin1$unique_patient_ID[ clin1$cancer_type %in% tumor1[l] ] ]
			c1 = clin1[ clin1$unique_patient_ID[ clin1$cancer_type %in% tumor1[l] ] , ]
			e2 = expr2[ , clin2$unique_patient_ID[ clin2$cancer_type %in% tumor1[l] ] ] 
			c2 = clin2[ clin2$unique_patient_ID[ clin2$cancer_type %in% tumor1[l] ] , ]

			if( ncol(e2) & !is.null( dim(e1) ) ){

				geneID = intersect( rownames(e1) , rownames(e2) )
				patient1 = colnames(e1)
				patient2 = colnames(e2)
				
				cor=NULL

				for( m in 1:length(patient1)){

					cor_p1 = NULL
					for( n in 1:length(patient2)){

						sex = ifelse( !is.na( c1[ patient1[m] , ]$sex ) & !is.na( c2[ patient2[n] , ]$sex ) , c1[ patient1[m] , ]$sex == c2[ patient2[n] , ]$sex , TRUE  )
						age = ifelse( !is.na( c1[ patient1[m] , ]$age ) & !is.na( c2[ patient2[n] , ]$age ) , abs( c1[ patient1[m] , ]$age - c2[ patient2[n] , ]$age ) <= 5 , TRUE  )
						
						if( sex & age ){

							cor.test( e1[ geneID , patient1[m] ] , e2[ geneID , patient2[n] ] , method='s' )$estimate 

						} else{

							cor_p1 = c( cor_p1 , 0 )
						}
					}
					cor = c( cor , max( cor_p1 , na.rm=TRUE ) )
				}

				names(cor) = patient1

				remove.patient = names( cor[ cor >= .98 ] )

				if( length( remove.patient ) ){
					print( paste( "Removed " , length( remove.patient ) , " (" , study[i] , "|" , tumor1[l] , ")/" , ncol(e1) , " (" , study[j] , ")" , sep="" ))
		  		flush.console()
		  		removed_patients <- c(removed_patients, remove.patient)
				}	
			}
		}
	}
  if( length(removed_patients) > 0 ){
    duplicated_patients[[study[i]]] <- removed_patients
  }	
}

saveRDS(duplicated_patients, file.path(work_dir, 'duplicated_patients_expr.rds'))

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

write.table( nb , file=file.path(work_dir, "duplicated_filtered_expr.csv") , quote=FALSE , sep=";" , col.names=TRUE , row.names=TRUE )

