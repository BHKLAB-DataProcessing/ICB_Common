library(readxl) 
library(data.table)

format_excel_and_save <- function(input_path, output_path, sheetname, compress=FALSE, header_row_num=1){
  df <- read_excel(input_path, sheet=sheetname)
  colnames(df) <- df[header_row_num, ]
  df <- df[-header_row_num, ]
  if(compress){
    gz <- gzfile(output_path, "w")
    write.table( df , file=gz , quote=FALSE , sep="\t" , col.names=TRUE , row.names=FALSE )
    close(gz)
  }else{
    write.table( df , file=output_path , quote=FALSE , sep="\t" , col.names=TRUE , row.names=FALSE )
  }
  rm(df)
}

read_and_format_excel <- function(input_path, sheetname, header_row_num=1){
  df <- read_excel(input_path, sheet=sheetname)
  colnames(df) <- df[header_row_num, ]
  df <- df[-header_row_num, ]
  return(df)
}