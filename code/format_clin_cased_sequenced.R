library(data.table)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
work_dir <- args[1]

case <- read.csv( file.path(work_dir, 'cased_sequenced.csv') , sep=";" , stringsAsFactors=FALSE )
case$patient <- str_replace_all(case$patient, '[-\\.]', '_')
clin <- read.csv( file.path(work_dir, 'CLIN.csv') , sep=";" , stringsAsFactors=FALSE )
clin$patient <- str_replace_all(clin$patient, '[-\\.]', '_')
clin$tissueid[is.na(clin$tissueid)] <- ""
clin$treatmentid[is.na(clin$treatmentid)] <- ""

saveRDS(case, file.path(work_dir, 'cased_sequences.rds'))
saveRDS(case, file.path(work_dir, 'CLIN.rds'))