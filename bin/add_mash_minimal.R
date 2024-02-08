message('=== Adding Mash group to report ===')
in_report_path <- commandArgs(TRUE)[1]
threshold <- as.numeric(commandArgs(TRUE)[2])

report_types <-
  rep('character', 6)
report <- 
  read.delim(in_report_path, sep = '\t', 
             header = FALSE, colClasses = report_types)

if (nrow(report)[1] > 1) {
  stop('multiple result are not supported with the minimal option')
}

mash_result_path <- paste0(dirname(in_report_path),'/', report$V6)

mash_names <- 
  c('identity', 'shared_hashes', 'median_multiplicity','p_value', 'query_ID', 'query_comment')
mash_types <- 
  c('numeric', 'character', 'numeric', 'numeric', 'character', 'character')
mash_result <- 
  read.delim(mash_result_path, sep = '\t', header = FALSE, 
             col.names = mash_names, colClasses = mash_types)


if(max(mash_result$identity) < threshold){
  mash_group <- 'Undetermined'
} else {

  mash_group <- mash_result$query_ID[order(mash_result$identity, decreasing = TRUE)][1]
  
  mash_group <- gsub(pattern = '.*_(.*)_.*', replacement = '\\1', x = mash_group)
}


report$V6 <- mash_group

write.table(x = report, file = in_report_path, sep = '\t', 
            col.names = FALSE, row.names = FALSE, quote = FALSE)


