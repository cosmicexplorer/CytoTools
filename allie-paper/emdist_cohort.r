#! /usr/local/bin/Rscript
source(file="flow_emd.r")

#To run from the command line, pass an argument that is the path to the cohort with patients as sub directories
# containing sample files.
args = commandArgs(trailingOnly=TRUE)
cohort_dir <- args[1]
output_file <- paste(cohort_dir,"/cohort.csv", sep="")
graphics_output_dir <- paste( cohort_dir, "/graphics", sep="")

count = 0;

run_emd <- function( sample1, sample2 ) {

  count <<- count + 1
  
  output_graphic_file <-  paste( graphics_output_dir, "/pairing_", count, ".png", sep="")
  
  results <- execute_emd( sample1, sample2, output_graphic_file )

  write.table( data.frame(as.list(results)), file = output_file, append = T, sep = ",", col.names = F, row.names = F) 
  
  return( results["result"])

}

#clear out the csv we're going to build
unlink(paste(cohort_dir, "/cohort.csv", sep=""))

# setup the graphics output directory
dir.create(graphics_output_dir, showWarnings = FALSE)
unlink( paste(graphics_output_dir, "/*.png", sep=""))
# for each sample
samples <- list.files(path = cohort_dir, pattern = "*.txt", full.names = T)

matrix_with_results <- outer( samples, samples, FUN= Vectorize(run_emd) )

print( sprintf("results are stored in %s", output_file))

