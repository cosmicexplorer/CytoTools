#! /usr/local/bin/Rscript

library(emdist)

# read configuration values
source( "config.r")

read_sample_file_to_dataframe <- function ( sample_file ) {
 
  print( sprintf("Loading %s", sample_file))
  
  #figure out which column has the tSNE1 and tSNE2 data, don't assume it because the column number isn't always the same.
  fileColumnsFrame <- read.table(sample_file, sep="\t", header = T, skip = 1, comment.char = "", nrows = 1)
  tSNEcolumns <- grep( "tSNE.", colnames(fileColumnsFrame))
  
  columnFilter <- rep( "NULL", ncol(fileColumnsFrame))
  columnFilter[tSNEcolumns] <- "numeric"
  
  return( read.table(sample_file, sep="\t", header = T, skip = 1, comment.char = "", colClasses = columnFilter))
}


# Making the EMD world (2D matrix) from the tSNE data frame
populate_matrix_from_frame <- function( dimension, fromFrame ) {

  v <- vector( mode="numeric", length = ( dimension^2 ))
  fillMatrix <- matrix( data=v, nrow = dimension)
  
  midpoint <- round(nrow(fillMatrix)/2)
  
  fromFrame <- round(fromFrame) + midpoint

  # should be able to do this without the loop somehow...
  for( selected_row in 1 : nrow(fromFrame)) {
    fillMatrix[ fromFrame$tSNE1[selected_row], fromFrame$tSNE2[selected_row] ] <- fillMatrix[ fromFrame$tSNE1[selected_row], fromFrame$tSNE2[selected_row] ] + 1
  }
  
  return( fillMatrix )
}

# Used to convert the binned matrix back to a frame for plotting
# This is for graphical and testing purposes for now... 
create_frame_from_matrix <- function( fromMatrix ) {
  
  matrix_row <- row(fromMatrix)[which(!fromMatrix == 0)]
  matrix_col <- col(fromMatrix)[which(!fromMatrix == 0)]
  matrix_vals <- fromMatrix[cbind(matrix_row, matrix_col)]
    
  return ( data.frame( matrix_row, matrix_col, matrix_vals ))
}

plot_binned_frame <- function( binnedFrame, my_color, my_ylim, my_xlim ) {

  symbols(  binnedFrame[[1]], 
            binnedFrame[[2]], 
            sqrt(binnedFrame[[3]]/pi), 
            inches=1/3, 
            ann=F, 
            bg=my_color, 
            fg=NULL, 
            xlim = my_xlim,
            ylim=my_ylim ) 
  
}

graph_binned_plots <- function (matrix1, 
                                matrix2, 
                                sample1_file, 
                                sample2_file, 
                                raw_tsne1_frame,
                                raw_tsne2_frame,
                                output_graphic_file,
                                results) {
  
  if( is.null(output_graphic_file)) {
    return()
  }

  png(output_graphic_file, width = 1200, height = 800)
  par(mfrow=c(2,2))

  smoothScatter(raw_tsne1_frame, colramp = colorRampPalette(c("white", "steelblue2")))
  
  smoothScatter(raw_tsne2_frame, colramp = colorRampPalette(c("white", "green")))
  
  my_ylim <- c(0,ncol(matrix2))
  my_xlim <-  c(0, nrow(matrix2))
  
  plot_binned_frame( create_frame_from_matrix(matrix2), "green", my_ylim, my_xlim)
  par(new=T)
  plot_binned_frame( create_frame_from_matrix(matrix1), "steelblue2", my_ylim, my_xlim)
  title( main = sprintf( "%s (blue)", sample1_file ), 
         sub = sprintf("Bin Factor = %d, Max Iterations=%d", as.integer( results["binning_factor"]), as.integer(results["max_iterations"])))

  plot_binned_frame( create_frame_from_matrix(matrix1), "steelblue2", my_ylim, my_xlim)
  par(new=T)
  plot_binned_frame( create_frame_from_matrix(matrix2), "green", my_ylim, my_xlim)
  title( main = sprintf( "%s (green)", sample2_file ), 
         sub = sprintf("EMD energy = %f", as.numeric(results["result"])))

  dev.off()
}

#get_max_dimension <- function( aFrame ) {
#  max_dimension <-  2 * ceiling(max(abs(range(aFrame)))) + 1
#}

execute_emd <- function( sample1_file, sample2_file, output_graphic_file = NULL) {

  process_start_time = date()
  t1 <- proc.time()
  
  rawsample1frame <- read_sample_file_to_dataframe( sample1_file )
  rawsample2frame <- read_sample_file_to_dataframe( sample2_file )

  # adjust for binning factor...
  binfactoredsample1frame <- rawsample1frame / binning_factor
  binfactoredsample2frame <- rawsample2frame / binning_factor
  
  # get the biggest dimension in both frames considering x and y  
  matrix_size <-  2 * ceiling( max(abs(range( c( binfactoredsample1frame, binfactoredsample2frame) )))) + 1

  print( sprintf("matrix dimension is %f", matrix_size) )
  print( "Converting sample 1 to EMD matrix...")
  matrix1 <- populate_matrix_from_frame(matrix_size, binfactoredsample1frame)
  normalized_matrix1 <- matrix1 / sum(matrix1)

  print( "Converting sample 2 to EMD matrix...")
  matrix2 <- populate_matrix_from_frame(matrix_size, binfactoredsample2frame)
  normalized_matrix2 <- matrix2 / sum(matrix2)

  print("Running EMD on prepared data...")  
  minimum_energy <- emd2d(normalized_matrix1, normalized_matrix2, max.iter=max_iterations, xdist=binning_factor, ydist=binning_factor)

  results <- c( file1 = sample1_file, 
                file2 = sample2_file, 
                result = minimum_energy,
                max_iterations = max_iterations,
                binning_factor = binning_factor,
                process_start_time = process_start_time,
                processing_duration_seconds = (proc.time() - t1)["elapsed"],
                graphic_output_file = output_graphic_file,
                machine = Sys.info()["nodename"]
  )

  graph_binned_plots(matrix1, matrix2, sample1_file, sample2_file, rawsample1frame, rawsample2frame, output_graphic_file, results)
  
  return( results )
}

# ways to graph the original vsSNE data 
# plot(frame1)
# smoothScatter(frame2)

# After converting to bins, and using create_frame_from_matrix() to to back to a frame from the binned matrix
# use the following to show views of the binned data.

#fback <- create_frame_from_matrix(matrix1)
# plot(fback$matrix_row, fback$matrix_col, pch=16, cex=2)
# with(fback, symbols(fback$matrix_row, fback$matrix_col, sqrt(fback$matrix_vals/pi), inches=1/3, ann=F, bg="steelblue2", fg=NULL) )
# plot_binned_frame( create_frame_from_matrix(matrix1), "steelblue2", c(0,ncol(matrix1))

# Keep a graph and add to it with next command
# par(new=T)

