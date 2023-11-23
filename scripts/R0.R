#====================

args <- commandArgs()
organism <- args[6]
sample <- args[7]

#====================

# Step 1: Load and Merge Bedgraph Files
merge_bedgraph_file <- function(forward_file, reverse_file, output_filename) {
  
  # Detect the separator in the bedgraph file
  forward_separator <- readLines(forward_file, n = 1)
  if (grepl("\t", forward_separator)) {
    forward_separator <- "\t"
  } else if (grepl(" ", forward_separator)) {
    forward_separator <- " "
  } else {
    stop("Unable to detect separator in the file.")
  }
  
  # Read the bedgraph file into a data frame
  forward_data <- read.table(forward_file, sep = forward_separator, header = FALSE,
                         col.names = c("chrom", "start", "end", "forward_intensity"))
  
  # Detect the separator in the bedgraph file
  reverse_separator <- readLines(reverse_file, n = 1)
  if (grepl("\t", reverse_separator)) {
    reverse_separator <- "\t"
  } else if (grepl(" ", reverse_separator)) {
    reverse_separator <- " "
  } else {
    stop("Unable to detect separator in the file.")
  }
  
  # Read the bedgraph file into a data frame
  reverse_data <- read.table(reverse_file, sep = reverse_separator, header = FALSE,
                         col.names = c("chrom", "start", "end", "reverse_intensity"))
  
  # Merge the two bedgraphs
  merged_data <- merge(forward_data, reverse_data, by=c("chrom", "start", "end"))
  
  # Calculate the difference in intensities to find divergent regions
  merged_data$divergent_intensity <- abs(merged_data$forward_intensity - merged_data$reverse_intensity)
  
  # Sort the dataframe
  merged_data <- merged_data[order(merged_data$chrom, merged_data$start), ]
  
  # Save the filtered data as a text file
  write.table(merged_data[, c("chrom", "start", "end", "divergent_intensity")], 
              file = output_filename, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}

bedgraph_pl <- paste0(organism, "/bigWig/", sample, "_allMap_unnorm_pl.bedgraph")
bedgraph_mn <- paste0(organism, "/bigWig/", sample, "_allMap_unnorm_mn.bedgraph")
output_folder <- paste0(organism, "/analysis/functionalGenomics_", sample)

# Process and merge the bedgraph files
merge_bedgraph_file(bedgraph_pl, bedgraph_mn, paste0(output_folder, "/merged.bedgraph"))

