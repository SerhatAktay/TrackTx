
#============================================================

args <- commandArgs()
organism <- args[6]
sample <- args[7]

path_to_bedgraphs <- paste0(organism, "/bigWig/")

#============================================================

## Required packages
packages = c("data.table", "dplyr", "progress")

suppressPackageStartupMessages({
  for (pkg in packages) {
    library(pkg, character.only = TRUE)
  }
})

## Load or install & load all
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

#-------------------------------------------------------------------

# Set the divergent transcription nt_window to a default value
nt_window <- 800

# Ask the user if they want to change the window value
cat("The default nt_window value is set to:", nt_window, "\n")
change_nt_window <- readline(prompt = "Do you want to change the nucleotide distance? (y/n): ")

# If the user wants to change the nt_window value, prompt them for a new value
if (tolower(change_nt_window) == "y") {
  nt_window <- as.numeric(readline(prompt = "Please enter the new value for nucleotide distance between the forward and reverse strand: "))
}

#-------------------------------------------------------------------

# Read the input files
positive_strand.data <- read.table(paste0(path_to_bedgraphs, sample, "_allMap_unnorm_pl.bedgraph"), header = FALSE, col.names = c("chromosome", "start", "end", "polymerase_signal"))
negative_strand.data <- read.table(paste0(path_to_bedgraphs, sample, "_allMap_unnorm_mn.bedgraph"), header = FALSE, col.names = c("chromosome", "start", "end", "polymerase_signal"))

#-------------------------------------------------------------------

# Convert the data to data.table format
positive_strand <- setDT(positive_strand.data)
negative_strand <- setDT(negative_strand.data)

# Sort the data by chromosome and start coordinate
setorder(positive_strand, chromosome, start)
setorder(negative_strand, chromosome, start)

#-------------------------------------------------------------------

# Function to identify peaks in polymerase engagement
find_peaks <- function(data, threshold = 2) {
  data[, abs_signal := abs(polymerase_signal)]
  peaks <- data[abs_signal >= threshold]
  setorder(peaks, -abs_signal)
  peaks[, abs_signal := NULL]  # Remove the temporary column
  return(peaks)
}

#-------------------------------------------------------------------

# Identify peaks in both strands
positive_peaks <- find_peaks(positive_strand)
negative_peaks <- find_peaks(negative_strand)

#-------------------------------------------------------------------

if (organism

# Function to find divergent transcription patterns within a predefined organism dependent nucleotide window
find_divergent_transcription <- function(pos_peaks, neg_peaks, max_window = nt_window) {
  divergent_regions <- data.table(chromosome = character(), start = integer(), 
                                  end = integer(), total_signal = numeric())
  
  # Iterate through chromosomes
  for (chr in unique(pos_peaks$chromosome)) {
    pos_chr <- pos_peaks[chromosome == chr]
    neg_chr <- neg_peaks[chromosome == chr]
    
    setkey(pos_chr, chromosome, start, end)
    setkey(neg_chr, chromosome, start, end)
    
    # Ensure the columns are in the correct order for foverlaps
    neg_chr[, end := start + max_window]
    
    joined <- foverlaps(neg_chr, pos_chr, by.x = c("chromosome", "start", "end"), by.y = c("chromosome", "start", "end"), nomatch = 0)
    
    # Filter to ensure the peaks are within a predefined nucleotide window
    joined <- joined[abs(start - i.start) <= max_window]
    
    if (nrow(joined) > 0) {
      joined[, `:=`(region_start = pmin(i.start, start),
                    region_end = pmax(i.end, end))]
      joined[, total_signal := sum(abs(polymerase_signal)) + sum(abs(i.polymerase_signal)), by = .(region_start, region_end)]
      divergent_regions <- rbind(divergent_regions, joined[, .(chromosome, start = region_start, end = region_end, total_signal)])
    }
  }
  
  setorder(divergent_regions, chromosome, start)
  return(divergent_regions)
}

#-------------------------------------------------------------------

# Find regions of divergent transcription
divergent_transcription <- find_divergent_transcription(positive_peaks, negative_peaks)

#-------------------------------------------------------------------

merge_overlapping_regions <- function(regions) {
  # Ensure the regions are ordered by chromosome and start position
  setorder(regions, chromosome, start)
  
  # Initialize variables to store merged regions
  merged_regions <- list()
  current_index <- 1
  
  # Initialize progress bar
  pb <- progress_bar$new(total = nrow(regions), format = "  Merging regions [:bar] :percent eta: :eta", clear = FALSE)
  
  # Iterate through regions
  while (current_index <= nrow(regions)) {
    # Initialize the current region
    current_chromosome <- regions$chromosome[current_index]
    current_start <- regions$start[current_index]
    current_end <- regions$end[current_index]
    current_total_signal <- regions$total_signal[current_index]
    
    # Initialize count of merged regions
    merged_count <- 1
    
    # Find the end of the merged region
    next_index <- current_index + 1
    while (next_index <= nrow(regions) &&
           regions$chromosome[next_index] == current_chromosome &&
           regions$start[next_index] <= current_end + 1) {
      current_end <- max(current_end, regions$end[next_index])
      current_total_signal <- current_total_signal + regions$total_signal[next_index]
      next_index <- next_index + 1
      merged_count <- merged_count + 1
    }
    
    # Save the merged region only if at least 3 regions are merged
    if (merged_count >= 3) {
      merged_regions[[length(merged_regions) + 1]] <- list(
        chromosome = current_chromosome,
        start = current_start,
        end = current_end,
        total_signal = current_total_signal
      )
    }
    
    # Update the progress bar
    pb$tick(next_index - current_index)
    
    # Move to the next non-overlapping region
    current_index <- next_index
  }
  
  # Convert to data.table
  if (length(merged_regions) > 0) {
    merged_regions <- rbindlist(merged_regions)
  } else {
    merged_regions <- data.table(chromosome = character(), start = integer(), end = integer(), total_signal = numeric())
  }
  
  return(merged_regions)
}

#-------------------------------------------------------------------

# Merge overlapping regions
merged_divergent_transcription <- merge_overlapping_regions(divergent_transcription)

#-------------------------------------------------------------------

# Specify the file path where you want to save the BED file
output_folder = paste0(organism, "/analysis/functionalGenomics_", sample)

# Write the data to file without headers
write.table(merged_divergent_transcription, file = paste0(output_folder, "/divergent_transcription.bed"), quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")

save.image()
q()
y
