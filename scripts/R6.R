
#============================================================

args <- commandArgs()
organism <- args[6]
sample <- args[7]
nt_window <- args[8]

path_to_bedgraphs <- paste0(organism, "/bigWig/")

#============================================================

## Required packages
packages = c("data.table", "dplyr")

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

# Read the input files
positive_strand.data <- read.table(paste0(path_to_bedgraphs, sample, "_allMap_unnorm_pl.bedgraph"), header = FALSE, col.names = c("chromosome", "start", "end", "polymerase_signal"))
negative_strand.data <- read.table(paste0(path_to_bedgraphs, sample, "_allMap_unnorm_mn.bedgraph"), header = FALSE, col.names = c("chromosome", "start", "end", "polymerase_signal"))

#-------------------------------------------------------------------

# Convert the data to data.table format
positive_strand <- setorder(setDT(positive_strand.data))
negative_strand <- setorder(setDT(negative_strand.data))

#-------------------------------------------------------------------

# Function to identify peaks in polymerase engagement
find_peaks <- function(data, threshold, signal_sum_threshold) {
  
  # Step 1: Preprocessing - Calculate absolute value of polymerase_signal
  data[, abs_signal := abs(polymerase_signal)]
  
  # Step 2: Identify peaks above the threshold
  peaks <- data[abs_signal >= threshold]
  setorder(peaks, chromosome, start)  # Order by chromosome and start
  
  # Step 3: Create a group for peaks within the window based on the genomic distance using vectorized operation to group based on proximity within the window
  peaks[, group := cumsum((start > shift(end, type = "lag", fill = first(end)) + 100) | (chromosome != shift(chromosome, type = "lag", fill = first(chromosome))))]
  
  # Step 4: Summing peaks and signals within each group
  peaks <- peaks[, .(
    start = min(start),
    end = max(end),
    signal_sum = sum(abs(polymerase_signal)),
    peak_count = .N
  ), by = .(chromosome, group)]
  
  # Step 5: Filter out peaks where the signal sum is below the threshold
  peaks <- peaks[peak_count > 1 & signal_sum >= signal_sum_threshold]
  
  # Step 6: Clean up and return results
  peaks[, group := NULL]  # Remove temporary grouping column
  
  return(peaks)
}

#-----------

# Identify peaks in both strands
positive_peaks <- find_peaks(positive_strand, threshold = 1, signal_sum_threshold = 5)
negative_peaks <- find_peaks(negative_strand, threshold = 1, signal_sum_threshold = 5)

#-------------------------------------------------------------------

# Function to find divergent transcription patterns within a predefined organism dependent nucleotide window
find_divergent_transcription <- function(pos_peaks, neg_peaks, max_window) {
  # Step 1: Preprocessing
  neg_peaks[, end := start + as.numeric(max_window)]  # Adjust negative peaks to set window size
  
  # Step 2: Set keys on full datasets (once, not in the loop)
  setkey(pos_peaks, chromosome, start, end)
  setkey(neg_peaks, chromosome, start, end)
  
  # Step 3: Apply foverlaps for all chromosomes at once (no loop)
  joined <- foverlaps(neg_peaks, pos_peaks, by.x = c("chromosome", "start", "end"), 
                      by.y = c("chromosome", "start", "end"), nomatch = 0)
  
  # Step 4: Filter peaks to ensure they are within the window
  joined <- joined[abs(start - i.start) <= as.numeric(max_window)]
  
  # Step 5: Summing and processing peaks
  if (nrow(joined) > 0) {
    # Step 5a: Calculate region boundaries
    joined[, `:=`(region_start = pmin(i.start, start), region_end = pmax(i.end, end))]
    
    # Step 5b: Calculate total signal
    joined[, total_signal := abs(signal_sum) + abs(i.signal_sum)]
    
    # Step 6: Reduce to necessary columns and return results
    divergent_regions <- joined[, .(chromosome, start = region_start, end = region_end, total_signal)]
    
    setorder(divergent_regions, chromosome, start)
    return(divergent_regions)
  } else {
    return(data.table(chromosome = character(), start = integer(), end = integer(), total_signal = numeric()))
  }
}

#-----------

# Find regions of divergent transcription
divergent_transcription <- find_divergent_transcription(positive_peaks, negative_peaks, nt_window)

#-------------------------------------------------------------------

# Function to merge overlapping regions
merge_overlapping_regions <- function(regions) {
  # Step 1: Ensure regions are ordered by chromosome and start position
  setorder(regions, chromosome, start)
  
  # Step 2: Create a vector to track whether regions are overlapping
  # Overlap occurs if the next region starts before or at the current region's end
  regions[, overlap := (shift(end, type = "lag", fill = first(end)) >= start & shift(chromosome, type = "lag", fill = first(chromosome)) == chromosome)]
  
  # Step 3: Use cumulative sum of non-overlap events to create unique group IDs for overlapping regions
  regions[, overlap_group := cumsum(!overlap)]
  
  # Step 4: Merge overlapping regions within the same group
  merged_regions <- regions[, .(
    start = min(start),
    end = max(end),
    total_signal = sum(total_signal)
  ), by = .(chromosome, overlap_group)]
  
  # Step 5: Remove the temporary grouping column
  merged_regions[, overlap_group := NULL]
  
  return(merged_regions)
}

#-----------

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
