args <- commandArgs(trailingOnly = TRUE)
organism <- args[1]
sample <- args[2]
nt_window <- as.numeric(args[3])

# Define paths
path_to_bedgraphs <- paste0(organism, "/bigWig/")
output_folder <- paste0(organism, "/analysis/functionalGenomics_", sample)

# Ensure output folder exists
if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
}

# File paths for positive and negative strands
pos_file <- paste0(path_to_bedgraphs, sample, "_allMap_unnorm_pl.bedgraph")
neg_file <- paste0(path_to_bedgraphs, sample, "_allMap_unnorm_mn.bedgraph")

# Check if files exist
if (!file.exists(pos_file)) {
    stop(paste("Error: File", pos_file, "does not exist."))
}
if (!file.exists(neg_file)) {
    stop(paste("Error: File", neg_file, "does not exist."))
}

# Set the CRAN mirror to avoid the "no mirror" error
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Ensure a personal library path is set (if it doesn't exist)
if (Sys.getenv("R_LIBS_USER") == "") {
  Sys.setenv(R_LIBS_USER = "~/R/x86_64-pc-linux-gnu-library/4.0")  # Set your personal library path
}

# Create the personal library directory if it doesn't exist
if (!dir.exists(Sys.getenv("R_LIBS_USER"))) {
  dir.create(Sys.getenv("R_LIBS_USER"), recursive = TRUE)
}

required_packages <- c("data.table", "dplyr")

# Helper function to install and load packages silently
install_and_load <- function(pkg) {
  # Capture output and suppress warnings/messages, but not errors
  suppressMessages(suppressWarnings(capture.output({
    if (!require(pkg, character.only = TRUE)) {
      install.packages(pkg, dependencies = TRUE, lib = Sys.getenv("R_LIBS_USER"), ask = FALSE)
      library(pkg, character.only = TRUE, lib.loc = Sys.getenv("R_LIBS_USER"))
    }
  })))
}

# Check and install packages if not available
lapply(required_packages, install_and_load)

# Read the input files
positive_strand.data <- read.table(pos_file, header = FALSE, col.names = c("chromosome", "start", "end", "polymerase_signal"))
negative_strand.data <- read.table(neg_file, header = FALSE, col.names = c("chromosome", "start", "end", "polymerase_signal"))

# Convert to data.table and sort
positive_strand <- setorder(setDT(positive_strand.data))
negative_strand <- setorder(setDT(negative_strand.data))

# Peak finding function
find_peaks <- function(data, threshold, signal_sum_threshold) {
    data[, abs_signal := abs(polymerase_signal)]
    peaks <- data[abs_signal >= threshold]
    setorder(peaks, chromosome, start)
    
    peaks[, group := cumsum((start > shift(end, type = "lag", fill = first(end)) + 100) | (chromosome != shift(chromosome, type = "lag", fill = first(chromosome))))]
    
    peaks <- peaks[, .(
        start = min(start),
        end = max(end),
        signal_sum = sum(abs(polymerase_signal)),
        peak_count = .N
    ), by = .(chromosome, group)]
    
    peaks <- peaks[peak_count > 1 & signal_sum >= signal_sum_threshold]
    peaks[, group := NULL]
    return(peaks)
}

# Identify peaks
positive_peaks <- find_peaks(positive_strand, threshold = 1, signal_sum_threshold = 5)
negative_peaks <- find_peaks(negative_strand, threshold = 1, signal_sum_threshold = 5)

# Find divergent transcription regions
find_divergent_transcription <- function(pos_peaks, neg_peaks, max_window) {
    neg_peaks[, end := start + as.numeric(max_window)]
    
    setkey(pos_peaks, chromosome, start, end)
    setkey(neg_peaks, chromosome, start, end)
    
    joined <- foverlaps(neg_peaks, pos_peaks, by.x = c("chromosome", "start", "end"), by.y = c("chromosome", "start", "end"), nomatch = 0)
    joined <- joined[abs(start - i.start) <= as.numeric(max_window)]
    
    if (nrow(joined) > 0) {
        joined[, `:=`(region_start = pmin(i.start, start), region_end = pmax(i.end, end))]
        joined[, total_signal := abs(signal_sum) + abs(i.signal_sum)]
        
        divergent_regions <- joined[, .(chromosome, start = region_start, end = region_end, total_signal)]
        setorder(divergent_regions, chromosome, start)
        
        return(divergent_regions)
    } else {
        return(data.table(chromosome = character(), start = integer(), end = integer(), total_signal = numeric()))
    }
}

divergent_transcription <- find_divergent_transcription(positive_peaks, negative_peaks, nt_window)

# Merging overlapping regions
merge_overlapping_regions <- function(regions) {
    setorder(regions, chromosome, start)
    
    regions[, overlap := (shift(end, type = "lag", fill = first(end)) >= start & shift(chromosome, type = "lag", fill = first(chromosome)) == chromosome)]
    
    regions[, overlap_group := cumsum(!overlap)]
    
    merged_regions <- regions[, .(
        start = min(start),
        end = max(end),
        total_signal = sum(total_signal)
    ), by = .(chromosome, overlap_group)]
    
    merged_regions[, overlap_group := NULL]
    
    return(merged_regions)
}

# Merge the overlapping regions
merged_divergent_transcription <- merge_overlapping_regions(divergent_transcription)

# Write the results to the BED file
write.table(merged_divergent_transcription, file = paste0(output_folder, "/divergent_transcription.bed"), quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")

# Save relevant objects
save(positive_peaks, negative_peaks, merged_divergent_transcription, file = paste0(output_folder, "/R6_workspace.RData"))

q("no")
