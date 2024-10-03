#-----------------------------------------

args <- commandArgs()
organism <- args[6]
sample <- args[7]

#-----------------------------------------

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

required_packages <- c("gtools", "dplyr")

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

#-----------------------------------------

col_names <- c("chr", "start", "end", "gene", "score", "strand", "start_2", "end_2", "color")

file <- paste0(organism, "/analysis/functionalGenomicRegions_", sample, ".bed")

# Read the .bed file into a data frame with specified column names
df <- read.delim2(file, col.names=col_names, sep="\t", header = FALSE, skip = 1)

# Remove columns "start_2", end_2", and "colour"
df <- df[, !colnames(df) %in% c("start_2", "end_2", "color")]

# Split the dataframe based on the "strand" column
df_minus <- df[df$strand == "-", ]
df_plus <- df[df$strand == "+", ]

# Collapse genes into one line per row per strand
df_minus <- df_minus %>%
  group_by(gene) %>%
  summarise(
    chr = first(chr),   # Keep the first chr value within each group
    strand = first(strand),  # Keep the first strand value within each group
    start = min(start),
    end = max(end),
    score_sum = sum(score)
    # Add other columns you want to summarize
  ) %>%
  ungroup()

df_plus <- df_plus %>%
  group_by(gene) %>%
  summarise(
    chr = first(chr),   # Keep the first chr value within each group
    strand = first(strand),  # Keep the first strand value within each group
    start = max(start),
    end = min(end),
    score_sum = sum(score)
    # Add other columns you want to summarize
  ) %>%
  ungroup()

#-----------------------------------------

data <- rbind(df_minus, df_plus)

data <- data[mixedorder(data$chr), ]

# Assuming 'your_dataframe' is your dataframe and 'column1' is the column you want to split by
unique_values <- unique(data$chr)

# Initialize an empty list to store individual dataframes
df_list <- list()

# Split the dataframe, sort each subset, and store them in the list
for (value in unique_values) {
  subset_df <- data[data$chr == value, ]
  sorted_subset <- subset_df[order(subset_df$start), ]
  df_list[[as.character(value)]] <- sorted_subset
}

# Concatenate the sorted dataframes back together and rearrange the columns
data <- do.call(rbind, df_list)
data <- data[, c("chr", "start", "end", "strand", "gene", "score_sum")]

out_path <- paste0(organism, "/analysis/differential_expression/DE_", sample, ".bed")

# save the dataframe as a .bed file
write.table(data, out_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
