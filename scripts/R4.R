args <- commandArgs()
organism <- args[6]
sample <- args[7]

library(dplyr)
library(gtools)

col_names <- c("chr", "start", "end", "score", "gene", "strand", "start_2", "end_2", "colour")

file <- paste0(organism, "/analysis/functionalGenomicRegions_", sample, ".bed")

# Read the .bed file into a data frame with specified column names
df <- read.delim2(file, col.names=col_names, sep="\t", header = FALSE, skip = 1)

# Remove columns "start_2", end_2", and "colour"
df <- df[, !colnames(df) %in% c("start_2", "end_2", "colour")]

# Split the dataframe based on the "strand" column
df_minus <- df[df$strand == "-", ]
df_plus <- df[df$strand == "+", ]

# Assuming 'your_dataframe' is your dataframe
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
