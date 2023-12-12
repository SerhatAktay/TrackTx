# Retrieve command-line arguments in R
args <- commandArgs(trailingOnly = TRUE)
organism <- args[1]
number_of_samples <- as.numeric(args[2])
samples_to_extract <- (3):(2 + number_of_samples)
sample_list <- args[samples_to_extract]
conditions_to_extract <- (3 + number_of_samples):(2 + (2 * number_of_samples))
sample_conditions <- args[conditions_to_extract]

suppressWarnings()

###########

## Required packages
packages = c("DESeq2", "dplyr", "ggrepel", "ggplot2")

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

###########

data_frames_list <- setNames(
  lapply(sample_list, function(sample_file) {
    full_path <- paste0(organism, "/analysis/differential_expression/DE_", sample_file, ".bed")
    df <- read.table(full_path, header = TRUE, col.names = c("chr", "start", "end", "strand", "gene", "signal"))
    df <- subset(df, select = c("gene", "signal"))
    return(df)
  }),
  sample_list
)

for (sample_name in sample_list) {
  current_df <- data_frames_list[[sample_name]]
  current_df <- current_df %>%
    group_by(gene) %>%
    summarise(signal = sum(signal))
}

for (sample_name in names(data_frames_list)) {
  colnames(data_frames_list[[sample_name]]) <- c('gene', sample_name) 
}

# Merge multiple data frames using Reduce
merged_df <- Reduce(function(x, y) merge(x, y, by = "gene", all = TRUE), data_frames_list)
merged_df <- merged_df %>%
  distinct(gene, .keep_all = TRUE)

coldata <- data.frame(
  sampleID = sample_list,
  condition = sample_conditions
)
row.names(coldata) <- coldata$sampleID

count_matrix <- merged_df
row.names(count_matrix) <- count_matrix$gene
count_matrix <- count_matrix[, -1]
count_matrix <- na.omit(count_matrix)

# Create a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = coldata,
                              design = ~ condition)

# Run DESeq2 analysis
dds <- DESeq(dds)

# Get differential expression results
results <- results(dds)

# Extract top up- and downregulated genes
top_upregulated <- head(results[order(-results$log2FoldChange), ], 50)
top_downregulated <- head(results[order(results$log2FoldChange), ], 50)

# Convert the DEseq-objects to DataFrames
upregulated_df <- data.frame(top_upregulated)
downregulated_df <- data.frame(top_downregulated)
result_df <- data.frame(results)
rows <- c(row.names(upregulated_df), row.names(downregulated_df))

# Merge and save dataFrames of DE-genes
merged_regulated_genes <- rbind(upregulated_df, downregulated_df)
genes_file <- paste0(organism, "/analysis/differential_expression/top_regulated_genes.txt")
write.table(merged_regulated_genes, file = genes_file, quote = FALSE, sep = "\t")

# Create a volcano plot
ggplot(result_df, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(color = ifelse(row.names(result_df) %in% row.names(upregulated_df), "Upregulated", 
                                ifelse(row.names(result_df) %in% row.names(downregulated_df), "Downregulated", "NotSig")))) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "NotSig" = "gray")) +
  geom_text_repel(data = bind_rows(upregulated_df, downregulated_df),
                  aes(label = rows), max.overlaps = 20, nudge_x = 0.3, nudge_y = 0.3, size = 2) +
  labs(title = "Volcano Plot",
       x = "log2 Fold Change",
       y = "-log10(p-value)",
       color = "Significance") +
  theme_minimal()

plot_file <- paste0(organism, "/analysis/differential_expression/volcano_plot.png")
ggsave(plot_file, bg = "white")



