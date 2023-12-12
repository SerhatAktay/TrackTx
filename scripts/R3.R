#-----------------------------------------

# R-script for comparing the gene expression

args <- commandArgs()
organism <- args[6]
sample_1 <- args[7]
sample_2 <- args[8]

path_1 <- paste0(organism, "/analysis/functionalGenomicRegions_", sample_1, "_norm.bed")
path_2 <- paste0(organism, "/analysis/functionalGenomicRegions_", sample_2, "_norm.bed")

#-----------------------------------------

## Required packages
packages = c("data.table")

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

#-----------------------------------------

# Read in the two .bed files as data.tables
file1 <- fread(path_1, sep = "\t", header=FALSE)
file2 <- fread(path_2, sep = "\t", header=FALSE)

# Identify the rows where the fourth column differs between the two files and neither value is 0
diff_rows <- merge(file1, file2, by = c("V1","V2","V3","V5","V6","V7","V8","V9"))
diff_rows <- diff_rows[diff_rows$V4.x != diff_rows$V4.y,]
diff_rows <- diff_rows[diff_rows$V4.x != 0]
diff_rows <- diff_rows[diff_rows$V4.y != 0]

# Add a new column with the difference between the two values in column 4 and the relative difference between the two
diff_rows$diff <- diff_rows$V4.y - diff_rows$V4.x
diff_rows$V4.x[diff_rows$V4.x == 0] <- 0.1
diff_rows$rel_diff <- round(diff_rows$diff/diff_rows$V4.x, digits = 4)

# Select the desired columns and rename them
diff_rows <- diff_rows[, c("V1", "V2", "V3", "V5", "V6", "V7", "V8", "V9", "V4.x", "V4.y", "diff", "rel_diff")]
colnames(diff_rows) <- c("chr", "start", "end", "gene", "strand", "start", "end", "colour", "control", "heat shock", "difference", "rel_diff")

# Write the result to a new file
output_results <- paste0(organism, "/analysis/result_", organism, ".txt")
write.table(diff_rows, output_results, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

# Sort the genes into two files, up- and downregulated, and sort from most to least.
top10pl <- diff_rows[diff_rows$rel_diff >= 0]
top10mn <- diff_rows[diff_rows$rel_diff <= 0]

top10pl <- top10pl[order(-rel_diff)]
top10mn <- top10mn[order(rel_diff)]

top10pl = top10pl[!duplicated(top10pl$gene),]
top10mn = top10mn[!duplicated(top10mn$gene),]

top10pl <- top10pl[, c("chr", "start", "end", "gene", "control", "heat shock", "difference", "rel_diff")]
top10mn <- top10mn[, c("chr", "start", "end", "gene", "control", "heat shock", "difference", "rel_diff")]

# Write the result to a new file
output_top <- paste0(organism, "/analysis/topUpRegulated_", organism, ".txt")
output_min <- paste0(organism, "/analysis/topDownRegulated_", organism, ".txt")
write.table(top10pl, output_top, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(top10mn, output_min, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

print('Number of differently expressed genes: ')
print(nrow(diff_rows[!duplicated(diff_rows$gene),]))

print('Number of upregulated genes: ')
print(nrow(top10pl))

print('Number of downregulated genes: ')
print(nrow(top10mn))

print("Comparison of gene expression done!")

#-----------------------------------------
