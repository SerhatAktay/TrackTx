# Optimized R1.R script

args <- commandArgs(trailingOnly = TRUE)
organism <- args[1]
sample <- args[2]

path_to_genes <- paste0("genes/genes_", organism, ".txt")

# Check if genes file exists
if (!file.exists(path_to_genes)) {
    stop(paste("Error: File", path_to_genes, "does not exist."))
}

# Read the reference gene data
refGene <- read.table(path_to_genes, header = TRUE, sep = "\t", quote = "")
names(refGene) <- c("chr", "txStart", "txEnd", "strand", "geneName")

# Filter out unwanted genes based on txStart and txEnd
refGene <- subset(refGene, txEnd >= 10499 & txStart >= 10499)

# Split into plus and minus strands
refGene_pl <- subset(refGene, strand == "+")
refGene_mn <- subset(refGene, strand == "-")

# Function to calculate regions based on strand
calculate_regions <- function(df, strand_type) {
    if (strand_type == "+") {
        df$TSS <- df$txStart
        df$CPS <- df$txEnd
        df$DIVs <- df$txStart - 750
        df$DIVe <- df$txStart - 251
        df$PPs <- df$TSS - 250
        df$PPe <- df$TSS + 249
        df$GBs <- df$TSS + 250
        df$GBe <- df$CPS - 501
        df$CPSs <- df$CPS - 500
        df$CPSe <- df$CPS + 499
        df$TWs <- df$CPS + 500
        df$TWe <- df$CPS + 10499
    } else if (strand_type == "-") {
        df$TSS <- df$txEnd
        df$CPS <- df$txStart
        df$DIVs <- df$txEnd + 251
        df$DIVe <- df$txEnd + 750
        df$PPs <- df$TSS - 249
        df$PPe <- df$TSS + 250
        df$GBs <- df$CPS + 501
        df$GBe <- df$TSS - 250
        df$CPSs <- df$CPS - 499
        df$CPSe <- df$CPS + 500
        df$TWs <- df$CPS - 10499
        df$TWe <- df$CPS - 500
    }
    return(df)
}

# Apply region calculations to both strands
refGene_pl <- calculate_regions(refGene_pl, "+")
refGene_mn <- calculate_regions(refGene_mn, "-")

# Combine the data for both strands
refGene <- rbind(refGene_pl, refGene_mn)

# Add promoter coordinates
refGene$promC1 <- refGene$TSS - 500
refGene$promC2 <- refGene$TSS + 500

# Define output folder and ensure it exists
output_folder <- paste0(organism, "/analysis/functionalGenomics_", sample)
if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
}

# Save the outputs
write.table(refGene, file = paste0(output_folder, "/refGene_allTranscripts.txt"), col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(refGene, file = paste0(output_folder, "/refGene_allTranscripts_withHeader.txt"), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(refGene[, c("chr", "promC1", "promC2", "geneName")], file = paste0(output_folder, "/refGenes_TSSpm500.txt"), col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

q("no")
