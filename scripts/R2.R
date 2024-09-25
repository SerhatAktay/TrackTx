args <- commandArgs(trailingOnly = TRUE)
organism <- args[1]
sample <- args[2]
output_folder <- paste0(organism, "/analysis/functionalGenomics_", sample)

# Ensure output folder exists
if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
}

# File paths
active_genes_path <- paste0(output_folder, "/activeGenes.bed")
ref_gene_path <- paste0(output_folder, "/refGene_allTranscripts_withHeader.txt")

# Check if files exist
if (!file.exists(active_genes_path)) {
    stop(paste("Error: File", active_genes_path, "does not exist."))
}
if (!file.exists(ref_gene_path)) {
    stop(paste("Error: File", ref_gene_path, "does not exist."))
}

# Read input files
Active <- read.table(active_genes_path, sep="\t", quote="")
refGene <- read.table(ref_gene_path, header=T, sep="\t", quote="")

# Sanity check for activeGenes.bed column count
if (ncol(Active) < 4) {
    stop("Error: activeGenes.bed does not have the expected number of columns.")
}

# Subset active genes
refGeneAct <- subset(refGene, geneName %in% Active[,4])

# Write promoter-proximal and divergent transcription regions
write.table(refGeneAct[,c("chr","PPs","PPe","geneName","geneName","strand")], 
            file=paste0(output_folder, "/ppPolII.txt"), col.names=F, row.names=F, quote=F, sep="\t")
write.table(refGeneAct[,c("chr","DIVs","DIVe","geneName","geneName","strand")], 
            file=paste0(output_folder, "/divTx.txt"), col.names=F, row.names=F, quote=F, sep="\t")

# Handle short genes
if (any(refGeneAct$txEnd < refGeneAct$txStart)) {
    warning("Some genes have txEnd < txStart. These genes will be excluded.")
    refGeneAct <- subset(refGeneAct, txEnd >= txStart)
}

shortGenes <- subset(refGeneAct, txEnd - txStart <= 750)
refGeneAct_ <- subset(refGeneAct, txEnd - txStart > 750)

shortGenes$geneLength <- shortGenes$txEnd - shortGenes$txStart

# Write short genes and their full information
write.table(shortGenes[,c("chr","txStart","txEnd","geneName","geneName","strand")], 
            file=paste0(output_folder, "/shortGenes.txt"), col.names=F, row.names=F, quote=F, sep="\t")
write.table(shortGenes[,c("chr","txStart","txEnd","strand","geneName","TSS","CPS","DIVs","DIVe","PPs","PPe","GBs","GBe","CPSs","CPSe","TWs","TWe","promC1","promC2")], 
            file=paste0(output_folder, "/shortGenesAllInfo.txt"), col.names=T, row.names=F, quote=F, sep="\t")

# Write CPS and termination window regions
write.table(refGeneAct_[,c("chr","CPSs","CPSe","geneName","geneName","strand")], 
            file=paste0(output_folder, "/CPS.txt"), col.names=F, row.names=F, quote=F, sep="\t")
write.table(refGeneAct_[,c("chr","TWs","TWe","geneName","geneName","strand")], 
            file=paste0(output_folder, "/TW.txt"), col.names=F, row.names=F, quote=F, sep="\t")

# Ensure gene body length is positive
if (!all(c("GBs", "GBe") %in% names(refGeneAct_))) {
    stop("Error: GBs or GBe columns are missing.")
}
refGeneAct_ <- subset(refGeneAct_, GBe - GBs > 1)

# Write gene body regions
write.table(refGeneAct_[,c("chr","GBs","GBe","geneName","geneName","strand")], 
            file=paste0(output_folder, "/geneBody.txt"), col.names=F, row.names=F, quote=F, sep="\t")

# Save relevant objects instead of the entire environment
save(refGeneAct, refGeneAct_, shortGenes, file = paste0(output_folder, "/R2_workspace.RData"))

q("no")
