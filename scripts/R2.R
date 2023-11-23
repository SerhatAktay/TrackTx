#====================

args <- commandArgs()
organism <- args[6]
sample <- args[7]
output_folder = paste0(organism, "/analysis/functionalGenomics_", sample)

Active = read.table(paste0(output_folder, "/activeGenes.bed"), sep="\t", quote="")
refGene = read.table(paste0(output_folder, "/refGene_allTranscripts_withHeader.txt"), header=T, sep="\t", quote="")

refGeneAct = subset(refGene, geneName %in% Active[,4])

write.table(refGeneAct[,c("chr","PPs","PPe","geneName","geneName","strand")], file=paste0(output_folder, "/ppPolII.txt"), col.names=F, row.names=F, quote=F, sep="\t")
write.table(refGeneAct[,c("chr","DIVs","DIVe","geneName","geneName","strand")], file=paste0(output_folder, "/divTx.txt"), col.names=F, row.names=F, quote=F, sep="\t")

shortGenes = subset(refGeneAct, txEnd-txStart<=750)
refGeneAct_ = subset(refGeneAct, txEnd-txStart>750)
shortGenes$geneLength <- shortGenes[,3] - shortGenes[,2]

write.table(shortGenes[,c("chr","txStart","txEnd","geneName","geneName","strand")], file=paste0(output_folder, "/shortGenes.txt"), col.names=F, row.names=F, quote=F, sep="\t")
write.table(shortGenes[,c("chr","txStart","txEnd","strand","geneName","TSS","CPS","DIVs","DIVe","PPs","PPe","GBs","GBe","CPSs","CPSe","TWs","TWe","promC1","promC2")], 
						file=paste0(output_folder, "/shortGenesAllInfo.txt"), col.names=T, row.names=F, quote=F, sep="\t")

write.table(refGeneAct_[,c("chr","CPSs","CPSe","geneName","geneName","strand")], file=paste0(output_folder, "/CPS.txt"), col.names=F, row.names=F, quote=F, sep="\t")
write.table(refGeneAct_[,c("chr","TWs","TWe","geneName","geneName","strand")], file=paste0(output_folder, "/TW.txt"), col.names=F, row.names=F, quote=F, sep="\t")

refGeneAct_ = subset(refGeneAct_,GBe-GBs>1)  #ensuring no negative gene body lengths remain.
write.table(refGeneAct_[,c("chr","GBs","GBe","geneName","geneName","strand")], file=paste0(output_folder, "/geneBody.txt"), col.names=F, row.names=F, quote=F, sep="\t")

save.image()
q()
y

