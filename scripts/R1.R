#====================

args <- commandArgs()
organism <- args[6]
sample <- args[7]

path_to_genes <- paste0(organism, "/genome/genes.txt")

refGene = read.table(path_to_genes, header=T, sep="\t", quote="")

names(refGene) =c("chr", "txStart", "txEnd", "strand", "geneName")

refGene = subset(refGene, txEnd>=10499 & txStart>=10499)

refGene_pl = subset(refGene, strand=="+")
refGene_mn = subset(refGene, strand=="-")

### Genes on the plus strand:

refGene_pl$TSS = refGene_pl$txStart       # TSS
refGene_pl$CPS = refGene_pl$txEnd         # CPS

refGene_pl$DIVs = refGene_pl$txStart-750  # region of divergent transcription
refGene_pl$DIVe = refGene_pl$txStart-251

refGene_pl$PPs = refGene_pl$TSS-250       # promoter-proximal region
refGene_pl$PPe = refGene_pl$TSS+249

refGene_pl$GBs = refGene_pl$TSS+250       # genebody
refGene_pl$GBe = refGene_pl$CPS-501

refGene_pl$CPSs = refGene_pl$CPS-500      # CPS region
refGene_pl$CPSe = refGene_pl$CPS+499

refGene_pl$TWs = refGene_pl$CPS+500       # termination window
refGene_pl$TWe = refGene_pl$CPS+10499

#### Genes on the minus strand:

refGene_mn$TSS = refGene_mn$txEnd         # TSS
refGene_mn$CPS = refGene_mn$txStart       # CPS

refGene_mn$DIVs = refGene_mn$txEnd+251    # divergent transcription region
refGene_mn$DIVe = refGene_mn$txEnd+750

refGene_mn$PPs = refGene_mn$TSS-249       # promoter-proximal region
refGene_mn$PPe = refGene_mn$TSS+250

refGene_mn$GBs = refGene_mn$CPS+501       # genebody
refGene_mn$GBe = refGene_mn$TSS-250

refGene_mn$CPSs = refGene_mn$CPS-499      # CPS region
refGene_mn$CPSe = refGene_mn$CPS+500

refGene_mn$TWs = refGene_mn$CPS-10499     # termination window
refGene_mn$TWe = refGene_mn$CPS-500

#### combine the data of plus and minus stands:

refGene = rbind(refGene_pl, refGene_mn)

refGene$promC1 = refGene$TSS-500
refGene$promC2 = refGene$TSS+500

#### generate data files:

output_folder = paste0(organism, "/analysis/functionalGenomics_", sample)

write.table(refGene, file=paste0(output_folder, "/refGene_allTranscripts.txt"), col.names=F, row.names=F, quote=F, sep="\t")
write.table(refGene, file=paste0(output_folder, "/refGene_allTranscripts_withHeader.txt"), col.names=T, row.names=F, quote=F, sep="\t")
write.table(refGene[,c("chr","promC1","promC2","geneName")], file=paste0(output_folder, "/refGenes_TSSpm500.txt"), col.names=F, row.names=F, quote=F, sep="\t")

save.image()
q()
y
