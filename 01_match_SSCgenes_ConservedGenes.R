# This R code was used to  match the Daphnia pulex single copy 
# genes with the conserved genes reported by Colbourne et al., (2011)
# Colbourne JK, et al.2011. The ecoresponsive genome of Daphnia pulex. Science 331:555-561.
id <- read.table("SingleCopyGenesFromEnsemble.txt", header =T)
daphGene <- read.table("ConservedGenes.txt", header = T)
head(id); head(daphGene)
singleCopyConservedGenes <- subset(id, (Gene_stable_ID %in% daphGene$Gene_stable_ID))
write.table(singleCopyConservedGenes, file = "listOfsingleCopyConservedGenes.txt", quote = FALSE, sep = " ", row.names = FALSE, col.names = TRUE)
