# This R code was used to run the haplotype analysis using the R function haplosep.
#===============================================================================#
library(devtools)
library(tidyverse)
library(haploSep)
library(adegenet)
library(polysat)
# Read prunned data
pruData <- read_csv("PrunedData.csv")
a1 <- select(pruData, contains(".1"))
m.new = floor(log2(ncol(a1)))
m.new
#======================================================#
# Sort by sample:
pruData <- pruData[order(pruData$sample),]
#======================================================#
data.a1=select(pruData,contains(".1"))
Alleles = colnames(data.a1)
data.a1 = cbind(sample=pruData$sample, data.a1)
#===============================================================================#
data.a1 <- as.data.frame(data.a1)
sampleNames = data.a1[,1]
data.a1 = data.a1[,2:ncol(data.a1)]
colnames(data.a1)
names(data.a1) <- NULL
matrix.a1 <- as.matrix(data.a1)
matrix.a1.t <- t(matrix.a1)
results.a1.new <- haploSep(data = matrix.a1.t, nHaplo=m.new, stabEval=FALSE, bias=FALSE)
haploFrq.a1.new <- as.data.frame(results.a1.new$haploFrq)
haploStr.a1.new <- as.data.frame(results.a1.new$haploStr)
write_csv(haploFrq.a1.new, file=paste0("haploFrq.a1.csv"))
write_csv(haploStr.a1.new, file=paste0("haploStr.a1.csv"))
#==========================================================================#
dat <- read_csv("haploStr.a1.csv")
dat$Sum = rowSums(dat)
dat$locus = Alleles
gcta=read_csv("../AllFilesCombinedCorrectClean.csv")
gcta=gcta[,c('locus','allele1','allele2')]
gcta = gcta[!duplicated(gcta$locus),]
dat2 = merge(dat,gcta)
dat2 <- subset(dat2, select=c(locus,allele1,allele2,Haplo1:Sum))
dat3=dat2
dat3T <- data.frame(t(dat3))
dat3T$locus = row.names(dat3T)
write_csv(dat3, file=paste0("HaplotypeGenotype.csv"))
write_csv(dat3T, file=paste0("HaplotypeGenotypeT.csv"), col_names = F)
#=================================================================
# Read the freqs of the 10 haplotypes:
hap10=read_csv("haploFrq.a1.csv")
hap10=t(hap10)
hap10=data.frame(hap10)
hapNmaes = c('Hap.1','Hap.2','Hap.3','Hap.4','Hap.5','Hap.6','Hap.7','Hap.8','Hap.9','Hap.10')
colnames(hap10) = hapNmaes[1:ncol(hap10)]
hap10$sample = sampleNames
hap10 <- subset(hap10, select=c(sample,1:(ncol(hap10)-1)))
#==========================================================================================================#
# Read copy number data
cn.data  <- read_csv("../../cn.csv")
# Pick up columns sample and r28S.2
cn.data <- cn.data[,c("sample","r28S.2")]
# Merge the two files
cnPruned <- merge(cn.data,hap10, by = 'sample')
#----------------------------------------
# Drop population column and rename copy number and sample columns
# Sample >>> Population & r28S.2 to Genomes
df.data3 <- cnPruned %>%
  rename(Genomes=r28S.2, Population=sample)
# Set any NAs to 0
df.data3[is.na(df.data3)] <- 0
# Make the population as a row 
df.popgen2.row <-df.data3 %>% column_to_rownames(var = "Population")
mode(df.popgen2.row)
write.table(df.popgen2.row, file = paste0("hapFile.txt"), quote = F, sep = " ",
            row.names = T,col.names = TRUE)
# Run polysat
df.popgen2.count <- polysat::freq.to.genpop(df.popgen2.row)
# Save data
save(df.popgen2.count, file=paste0("df.popgen2.count.Rda"))
#----------------------------------------------------------------------#
# Run genpop
df.genpop2 <-adegenet::genpop(df.popgen2.count, ploidy = '1', type = "codom")
save(df.genpop2, file = paste0("df.genpop2.Rda"))
df.genpop2
df.genpop2[1:nrow(data.a1)]@tab
df.genpop2[1:nrow(data.a1)]@loc.fac
df.genpop2[1:nrow(data.a1)]@loc.n.all
df.genpop2[1:nrow(data.a1)]@all.names

nPop(df.genpop2)
nLoc(df.genpop2)
tab(df.genpop2)
allele.number <- nAll(df.genpop2)
loc.names <- locNames(df.genpop2)
alleles <- alleles(df.genpop2)

df= df.genpop2[1:nrow(data.a1)]@tab
write.table(df, file = "df.genpop2HapTab.txt", quote = F, sep = " ",
            row.names = T,col.names = TRUE)