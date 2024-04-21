# This are code was used to clean and prune the variance calling data.
#========================================================================#
library(tidyverse)
library("stringr")
abir.data <- read.csv ('AllFilesCombinedCorrectClean.csv')

# Get allele fregs
abir.data$a1_freq = abir.data$a1_count/(abir.data$a1_count+abir.data$a2_count+abir.data$a3_count+abir.data$a4_count)
abir.data$a2_freq = abir.data$a2_count/(abir.data$a1_count+abir.data$a2_count+abir.data$a3_count+abir.data$a4_count)
abir.data$a3_freq = abir.data$a3_count/(abir.data$a1_count+abir.data$a2_count+abir.data$a3_count+abir.data$a4_count)
abir.data$a4_freq = abir.data$a4_count/(abir.data$a1_count+abir.data$a2_count+abir.data$a3_count+abir.data$a4_count)

# Check if allele freqs sum up to 1:
sumFreqs = abir.data$a1_freq+abir.data$a2_freq+abir.data$a3_freq+abir.data$a4_freq
max(sumFreqs); min(sumFreqs)

# Subset the data to pivot it so we can do the data pruning
abir.data = abir.data[,c("sample", "locus", "gene", "position", "a1_freq", "a2_freq", "a3_freq", "a4_freq")]

abir.data.wide <- abir.data %>%
  select(sample, locus, a1_freq, a2_freq, a3_freq, a4_freq) %>%
  pivot_wider(names_from = locus,
              values_from = c(a1_freq, a2_freq, a3_freq, a4_freq), names_sep = "-")
#-----------------------------------------------------------------------------#
w1=select(abir.data.wide,contains("a1_freq"))
w2=select(abir.data.wide,contains("a2_freq"))
w3=select(abir.data.wide,contains("a3_freq"))
w4=select(abir.data.wide,contains("a4_freq"))

w1[is.na(w1)] <- 1
w2[is.na(w2)] <- 0
w3[is.na(w3)] <- 0
w4[is.na(w4)] <- 0

m=abir.data.wide[,1]
mm1=data.frame("Sum");names(mm1) <- c("sample")
mm2=data.frame("Mean");names(mm2) <- c("sample")

prePrun=rbind(m,mm1)
prePrun=rbind(prePrun,mm2)

ww1 = rbind(w1, round(colSums (w1),3))
ww1 = rbind(ww1,round(colMeans(w1),3))
ww2 = rbind(w2, round(colSums (w2),3))
ww2 = rbind(ww2,round(colMeans(w2),3))
ww3 = rbind(w3, round(colSums (w3),3))
ww3 = rbind(ww3,round(colMeans(w3),3))
ww4 = rbind(w4, round(colSums (w4),3))
ww4 = rbind(ww4,round(colMeans(w4),3))

prePrun=cbind(prePrun,ww1)
prePrun=cbind(prePrun,ww2)
prePrun=cbind(prePrun,ww3)
prePrun=cbind(prePrun,ww4)

sample = prePrun[1:(nrow(prePrun)),1:1]

freqs = prePrun[1:((nrow(prePrun))-2),2:(ncol(prePrun))]
freqs[colSums(freqs)==0] <- NULL

t=freqs
t = select(freqs,contains('a1'))
t1=t
t[(round(colMeans (t),3) <= 0.990) & (round(colMeans (t),3) >= 0.01) ] <- NULL

cat(nrow(t1),ncol(t1),ncol(t),ncol(t1)-ncol(t),floor(log2(ncol(t1)-ncol(t))),   '\n')
cat('N_Sample','N_Loci','N_removed_Loci','N_KeptLoci','N_Haps\n', sep=",", file = 'Number_SNP_Hap_Summary.csv')  
cat(nrow(t1),ncol(t1),ncol(t),ncol(t1)-ncol(t),floor(log2(ncol(t1)-ncol(t))),sep=",", file = 'Number_SNP_Hap_Summary.csv',append = TRUE)

# Now remove them
freqs[names(freqs) %in% names(t)] <- NULL

# Replace allele a1 with a2 for the same list of alleles
names(t) <- str_replace(names(t), "a1", "a2")
freqs[names(freqs) %in% names(t)] <- NULL
# Replace allele a2 with a3 for the same list of alleles
names(t) <- str_replace(names(t), "a2", "a3")
freqs[names(freqs) %in% names(t)] <- NULL
# Replace allele a3 with a4 for the same list of alleles
names(t) <- str_replace(names(t), "a3", "a4")
freqs[names(freqs) %in% names(t)] <- NULL

sums = round(colSums (freqs),3)
means = round(colMeans(freqs),3)

freqs = rbind(freqs, sums)
freqs = rbind(freqs,means)

# Join the data back with sample IDs
newData = cbind(sample, freqs)
#==================================================================================================#
# Drop last two columns (sum and mean)
pruData <- newData[1:(nrow(newData)-2),]
# Rename alleles
out= c()
out <- c(out, "sample")
dfColNames <- colnames(pruData)
# Rename D_obt_U1 and D_obt_U2 to IGS1 and IGS2, rsp.
dfColNames <- gsub("D_obt_U1", "IGS1", dfColNames)
dfColNames <- gsub("D_obt_U2", "IGS2", dfColNames)
dfColNames <- gsub("D_obt_U3", "IGS3", dfColNames)
dfColNames <- gsub("Dobtusa-C11-tpase", "tpase", dfColNames)
dfColNames <- gsub("Dobtusa-C11-3end","end3", dfColNames)

for (i in 2: length(dfColNames)){
  out <- c(out, paste0(substr(dfColNames[i],9,30),".",substr(dfColNames[i],2,2)))
}

colnames(pruData) = out

#======================================================#
# Check if freqs sum to 1 after data pruning:
nof = 'CheckFreqSumTo_1.txt' # Name of out file
if (file.exists(nof)) {
  file.remove(nof)
}

# Drop last two columns (sum and mean)
#newDataCheck <- newData[1:(nrow(newData)-2),]
colNamesOfdf = colnames(pruData)
colNamesOfdf = colNamesOfdf[2:(length(colNamesOfdf))]
colNamesOfdf = gsub("\\..*","",colNamesOfdf)
colNamesOfdf = unique(colNamesOfdf)
for (i in colNamesOfdf){
  theColToCheck = select(pruData,contains(i))
  cat (i,rowSums(theColToCheck),'\n',file = nof,append=TRUE)
}
#======================================================#
for (i in colNamesOfdf){
  theColToCheck = select(pruData,contains(i))
  cat (i,rowSums(theColToCheck),'\n',file = nof,append=TRUE)
}

# Check if freqs sum to 1 after data pruning:
nof2 = 'CheckAlleleFreq.txt' # Name of out file
if (file.exists(nof2)) {
  file.remove(nof2)
}

k=0
for (i in colNamesOfdf){
  k=k+1
  aa= select(pruData,contains(i))
  if (ncol(aa) > 1){
    aa= aa[ , order(names(aa))]}
  meansOut=round(as.numeric(colMeans(aa)),3)
  if (ncol(aa) == 1){ meansOut=c(meansOut, 'none','none','none')}
  if (ncol(aa) == 2){ meansOut=c(meansOut, 'none','none')}
  if (ncol(aa) == 3){ meansOut=c(meansOut, 'none')}
  cat(k,gsub("\\..*","",colnames(aa)[1]),meansOut,":",sub(".*\\.", "", colnames(aa)),'\n',
      file = nof2,append=TRUE)
  
}

LociToDel = c('IGS1_0443','r28S_0628','r28S_0639')
for (i in LociToDel){
  pruData = select(pruData,-contains(i))
}

#======================================================#
# Sort by sample:
pruData <- pruData[order(pruData$sample),]
#======================================================#
write.csv(pruData,"PrunedData.csv", row.names = FALSE,quote = FALSE)


