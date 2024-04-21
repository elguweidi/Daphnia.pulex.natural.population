# This R code was used to plot the neighbor-joining tree.
#==========================================================================#
library(tidyverse)
library(adegenet)
library(polysat)
library(ape)
library(TreeTools)
#==========================================================================#
dataIn = read.table('square-distance-matrix.txt',header = FALSE)
dataIn = as.matrix(dataIn); colnames(dataIn) =NULL
min(dataIn);max(dataIn)
dataIn[dataIn <0] = 0
dataIn = sqrt(dataIn)
min(dataIn);max(dataIn)
dataIn[1:5,1:5]
M <- matrix(NA, nrow(dataIn), ncol(dataIn))
M[lower.tri(M)] <- dataIn[lower.tri(dataIn)]
M[1:5,1:5]
min(M, na.rm=T);max(M, na.rm=T)
newNmaesData <- read.table('NamesOfSamples.txt', header = TRUE)
row.names(M) <- newNmaesData$sample
colnames(M) <- newNmaesData$sample
write.csv(M, file = 'lowerTriangularMatrix.csv',na="")
#==========================================================================#
tr <- nj(M)
annColor = c()
n=0
for (i in unique(newNmaesData$population)){
  n=n+1
  tempS = subset(newNmaesData, newNmaesData$population == i)
  tex1 = substr(unique(tempS$sample),1,nchar(i))
  annColor = c(annColor, rep(n,length(tex1)))
}
#==========================================================================#
myPal <- colorRampPalette(c("red","#40B0A6","#E66100","#5D3A9B","#1A85FF","black","#D35FB7","#FFC20A","springgreen","slategray"))
png("phylogram.png", width=8000, height=13000, res=350)
plot(tr, show.tip=FALSE, edge.width=3)
tiplabels("   ", bg=transp(num2col(annColor, col.pal=myPal),1.0), cex=2, fg="transparent")
add.scale.bar(0,-3,cex = 5, font = 1, col = "black",lwd = 10)
temp <- unique(newNmaesData$population, 10)
legend("topright", fill=transp(num2col(c(1:10), col.pal=myPal),1.0), leg=temp, ncol=2, cex=3.0)
dev.off()
#==========================================================================#
