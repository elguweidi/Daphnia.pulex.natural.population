# this code is used to calculate rDNA copy number from read depth summaries.

library(tidyverse)
#-------------------------------------------------------------------#
df <- read_delim('Summary.txt',col_names = FALSE)
#-------------------------------------------------------------------#
df = df[,c(1,2,3)]
colnames(df) <- c("sample", "gene","mean_count")
#-------------------------------------------------------------------#
# See all regions included in file
unique(df$gene) 
#-------------------------------------------------------------------#
# Mean of exons read depth (Exons start with EFX)
df = df %>% pivot_wider(names_from = gene, values_from = mean_count)
df$mean_count = round(rowMeans(select(df,contains("EFX"))),3)
#-------------------------------------------------------------------#
# Calculate haploid CN
df$IGS1.1 = round(df$IGS1 / df$mean_count,0)
df$IGS2.1 = round(df$IGS2 / df$mean_count,0)
df$r18S.1 = round(df$r18S / df$mean_count,0)
df$r28S.1 = round(df$r28S / df$mean_count,0)
#-------------------------------------------------------------------#
# Calculate diploid CN
df$r18S.2 = df$r18S.1 * 2
df$r28S.2 = df$r28S.1 * 2
#-------------------------------------------------------------------#
# Order by sample ID
df <- df[order(df$sample),]
#-------------------------------------------------------------------#
# Save to file:
write_csv(df, file = "cn.csv")
#-------------------------------------------------------------------#

