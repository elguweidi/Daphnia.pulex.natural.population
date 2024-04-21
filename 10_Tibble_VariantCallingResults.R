# This R code was used to tabulate the variance calling results from the vcf files.
#=========================================================================#
library(tidyverse)
#=========================================================================#
# This step is to read the names and the number of all files to be analyzed 
listOfFiles <- list.files(path = "../2.VCF_Results/" , pattern = "\\.vcf$")
# Print number of files read 
cat ("Number of files read is" ,length(listOfFiles), "\n")
#=========================================================================#
# Create a container to store outputs of each loop.
LoopOutPutStore = c()
for (i in 1:length(listOfFiles)){
theSample = str_remove(listOfFiles[i], ".vcf")
theSample = str_remove(theSample, "Results_") 
cat ("Processing sample #", i, " >>> Sample ID is: ", theSample, "\n")
fileName = paste0("df.",theSample)
# Read data
fileName = read_delim(paste0("../2.VCF_Results/",listOfFiles[i]),col_names = TRUE)
# Select columns needed
fileName = fileName[!grepl("EFX", fileName$CHROM), ]
fileName <- fileName %>% 
        rename("gene" = "CHROM", "position" = "POS", "allele1" = "REF", "depth" = "DP")
# Add 0 to the left for loci that are less than 4 digits
fileName <- fileName %>%
mutate(
      position = as.character(position),
      position = str_pad(
        position,
        width = 4,
        side = "left",
        pad = "0"
      )
    )

# Add new columns for the sample and locus
fileName <- fileName %>%
    add_column(sample = theSample, .before = "gene") %>%
    unite('locus', gene:position, remove = FALSE, sep = "_")

# Split the alleles and their read depth in separate columns
fileName <- fileName %>%
    separate(ALT, c("allele2", "allele3", "allele4")) %>%
    mutate(across(starts_with("allele"), ~ replace_na(., "none"))) %>% 
    separate(AD, c("blank", "a1_count", "a2_count", "a3_count", "a4_count")) %>% 
    select(-blank) %>% 
    mutate(across(ends_with("count"), ~ replace_na(., "0"))) %>%
    mutate(
      a1_count = parse_number(a1_count),
      a2_count = parse_number(a2_count),
      a3_count = parse_number(a3_count),
      a4_count = parse_number(a4_count)
    )
LoopOutPutStore = rbind(LoopOutPutStore,fileName)
}

# Rename the samples:
NewNames = as.data.frame(sort(unique(LoopOutPutStore$sample)))
colnames(NewNames) = 'sample'
library(stringr)
NewNames$NewNames = paste0("CH", str_pad(c(1:nrow(NewNames)), 2, pad = "0"))
write_csv(NewNames, file= "SamplesOldNewNames.csv")

LoopOutPutStore = merge(LoopOutPutStore,NewNames, by = 'sample')
LoopOutPutStore$sample = LoopOutPutStore$NewNames
LoopOutPutStore = subset(LoopOutPutStore, select = -NewNames )

write_csv(LoopOutPutStore, file= "Combined_VC.csv")
#=====================================================================================#






