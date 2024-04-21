# This R code was used to correct third and fourth alleles if they are filliped. 
# Read data in
dataIn <- read.csv(file ="Combined_VC.csv")

# Unique loci available in the data
a <- unique(dataIn$locus)

# Create a container to store
DataStore <- c()
# Start looping over loci
for (i in 1:length(a)){
  b <- subset(dataIn, dataIn$locus == a[i])
  
  # Pass in loci that do not need modification
  if (nrow(b) < 2){
    b$New_allele3 <- b$allele3
    b$New_allele4 <- b$allele4
    DataStore <- rbind(DataStore,b)
    next
  }
  
  # Determine which allele is the most popular
  # and use it to be the most allele
  # to appear  
  e1 <- data.frame(table(b$allele3))
  e2 <- data.frame(table(b$allele4))
  
  k1= subset(e1, e1$Var1 == "none")
  k2= subset(e2, e2$Var1 == "none")
  
  if (nrow(k1) == 0 | nrow(k2) == 0){
    Alleles = unique(c(b$allele3, b$allele4))
    b$New_allele3 <- Alleles[1]
    b$New_allele4 <- Alleles[2]
    DataStore <- rbind(DataStore,b)
    next
  }

  if (k1$Freq > k2$Freq) {
    e = e2
    e= subset(e, e$Var1 != "none")
    if (nrow(e) < 2){
      b$New_allele3 <- b$allele3
      b$New_allele4 <- b$allele4
      DataStore <- rbind(DataStore,b)
      next
    }
    e=e[order(e$Freq),]
    Allele1 = as.character(e$Var1[[2]])
    Allele2 = as.character(e$Var1[[1]])
    b$New_allele3 <- (ifelse(b$allele3 != "none", Allele1, "none"))
    b$New_allele4 <- (ifelse(b$allele4 != "none", Allele2, "none"))
    
  } else {
    e = e1
    e= subset(e, e$Var1 != "none")
    if (nrow(e) < 2){
      b$New_allele3 <- b$allele3
      b$New_allele4 <- b$allele4
      DataStore <- rbind(DataStore,b)
      next
    }
    e=e[order(e$Freq),]
    Allele1 = as.character(e$Var1[[2]])
    Allele2 = as.character(e$Var1[[1]])
    b$New_allele3 <- (ifelse(b$allele3 != "none", Allele1, "none"))
    b$New_allele4 <- (ifelse(b$allele4 != "none", Allele2, "none"))
  }
  
  DataStore <- rbind(DataStore,b)
}

write.csv(DataStore,"AllFilesCombinedCorrect.csv",row.names = F,quote = FALSE)





