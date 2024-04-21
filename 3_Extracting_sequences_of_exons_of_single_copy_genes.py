# This Python code was used to extract the sequences of 
# single copy exons of interest withdesired length (300 to 1500pb)
# from the mart_export.txt.gz file, which contains sequences of all exons.
#===========================================================================#
ListofLookUpExons = []
with open("listOfSCgenes.txt","r") as inFile:
    for exonNmae in inFile:
        ListofLookUpExons.append(exonNmae.strip())
#===========================================================================#
#************************************#
startLength = 300
endLength = 600
#************************************#
scgFile = open("scgFile.fa", 'w')
StoreSeq = ""
exonIDs = []
k=0
with open("mart_export.txt","r") as fi:
    for line in fi:
        if not line.startswith(">"):
            StoreSeq = StoreSeq + line.strip()
        if line.startswith(">"):
            k=k+1
            exonIDs.append(line.strip())
            a = StoreSeq
            StoreSeq = ""
            if ((len(a)>= startLength) &  (len(a) <= endLength) & (exonIDs[k-2].strip() in ListofLookUpExons)):
                print(exonIDs[k-2].strip(), file = scgFile)
                print(a, file = scgFile)
scgFile.close()
#===========================================================================#
