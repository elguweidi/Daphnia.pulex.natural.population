# This bash code was used to extract the exons from the GTF file
gunzip Daphnia_pulex.V1.0.52.gtf.gz
awk '$3 == "exon"' Daphnia_pulex.V1.0.52.gtf > temp1
awk '{print $3, $4 ,$5 ,$10 ,$12 ,$14}' temp1 > temp2
sed -i 's/";//g' temp2
sed -i 's/"//g' temp2

# Add a header to the temp2:
sed -e -i '1i\Exon Start End Gene_id Transcript_id Exon_number' temp2

# Now, provide a list of single copy genes, which is saved in file "listOfsingleCopyConservedGenes.txt",
# then sort it and join it with file temp2, as follows: 
dos2unix listOfsingleCopyConservedGenes.txt
LANG=en_EN sort -k 1,1 listOfsingleCopyConservedGenes.txt > listOfsingleCopyConservedGenesSorted.txt
LANG=en_EN sort -k 4,4 temp2 > temp3
LANG=en_EN join -1 1 -2 4 listOfsingleCopyConservedGenesSorted.txt temp3 > temp4
