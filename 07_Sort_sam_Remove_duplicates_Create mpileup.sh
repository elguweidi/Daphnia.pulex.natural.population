# This bash code was used to sort the same file, creates a bam file 
# and then drop duplicate reads, using [picard.jar MarkDuplicates] 
# and then creates mpileup files
#=================================================================#
# Path of input folder
PATH_to_INPUTS=path
# Reference
TheReference=Reference.fa
# Path to file with list of sample names
FileListOfSamples=file.txt
#=================================================================#
k=0
while read line
do 
k=$((k+1))
echo "============================================================"
echo "************************************************************"
echo Actions: 
echo 1. Sort sam files and generate bam files  
echo 2. Remove duplicates 
echo 3. Create mpileup file and variant calling 
echo The total number of files in list is: $(wc -l < "${FileListOfSamples}")
echo The current run is for sample number: $k
echo The sample ID is: $line 
echo "************************************************************"
echo "============================================================"
# Sorting
#=================================================================#
# Sort the sam and get bam
samtools view -bS -q 20 \
${PATH_to_INPUTS}${line}.sam | samtools sort > \
${PATH_to_INPUTS}${line}.bam
#=================================================================#
# Remove duplicates:
#=================================================================#
java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
      REMOVE_DUPLICATES=true \
      I=${PATH_to_INPUTS}${line}.bam \
      O=${PATH_to_INPUTS}${line}.dupRem.bam \
      M=${PATH_to_INPUTS}${line}.dup.txt
#=================================================================#
# Create mpileup file and variant calling (Variant identification):
bcftools mpileup -a FORMAT/AD,FORMAT/DP \
-d 10000 \
-f $TheReference ${PATH_to_INPUTS}${line}.dupRem.bam \
 | bcftools \
call -A -m -mv -Ov -o ${line}.vcf
#=================================================================#
# Get needed results
bcftools query -f '%CHROM %POS %REF %ALT %DP [ %AD]\n' ${line}.vcf > Results_${line}.vcf
#=================================================================#
# Table the results
sed  -i '1i CHROM POS REF ALT DP AD' Results_${line}.vcf
#=================================================================#
done < ${FileListOfSamples}


