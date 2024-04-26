# This bash code is for mapping samples sequences against a reference genome.
# The mapping/alignment was done by BWA mem software package
#=====================================================================#
# Path to input files
PATH_to_INPUTS=path
# Reference
TheReference=Reference.fa
# Path to file with list of sample names
FileListOfSamples=${PATH_to_INPUTS}file.txt
#=====================================================================#
k=0
while read line
do 
k=$((k+1))
echo "================================================================"
echo "****************************************************************"
echo Action: Alignment 
echo The total number of files in list is: $(wc -l < "${FileListOfSamples}")
echo The current run is for sample number: $k
echo The sample ID is: $line 
echo "****************************************************************"
echo "================================================================"
#=====================================================================#
bwa mem -v 1 -t 12 \
$TheReference \
${PATH_to_INPUTS}${line}_1.trim.fastq.gz \
${PATH_to_INPUTS}${line}_2.trim.fastq.gz > \
${PATH_to_INPUTS}${line}.sam
#=====================================================================#
done < ${FileListOfSamples}