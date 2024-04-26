# This bash code was used to trim and crop FASTQ files as 
# well as to remove adapters using Trimmomatic software. 
#=======================================================================#
PATH_to_INPUTS=path
FileListOfSamples=file.txt
#=======================================================================#
k=0
while read line
do 
k=$((k+1))
echo "=================================================================="
echo "******************************************************************"
echo Action: Trimming 
echo The total number of samples = $(wc -l < "${FileListOfSamples}") 
echo The run is for sample: $k 
echo Sample ID: $line
echo "******************************************************************"
echo "=================================================================="
#=======================================================================#
java -jar trimmomatic-0.39.jar PE -threads 50 \
${line}_1.fastq.gz \
${line}_2.fastq.gz \
${line}_1.trim.fastq.gz \
${line}.trimUnp.fastq.gz \
${line}_2.trim.fastq.gz 
#=======================================================================#
done < ${FileListOfSamples}