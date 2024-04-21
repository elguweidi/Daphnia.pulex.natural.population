# This bash code was used to obtain copy number based on single copy genes.
#========================================================================================#
# Path of input folder
PATH_to_INPUTS=path
# Path of OUTPUT folder
PATH_to_OUT=path
# Path to file with list of samples
FileListOfSamples=file.txt
#========================================================================================#
if [ -d "Summaries_CNV" ]; then rm -Rf Summaries_CNV; fi 
mkdir Summaries_CNV
k=0
while read line
do 
k=$((k+1))
echo "=================================================================================="
echo "**********************************************************************************"
echo Action: Copy number variation
echo The total number of files in list is: $(wc -l < "${FileListOfSamples}")
echo The current run is for sample number: $k
echo The sample ID is: $line 
echo "**********************************************************************************"
echo "=================================================================================="

#=====================================================================#
# Copy number variation:
if [ -d "${line}" ]; then rm -Rf ${line}; fi 
mkdir ${line}
#cp ${line}.coverage.txt ${line}
cd ${line}
bedtools genomecov -d -ibam ${PATH_to_INPUTS}${line}.dupRem.bam > ${line}.coverage.txt
#========================================================================================#
awk -F '\t' '{print>$1}' *.coverage.txt
mv -f *.coverage.txt ${PATH_to_OUT}Summaries_CNV
for f in *; do sh stats2.sh $f >> summary_${line}.txt; done
cat summary_${line}.txt >> ${PATH_to_OUT}Summaries_CNV/All_Summary.txt
#========================================================================================#
cp summary_${line}.txt ${PATH_to_OUT}Summaries_CNV
echo Done the sample: $line 
echo "=================================================================================="
#========================================================================================#
cd ${PATH_to_OUT}
done < ${FileListOfSamples}



