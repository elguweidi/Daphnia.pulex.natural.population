SCRIPTS and CODES used for the analysis of rDNA variation in Daphnia pulex

1. Obtain single copy genes from EnsemblMetazoa (https://metazoa.ensembl.org/index.html):
We obtained a list of 12,295 Daphnia pulex single copy genes, which was saved in file [SingleCopyGenesFromEnsemble.txt]. These genes were then matched with a list of conserved genes (ConservedGenes.txt) obtained from Colbourne et al., (2011), which are genes conserved as single-copy orthologs across eukaryotic genomes. The R code [1_match_SSCgenes_ConservedGenes.R] was used for this step. This resulted in 256 genes, which we considered to be the single copy genes.

2. Extract exons of single copy genes from EnsemblMetazoa (https://metazoa.ensembl.org/index.html):
From EnsemblMetazoa, we downloaded a GTF file [Daphnia_pulex.V1.0.52.gtf] containing all information needed to extract the exons. The Linux code [2_extraction_exons.sh] was used to extract the exons from the GTF file. The final output was saved in an Excel file for further analysis.

3. Extract the sequences of exons of single copy genes:
The file named [mart_export.txt.gz] was downloaded from EnsemblMetazoa and contains sequences of all exons. First, a file [listOfSCgenes.txt], which includes a list of the name of the genes and the exon IDs separated by '|', was generated. Next, a Python code [3_Extracting_sequences_of_exons_of_single_copy_genes.py] was generated to extract the sequences of single copy exons of interest with desired length (e.g., length from 300 to 1500 bp). The outputs were saved in a fasta file [scgFile.fa]. The list in this file contained 228 exons that had length between lengths specified in the code (i.e., 300 to 1500). However, some of them are from the same genes. We only kept the longest exon from each gene.  

4. Remove additional exons: 
Additional exons were removed from the file [scgFile.fa] using the Python code [4_Drop_exons_from_seme_gene.py]. The output file [scgenes.fa] contains the final unique list, which was 167 exons.

5. Trim and crop FASTQ files:
FASTQ files (sample.fastq.gz) containing the sequences of Daphnia samples were trimmed and adapters were removed using Trimmomatic software, which was run using the bash code [Trimming.sh]. The output files are the trimmed files (samples.trim.fastq.gz).  

6. Map genome sequences against a reference genome:
The [Alignment.sh] bash code was used for mapping trimmed samples sequences (samples.trim.fastq.gz) against a reference genome. The mapping/alignment was done by the BWA mem software package and the outputs were SAM files (sample.sam)

7. Sorting SAM files, remove duplicate reads and create mpileup files:
The code [7_Sort_sam_Remove_duplicates_Create mpileup.sh] takes the SAM files from the alignment step and performs sorting of the SAM files, creates BAM files, drops duplicated reads (using [picard.jar MarkDuplicates]) and then generates mpileup files (using [bcftools]) for each sample in the list of samples [file.txt]. The first output files are [sample.dupRem.bam] and the second output files are the [Results_sample.vcf] files.

8. Obtain results needed to calculate rDNA regions copy number: 

The Bash code [8_CopyNumberVariation.sh] takes [sample.dupRem.bam] files as inputs and using bedtools genomecov generates the genome coverage for each sample and saved in files [sample.coverage.txt]. Then this coverage files are used as inputs to obtain outputs needed for rDNA copy number following the method proposed by Colbourne at al., (2011). The final outputs are saved in files [summary_sample.txt] for each sample and in [All_Summary.txt] for all samples together in a single file.
 
9. Calculation of rDNA copy number:
The R code [9_CopyNumberCalculation.R] takes the [All_Summary.txt] files and calculates rDNA haploid and diploid copy number. The output file is [cn.csv].

10. Tabulate variant calling results: 
The R code [10_Tibble_VariantCallingResults.R] was used to tabulate the variant calling results based on the vcf files obtained from step # 7 above. The outputs are saved in a single file named [Combined_VC.csv].

11. Correction alleles 3 and 4:
The R code [11_CorrectAlleles_3_and_4.R] was used to correct third and fourth alleles at SNP loci if they differ between individuals. The input file is [Combined_VC.csv], and the output file is [AllFilesCombinedCorrect.csv].



12. Variant calling data pruning:
The R code [12_dataPruning.R] was used to remove SNPs whose mean frequency of allele 1 within a population was < 0.010 or > 0.990 from the variance calling data [AllFilesCombinedCorrect.csv]. The output file is [PrunedData.csv].

13. Haplotype analysis: 
The R code [13_haploSep.R] was used to run the haplotype analysis using the R function haplosep. The input file is [PrunedData.csv] and the output files are: 
I. haploFrq.a1.csv
II. haploStr.a1.csv
III. HaplotypeGenotype.csv
IV. HaplotypeGenotypeT.csv
V. hapFile.txt
VI. df.genpop2HapTab.txt

14. Plot a neighbor-joining tree: 
The R code [14_Neighbor-Joining.R] was used to plot the neighbor-joining tree. The input file is [square-distance-matrix.txt] and the output file is [phylogram.png].

References
Colbourne JK, et al.2011. The ecoresponsive genome of Daphnia pulex. Science 331:555-561.
