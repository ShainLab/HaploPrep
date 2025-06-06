# HaploPrep

HaploPrep_v1.2.2

## What does this app do?
1. Initial filtering : Processes FreeBayes variant calls to identify heterozygous SNPs. Filters for SNPs with a population allele frequency ≥ 1% in the 1000 Genomes database and read allele frequency between 30-70% (heterozygous) based on FreeBayes estimates.  
2. Depth analysis : Runs samtools mpileup on filtered SNPs to check read depth at each variant position and count reference and alternate allele reads.  
3. Refinement : Further filters SNPs based on actual read data; filters for variants with read allele frequency between 40-60% and read depth ≥ 10. These are the "Informative SNPs".  
4. Phasing : Uses GATK ReadBackedPhasing tool to phase these "Informative SNPs", determining which variants are on the same haplotype based on read information. Extracts chromosome position, reference/alternate alleles, and phased allele assignments indicating which allele is on each haplotype.  

- - -
Note: Other variant callers aside from FreeBayes could be used, as long as they produce a vcf of candidate SNPs with similar column headers and structure. We provide sample input files to mimic for the user's convenience; please refer to test-dataset.
- - -

## What does this app output?
•	InformativeSNPs.txt : txt file of the FreeBayes data after initial filtering (output from step 1)  
•	MpileupOutput.txt : txt file containing results from samtools mpileup analysis (output from step 2)  
•	InformativeSNPs.vcf : vcf file created after the mpileup depth analysis and refinement containing only SNPs that passed additional filtering (output from step 3)  
•	MpileupInformativeSNPInput.txt : txt file containing the filtered SNPs used for phasing (output from step 3)  
•	phased_SNPs.vcf : vcf file output from GATK ReadBackedPhasing tool (output from step 4)  
•	InputVariantsforphasing.txt : txt file containing haplotype-resolved SNPs distinguishing maternal and paternal allelic combinations (output from step 4)

## What data is required for this app to run?
- DNAbam - A BAM file containing the DNA sequencing reads
- FBvcf - A VCF file from FreeBayes containing the initial variant calls
- FBtxt - A text file containing annovar annotation of the FreeBayes vcf

## What are typical use cases for this app?
Please refer to the following publication and preprints for use cases:
1. Tang J, Fewing, Chang D, Zeng H, Liu S, Jorapur A, Belote RL, McNeal AS, Tan T, Yeh I, Arron ST, Judson-Torres RL, Bastian BC, Shain AH. The genomic landscapes of individual melanocytes from human skin. Nature. October 2020. https://doi.org/10.1038/s41586-020-2785-8.
2. Tandukar B\*, Deivendran D\*, Chen L, Cruz-Pacheco N, Sharma H, Xu A, Bandari AK, Chen DB, George C, Marty A, Cho RJ, Cheng J, Saylor D, Gerami P, Arron ST, Bastian BC, Shain AH. Genetic evolution of keratinocytes to cutaneous squamous cell carcinoma. bioRxiv. July 2024. https://doi.org/10.1101/2024.07.23.604673.
3. Tandukar B, Deivendran D, Chen L, Bahrani N, Weier B, Sharma H, Cruz-Pacheco N, Hu M, Marks K, Zitnay RG, Bandari AK, Nekoonam R, Yeh I, Judson-Torres R, Shain AH. Somatic mutations distinguish melanocyte subpopulations in human skin. bioRxiv. February 2025. https://doi.org/10.1101/2025.02.07.637114.

_*jointly led the project._
