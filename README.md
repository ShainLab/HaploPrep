# HaploPrep (DNAnexus Platform App)
This app identifies high-quality heterozygous SNPs (Single Nucleotide Polymorphisms) and prepares them for haplotype phasing by filtering variants based on population frequency, allelic balance, and sequencing depth from DNA-seq data.

_Note: Two outputs generated by this app (MpileupInformativeSNPInput.txt and InputVariantsforphasing.txt) are required inputs to run our [Single_Cell_Somatic_Mutation_Caller](https://github.com/ShainLab/Single_Cell_Somatic_Mutation_Caller) app._

- - -
_**This repository is for an app that runs on the DNAnexus Platform.**  
This app has been made publicly available on DNAnexus. To run, navigate to the Tools Library and select the app (HaploPrep). For more information about how to run or modify it, see https://documentation.dnanexus.com/.  
If you would like you clone the app and build it on your own DNAnexus account, please refer to [this section](https://github.com/ShainLab/HaploPrep/blob/main/README.md#how-to-clone-and-build-this-app-to-run-on-dnanexus)._
- - -

## What does this app do?
1.	<ins>Initial filtering</ins> : Processes FreeBayes variant calls to identify heterozygous SNPs. Filters for SNPs with a population allele frequency ≥ 1% in the 1000 Genomes database and read allele frequency between 30-70% (heterozygous) based on FreeBayes estimates.  
2.	<ins>Depth analysis</ins> : Runs samtools mpileup on filtered SNPs to check read depth at each variant position and count reference and alternate allele reads.  
3.	<ins>Refinement</ins> : Further filters SNPs based on actual read data; filters for variants with read allele frequency between 40-60% and read depth ≥ 10. These are the "Informative SNPs".  
4.	<ins>Phasing</ins> : Uses GATK ReadBackedPhasing tool to phase these "Informative SNPs", determining which variants are on the same haplotype based on read information. Extracts chromosome position, reference/alternate alleles, and phased allele assignments indicating which allele is on each haplotype.  

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

- - -
Note: Other variant callers aside from FreeBayes could be used, as long as they produce a vcf of candidate SNPs with similar column headers and structure. We have tested our pipeline with HaplotypeCaller; please refer to [benchmarking_data](https://github.com/ShainLab/HaploPrep/tree/main/benchmarking_data) for the comparison metrics.  
We provide sample input files to mimic for the user's convenience; please refer to [test-dataset](https://github.com/ShainLab/HaploPrep/tree/main/test-dataset).
- - -

## What are typical use cases for this app?
Please refer to the following publication and preprints for use cases:
1. Tang J, Fewing, Chang D, Zeng H, Liu S, Jorapur A, Belote RL, McNeal AS, Tan T, Yeh I, Arron ST, Judson-Torres RL, Bastian BC, Shain AH. The genomic landscapes of individual melanocytes from human skin. Nature. October 2020. https://doi.org/10.1038/s41586-020-2785-8.
2. Tandukar B\*, Deivendran D\*, Chen L, Cruz-Pacheco N, Sharma H, Xu A, Bandari AK, Chen DB, George C, Marty A, Cho RJ, Cheng J, Saylor D, Gerami P, Arron ST, Bastian BC, Shain AH. Genetic evolution of keratinocytes to cutaneous squamous cell carcinoma. bioRxiv. July 2024. https://doi.org/10.1101/2024.07.23.604673.
3. Tandukar B, Deivendran D, Chen L, Bahrani N, Weier B, Sharma H, Cruz-Pacheco N, Hu M, Marks K, Zitnay RG, Bandari AK, Nekoonam R, Yeh I, Judson-Torres R, Shain AH. Somatic mutations distinguish melanocyte subpopulations in human skin. bioRxiv. February 2025. https://doi.org/10.1101/2025.02.07.637114.

_*jointly led the project._

## How to clone and build this app to run on DNAnexus?

1. Clone this repository:
  ```
  git clone https://github.com/ShainLab/HaploPrep.git
  ```
2. Move and rename "DNAnexus_app" subdirectory to "HaploPrep" (test-dataset should not be included). _Note: This is to ensure the proper sturucture for DNAnexus app build_:
  ```
  cp -r HaploPrep/DNAnexus_app ./HaploPrep
  ```
3. **[IMPORTANT!]** Due to GitHub file size limits, please do the following before building the app:

&nbsp;&nbsp;&nbsp;&nbsp; a. Move into this subdirectory:
 ```
  cd HaploPrep/resources/home/dnanexus/hg19
  ```
&nbsp;&nbsp;&nbsp;&nbsp; b. Uncompress ucsc.hg19.fasta and remove the tarball:
  ```
  tar -xzf ucsc.hg19.fasta.tar.gz &&  rm ucsc.hg19.fasta.tar.gz
  ```
&nbsp;&nbsp;&nbsp;&nbsp; c. Uncompress hg19.fa and remove the tarball:
  ```
  tar -xzf hg19.fa.tar.gz && rm hg19.fa.tar.gz
  ```
&nbsp;&nbsp;&nbsp;&nbsp; d. Move out of the current directory, create a tarball of the hg19 directory, then remove the directory:
  ```
  cd ../ ; tar -czf hg19.tar.gz hg19 && rm -rf hg19
  ```
4. You can now begin building this app on DNAnexus using the following commands:

&nbsp;&nbsp;&nbsp;&nbsp; a. Install dx-toolkit (if not already installed): [dx-toolkit installation](https://documentation.dnanexus.com/downloads)  
&nbsp;&nbsp;&nbsp;&nbsp; b. Move into the root directory:
```
cd HaploPrep
```  
&nbsp;&nbsp;&nbsp;&nbsp; c. Run the build command. _Replace ProjectName:path/to/your/directory/ with your DNAnexus project name and path for app to build._
```
dx build -d ProjectName:path/to/your/directory/
``` 

You should see the app in this directory when done. Select it, add the input files to their respective input boxes, under Workflow Actions set your output folder, then click start as analysis.


