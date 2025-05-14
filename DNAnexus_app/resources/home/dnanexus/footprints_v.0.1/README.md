# Footprints_v1

### Footprints runs a c based tool, foots, that uses samtools to count the number of bases in your bam file with a minimum coverage of your choice.

##            Welcome to footprints

Thank you for downloading footprints v.0.1.

Footprints counts the number of bases in your bam file with a minimum coverage of your choice.

The sequencing footprint of your bam can then be used to accurately calculate mutation burden.

We recommend you select a minimum coverage value `-n` that reflects the minimum number of reads your 
pipeline requires to call variants.


## Installation

To install please make sure you have the following requirements


### Requirements

1. samtools-0.1.19
2. bedtools

If you don't have this version of samtools installed please follow steps below to install and add to path. 
Otherwise skip to the next step.

### Install Samtools v.0.1.19

You make require --no-check-certificates here. 
If you are not comfortable with using this link please find an alternative download online.

```
$ wget https://sourceforge.net/projects/samtools/files/samtools/0.1.19/samtools-0.1.19.tar.bz2

$ tar -xf samtools-0.1.19.tar.bz2

$ cd samtools-0.1.19

$ make

$ export PATH=/path/to/samtools-0.1.19:$PATH
```

### Install Footprints v.0.1

```
$ ./configure
$ make
```

## Run Footprints v.0.1

```
$ ./footprints -h

Program: Footprint calculation (for accurately calculating sequencing footprint)

Version: 0.1

Usage: ./footprints -b <bam> -i <intervals bed file> -n <minimum coverage(default=5)> -o <output file>

Options:
        -b      Bam: Bam file or directory containing bam files on which footprints are to be calculated [required]
        -i      Intervals bed file: Genomic intervals on which to operate [required]
        -n      Minimum coverage: Minimum number of reads a base needs to be covered by to count as part of callable footprint [default=5]
        -o      Output file: Text output for results [default=STDOUT]
        -h      Help: Does what it says on the tin
```

## Input

Any input bams must be indexed.

You can input one bam file: `-b /path/to/sample1.bam`

Or a directory containing multiple bams: `-b /path/to/dir`

## Output

By default, each bamfiles footprint is written to standard out:
    Calculating footprint on file: sample1.bam
    Minimum coverage selected: 5
    Number of base pairs in file with greater than 5 coverage: 4910521
Or you can output the results into a tab seperated file of your choice.
You can select this option by adding `-o output.txt`

BAM | BP
--- | --
sample1.bam | 4910521
sample2.bam | 4906269
sample3.bam | 4930272

