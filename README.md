# Analysis of RNA-seq data

All code used to analyse RNA-seq data for $PAPERNAME are contained within this repository.

Running the script `wrapper.sh` on raw-read `.fastq.gz` data in an environment with all dependencies installed will be sufficient to recreate analysis.

All programs used need to be added to the PATH

## Usage

`./wrapper.sh "/path/to/data" "/path/to/results/dir" "number of cores"`

* The first arg so point to a directory containing condition-named directories with the paired end sequencing reads
    * e.g. `data/` may contain `25uM/ DMSO/ Untreated/` and then `data/25uM/` has `25uM_1_R1.fastq.gz  25uM_1_R2.fastq.gz  25uM_2_R1.fastq.gz  25uM_2_R2.fastq.gz  25uM_3_R1.fastq.gz  25uM_3_R2.fastq.gz`
* The second arg is the location you want results to be saved to
* The third arg is the number of threads to allocate to this pipeline. The more the better: vroom vroom

The call to generate the data for the paper was

`./wrapper.sh "/home/slee/data/bcl6_raw" "/home/slee/outputs/bcl6_paper" "18" | tee log.txt`

## Required meta-data

### extra_paths.txt

This file should contain the following paths with trailing `/`

* path to the trimmomatic directory
* path to the reference transcriptome
    *  A reference transcriptome can be downloaded from the UCSC genome browser

### sample_table.txt

This is a TSV file that contains the different comparisons that should be made within DESeq2 for the Salmon data. Column one is the treatment and column two the control. It should look something like this. Each row of the file will result in a new set of calls to DESeq2 for the Salmon results to be made.

```
Treatment   Control Replicates
25uM    DMSO    3
DMSO    Untreated   3
Dox Untreated   3

```

## Outputs

* FastQC Reports for before and after read trimming
* Trimmed `.fastq.gz` files for each pair of raw-reads.
* Kallisto quantifications for each sample
* Salmon quantifications for each sample
* DESeq2 outputs, a number of plots are generated as well as the results file in `.csv` for further downstream visualisations


## Dependencies

If you can run all these programs you should be good to go..

* FastQC
* Trimmomatic
* salmon

The following R packages will be needed

* biomaRt
* tidyverse
* tximport
* DESeq2
* stringR

Python 3.XX is required with the following non-standard packages

* joblib
