# Analysis of RNA-seq data

All code used to analyse RNA-seq data for $PAPERNAME are contained within this repository.

Running the script `wrapper.sh` on raw-read `.fastq.gz` data in an environment with all dependencies installed will be sufficient to recreate analysis.

All programs used need to be added to the PATH

### Usage

The call to generate the data for the paper was

`wrapper.sh "/home/slee/raw_zipped" "/home/slee/outputs/bcl6_paper" "18"`

### extra_paths.txt

This file should contain the following paths with trailing `/`

* path to trimmomatic directory
* path to the reference transcriptome
