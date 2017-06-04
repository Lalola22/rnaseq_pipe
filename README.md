# Analysis of RNA-seq data

All code used to analyse RNA-seq data for $PAPERNAME are contained within this repository.

Running the script `wrapper.sh` on raw-read `.fastq.gz` data in an environment with all dependencies installed will be sufficient to recreate analysis.

All programs used need to be added to the PATH

### Usage

`./wrapper.sh "/path/to/data" "/path/to/results/dir" "number of cores"`

The call to generate the data for the paper was

`./wrapper.sh "/home/slee/data/test_data" "/home/slee/outputs/bcl6_paper_test" "18" | tee log.txt`

### extra_paths.txt

This file should contain the following paths with trailing `/`

* path to trimmomatic directory
* path to the reference transcriptome

### sleuth_samples.txt

This is a TSV file that contains the different comparisons that should be made within sleuth for both of the Salmon and kallisto data. Column one is the treatment and column two the control.
