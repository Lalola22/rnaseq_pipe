# Analysis of RNA-seq data

A basic pipeline for analysis of RNA-seq data. This provides a framework for reference guided quantification and differential expression testing including both quality control steps on read files and trimming of adapters / low quality reads.

 This pipeline has only been tested on RNA-seq data generated on the Illumina HiSeq platform.

**Note:** All programs used need to be added to `$PATH`

## Usage

`./wrapper.sh "/path/to/data" "/path/to/results/dir" "number of cores"`

* The first argument points to a directory containing condition-named directories with the paired end sequencing reads
    * e.g. `data/` may contain `{foo/,bar/}` and then `data/foo/` has `foo_1_R1.fastq.gz  2foo_1_R2.fastq.gz  foo_2_R1.fastq.gz  foo_2_R2.fastq.gz  foo_3_R1.fastq.gz  foo_3_R2.fastq.gz`
* The second arguement is the location you want results to be saved to. If the directory does not exist it will be created
* The third arguement is the number of threads to allocate to this pipeline. The more the better: vroom vroom

An example call that also uses piping to `tee` to save a log of commands run is:

`./wrapper.sh "/home/slee/data/bcl6_raw" "/home/slee/outputs/bcl6_paper" "18" | tee 2017-06-10-log.txt`

## Required meta-data

### `extra_paths.txt`

This file should contain the following paths with trailing `/`

* path to the trimmomatic directory
* path to the reference transcriptome
    *  A reference transcriptome can be downloaded from the UCSC genome browser or from ensembl. The pipeline has been tested using ensembl GRCh38 release 88/89

### `TruSeq3-PE.fa`

This fasta files contains the adapter sequences for Trimmomatic to remove. By default contains the Illumina adapters.

### `condition_table.txt`

This is a csv table that lists each of the conditions and how many replicates there are. Used for creating the DESeq2 sample table

### `de_tests_table.txt`

This is a csv table that lists each of the log2 fold change tests to be run by DESeq2. Each line of this file is a test that will be conducted and saved to its own folder.


## Outputs

* FastQC Reports for before and after read trimming
* Trimmed `.fastq.gz` files for each pair of raw-reads.
* Kallisto quantifications for each sample
* Salmon quantifications for each sample
* DESeq2 outputs, a number of plots are generated as well as the results file in `.csv` for further downstream visualisations


## Dependencies

While this pipeline does not install required software there are two file `checks/check-packages.R` and `checks/check-programs.sh` that can be run as a quick test that all programs are installed. If either of these fail look at the error message to work out what packages or programs are stil required.
