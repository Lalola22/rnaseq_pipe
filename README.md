# Analysis of RNA-seq data

Written by Samuel Lee: samleenz@me.com
Forked version by Lalola22 for McConnell group

A basic pipeline for analysis of RNA-seq data. This provides a framework for reference guided quantification and differential expression testing including both quality control steps on read files and trimming of adapters / low quality reads.

* `FastQC` for fastq file quality checking
* `trimmomatic` for fastq read trimming
* `Kallisto` for pseudo-alignment and quantifications
* `tximport` to prep quantifications for DE testing
* `DESeq2` for testing of differential expression

* investigating incorportaton of additional qc to control for sample bias

 This pipeline has only been tested on RNA-seq data generated on the Illumina HiSeq platform.

While every effort has been made to ensure this works as intended there are no guarantees of accuracy / vailidity of results generated etc.

**Notes:**

* All programs used need to be added to `$PATH` (the dependencies check will pick up if they're not)
* At this stage assumes your RNA-seq reads are Human... will update with an option to set the species for the transcript -> gene mapping soon. The transcript level data should be fine at this point if it manages to run through despite the species issue.

## Usage

`./run_bioinformatics.sh "/path/to/data" "/path/to/results/dir" "number of cores"`

* The first argument points to a directory containing condition-named directories with the paired end sequencing reads
    * e.g. `data/` may contain `{foo/,bar/}` and then `data/foo/` has `foo_1_R1.fastq.gz  foo_1_R2.fastq.gz  foo_2_R1.fastq.gz  foo_2_R2.fastq.gz  foo_3_R1.fastq.gz  foo_3_R2.fastq.gz`
* The second arguement is the location you want results to be saved to. If the directory does not exist it will be created
* The third arguement is the number of threads to allocate to this pipeline. The more the better: vroom vroom

One approach to keeping the scripts used to create results associated with them will be to create a fresh clone of this repository prior to running the pipeline and saving the results within it e.g.

```bash
# clone the git repo and give a descriptive name
git clone https://github.com/samleenz/rnaseq_pipe.git
mv rnaseq_pipe YYYY-MM-DD-rnaseq_pipe-experiment-details
cd YYYY-MM-DD-rnaseq-pipe-experiment-details

# STOP: edit the aux files to match your experimental design

# check all dependencies

Rscript checks/check-packages.R

./checks/check-programs.sh

# fix issues and re-run above to resolve dependencies

datDir="/home/slee/data/bcl6_raw"
outDir="outputs"
cores="18"

./run_bioinformatics.sh $datDir $outDir $cores | tee log.txt
```


## Required meta-data

All meta-data files have example content within.

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
    * Forward / reverse paired and forward / reverse unpaired
* kallisto quantifications with bootstrapping (n = 100) for each sample
* DESeq2 outputs at gene and transcript level
    * a number of plots are generated as well as the results file in `.csv` format for further downstream visualisations


## Dependencies

While this pipeline does not install required software there are two file `checks/check-packages.R` and `checks/check-programs.sh` that can be run as a quick test that all programs are installed. If either of these fail look at the error message to work out what packages or programs are stil required and refer to google or contact me. Currently trimmomatic install status is not tested by `check-programs.sh`
