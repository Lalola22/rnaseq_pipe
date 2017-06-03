# /bin/bash
set -e
set -u
set -o pipefail

# usage: scriptname.sh "/path/to/fastq/dirs" "/path/to/output/dir"
# run from top of repository
# $cores needs to be adjusted for individual machines

rawDir="$1"
outDir="2"
cores="16"

# test fastq files are readable

if [ ! -d $rawDir ]
then
    echo "${rawDir} is not a directory"
    exit 1
fi


# environment prep
if [ ! -d $outDir ]
then
    echo mkdir $outDir
fi

fastqArray=($(find $rawDIR -type f -name "*.fastq.gz"))

# Pre-trim QC

fastqc -o "${outDir}/fastqc/raw" -t $cores ${fastqArray[@]]}


# Adaptor trimming



# Trimmed read QC


# Pseudo
