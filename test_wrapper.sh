# /bin/bash
set -e
set -u
set -o pipefail

# usage: scriptname.sh "/path/to/fastq/dirs" "/path/to/output/dir" "n_cores"
# run from top of repository
# $cores needs to be adjusted for individual machines

# Parameters

rawDir="$1"
outDir="$2"
cores="$3"

# test fastq files are readable

if [ ! -d $rawDir ]
then
    echo "${rawDir} is not a directory"
    exit 1
fi


# -- environment prep

# make results directory if required

if [ ! -d $outDir ]
then
    echo mkdir $outDir
fi

# Initialise bash array with raw-read paths

fastqArray=($(find $rawDir -type f -name "*.fastq.gz"))

# Read the two auxillary file into bash arrays

oldIFS="$IFS"
IFS=$'\n' extraPaths=($(<"extra_paths.txt"))
IFS="$oldIFS"

oldIFS="$IFS"
IFS=$'\n' samples=($(<"sleuth_samples.txt"))
IFS="$oldIFS"

# Set paths

trimmomaticPath=${extraPaths[0]}
GRCh38trans=${extraPaths[1]}

export PATH=${PWD}:$PATH

# -- Adaptor trimming

echo "Trimming adaptors..."

python3 batch_trim.py "${outDir}/" "${rawDir}/" "${trimmomaticPath}/" $cores


# -- Trimmed read QC

echo "Post-trimming read QC"

# Initialise bash array for trimmed + paired reads

fastqArrayTrim=($(find "${outDir}/batch_trim" -type f -name "*P.fastq.gz"))

mkdir -p "${outDir}/fastqc/trimmed"

fastqc -o "${outDir}/fastqc/trimmed" -t $cores ${fastqArrayTrim[@]]}
