# /bin/bash
set -e
set -u
set -o pipefail

# usage: scriptname.sh "/path/to/fastq/dirs" "/path/to/output/dir" "trimmomaticPath" "n_cores"
# run from top of repository
# $cores needs to be adjusted for individual machines

# Parameters

rawDir="$1"
outDir="$2"
trimmomaticPath="$3"
cores="$4"

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

python3 batch_trim.py "${outDir}/" "$rawDIR" "$trimmomaticPath"
### trimmomatic .fa files needs to be added to aux_files


# Trimmed read QC

fastqArrayTrim=($(find "${outDir}/batch_trim/processed" -type f -name "*P.fastq.gz"))

fastqc -o "${outDir}/fastqc/trimmed" -t $cores ${fastqArrayTrim[@]]}

# kallisto

kallisto_quant.py "${outDir}/" "${outDir}/batch_trim/processed/"

# Salmon



# Salmon prep with Wasabi


# Sleuth analysis


# Create heatmaps
