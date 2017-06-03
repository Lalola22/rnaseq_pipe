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


# environment prep
if [ ! -d $outDir ]
then
    echo mkdir $outDir
fi

fastqArray=($(find $rawDIR -type f -name "*.fastq.gz"))

oldIFS="$IFS"
IFS=$'\n' extraPaths=($(<"extra_paths.txt"))
IFS="$oldIFS"

export PATH=${PWD}:$PATH

# Pre-trim QC

echo "Inital read QC..."

fastqc -o "${outDir}/fastqc/raw" -t $cores ${fastqArray[@]]}


# Adaptor trimming

echo "Trimming adaptors..."
echo "Outputs will be sent to log"

echo python3 batch_trim.py "${outDir}/" "$rawDIR" "$trimmomaticPath"
### trimmomatic .fa files needs to be added to aux_files


# Trimmed read QC

echo "Post trim read QC"

echo fastqArrayTrim=($(find "${outDir}/batch_trim/processed" -type f -name "*P.fastq.gz"))

echo fastqc -o "${outDir}/fastqc/trimmed" -t $cores ${fastqArrayTrim[@]]}

# kallisto index

echo "Creating the kallisto index..."

echo kallisto index -i "aux_files/GRCh38transcriptome_kal.idx" "${extraPaths[1]}"

# kallisto quant

echo "Kallisto quantifications..."
echo "Outputs will be sent to log"

echo kallisto_quant.py "${outDir}/" "${outDir}/batch_trim/processed/" "$cores"

# Salmon index

echo "Creating the Salmon index..."

echo salmon index -t "${extraPaths[1]}" -i "aux_files/GRCh38transcriptome_sal.idx"

# Salmon quant

echo "Salmon quantifications..."
echo "Outputs will be sent to log"

echo salmon_quant.py "${outDir}/" "${outDir}/batch_trim/processed/" "$cores"

# Salmon prep with Wasabi

echo "Wasabi-ing that Salmon"

echo Rscript wasabi.R "${outDir}/salmon/processed"

# Sleuth analysis for kallisto

echo sleuth_wrapper.sh "in" "out"

# Sleuth analysis for Salmon

echo sleuth_wrapper.sh "in" "out"

# Create heatmaps
