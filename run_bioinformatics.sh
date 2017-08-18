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
    mkdir -p $outDir
fi

# Initialise bash array with raw-read paths

fastqArray=($(find $rawDir -type f -name "*.fastq.gz"))

# Read the two auxillary file into bash arrays

cd aux_files

oldIFS="$IFS"
IFS=$'\n' extraPaths=($(<"extra_paths.txt"))
IFS="$oldIFS"

# oldIFS="$IFS"
# IFS=$'\n' samples=($(<"sample_table.txt"))
# IFS="$oldIFS"

cd -

# Set paths

trimmomaticPath=${extraPaths[0]}
transcriptome=${extraPaths[1]}

export PATH=${PWD}:$PATH

# -- Pre-trim QC

echo "Inital read QC..."

mkdir -p "${outDir}/fastqc/raw"

fastqc -o "${outDir}/fastqc/raw" -t $cores ${fastqArray[@]]}


# -- Adaptor trimming

echo "Trimming adaptors..."

python3 batch_trim.py "${outDir}/" "${rawDir}/" "${trimmomaticPath}/" $cores


# -- Trimmed read QC

echo "Post-trimming read QC"

# Initialise bash array for trimmed + paired reads

fastqArrayTrim=($(find "${outDir}/batch_trim" -type f -name "*P.fastq.gz"))

mkdir -p "${outDir}/fastqc/trimmed"

fastqc -o "${outDir}/fastqc/trimmed" -t $cores ${fastqArrayTrim[@]]}

# -- Salmon index
#
# echo "Creating the Salmon index..."
#
# salmon index -t "$GRCh38trans" -i \
# "${outDir}/reference_files/GRCh38transcriptome_sal.idx"
#
# # -- Salmon quant
#
# echo "Salmon quantifications..."
#
# python3 salmon_quant.py "${outDir}/" "${outDir}/batch_trim/" "$cores"

# -- Kallisto index

echo "Creating the Kallisto index for " $(basename $transcriptome)

kallisto index -i "${outDir}/transcriptome_kallisto.idx" \
"$transcriptome"

# -- Kallisto quants

echo "Running Kallisto to quantify reads"

python3 kallisto_quant.py "${outDir}/" "${outDir}/batch_trim/" "$cores"

# -- Prep for DESeq2

Rscript create_tx2g.R "${outDir}/reference_files" "$cores"


# -- DESeq2 for salmon

# make array with comparisons as each element

# ((n_elements=${#samples[@]}, max_index=n_elements - 1))

echo "DESeq2 analysis for differential expression at the \
gene and transcript level"

Rscript deseq2_analysis.R "$outDir" "aux_files" "$cores"


# for ((i = 1; i <= max_index; i++)); do
#   echo "Running deseq2 for: " "${samples[i]}"
#   Rscript dseq2_analysis.R "/$outDir" \
#   "salmon" "$cores" ${samples[i]};
#  done


echo "Finished at" $(date)

echo "Results are in " ${outDir}

echo "Check the FastQC results to ensure read quality was appropriate."
