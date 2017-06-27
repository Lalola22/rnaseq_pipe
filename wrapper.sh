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
    mkdir $outDir
fi

# Initialise bash array with raw-read paths

fastqArray=($(find $rawDir -type f -name "*.fastq.gz"))

# Read the two auxillary file into bash arrays

cd aux_files

oldIFS="$IFS"
IFS=$'\n' extraPaths=($(<"extra_paths.txt"))
IFS="$oldIFS"

oldIFS="$IFS"
IFS=$'\n' samples=($(<"sample_table.txt"))
IFS="$oldIFS"

cd -

# Set paths

trimmomaticPath=${extraPaths[0]}
GRCh38trans=${extraPaths[1]}

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

# -- kallisto index

# echo "Creating the kallisto index..."
#
# mkdir -p "${outDir}/reference_files"
#
# kallisto index -i "${outDir}/reference_files/GRCh38transcriptome_kal.idx" \
# "$GRCh38trans"

# -- kallisto quant
#
# echo "Kallisto quantifications..."

# python3 kallisto_quant.py "${outDir}/" "${outDir}/batch_trim/" "$cores"

# -- Salmon index

echo "Creating the Salmon index..."

salmon index -t "$GRCh38trans" -i \
"${outDir}/reference_files/GRCh38transcriptome_sal.idx"

# -- Salmon quant

echo "Salmon quantifications..."

python3 salmon_quant.py "${outDir}/" "${outDir}/batch_trim/" "$cores"

# # -- Salmon prep with Wasabi
#
# echo "Wasabi-ing that Salmon"
#
# Rscript wasabi.R "${outDir}/salmon"
#
# # -- Sleuth analysis for kallisto
#
# echo "Sleuthing around kallisto..."
#
# # make array with comparisons as each element
#
# ((n_elements=${#samples[@]}, max_index=n_elements - 1))
#
# # i starts at one to avoid header row
# for ((i = 1; i <= max_index; i++)); do
#   echo "Running sleuth for: " "${samples[i]}"
#   Rscript sleuth_analysis.R "${outDir}/kallisto" "${outDir}/sleuth_kallisto" \
#   "$cores" ${samples[i]}
# done
#
# # -- Sleuth analysis for Salmon
#
# echo "Sleuthing around Salmon..."
#
# for ((i = 1; i <= max_index; i++)); do
#   echo "Running sleuth for: " "${samples[i]}"
#   Rscript sleuth_analysis.R "${outDir}/salmon" "${outDir}/sleuth_salmon" \
#   "$cores" ${samples[i]}
# done

# -- Prep for DESeq2

Rscript create_tx2g.R "${outDir}/reference_files" "$cores"


# -- DESeq2 for salmon

# make array with comparisons as each element

((n_elements=${#samples[@]}, max_index=n_elements - 1))

echo "DeSeq2 analysis of salmon quants happening..."


for ((i = 1; i <= max_index; i++)); do
  echo "Running deseq2 for: " "${samples[i]}"
  Rscript dseq2_analysis.R "/$outDir" \
  "salmon" "$cores" ${samples[i]};
 done


echo "Finished at" $(date)

echo "Results are in " ${outDir}

echo "Check the FastQC results to ensure read quality was appropriate."
