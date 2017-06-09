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
    mkdir "$outDir"
    mkdir "$outDir/reference_files"
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


# -- Prep for DESeq2


Rscript create_tx2g.R "${outDir}/reference_files" "$cores"


# -- DESeq2 for salmon

# make array with comparisons as each element

((n_elements=${#samples[@]}, max_index=n_elements - 1))

echo "DeSeq2 analysis of salmon quants happening..."

echo "the quant dir is still hard-coded due to not running the full pipeline"

for ((i = 1; i <= max_index; i++)); do
  echo "Running deseq2 for: " "${samples[i]}"
  echo Rscript dseq2_analysis.R "/home/slee/outputs/test/june_09" \
  "salmon" "$cores" ${samples[i]};
 done


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
