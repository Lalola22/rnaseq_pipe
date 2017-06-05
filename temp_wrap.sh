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

fastqArray=($(find $rawDir -type f -name "*.fastq.gz"))

oldIFS="$IFS"
IFS=$'\n' extraPaths=($(<"extra_paths.txt"))
IFS="$oldIFS"

oldIFS="$IFS"
IFS=$'\n' samples=($(<"sleuth_samples.txt"))
IFS="$oldIFS"

trimmomaticPath=${extraPaths[0]}
GRCh38trans=${extraPaths[1]}

export PATH=${PWD}:$PATH


# Salmon quant

echo "Salmon quantifications..."

python3 salmon_quant.py "${outDir}/" "${outDir}/batch_trim/" "$cores"

# Salmon prep with Wasabi

echo "Wasabi-ing that Salmon"

Rscript wasabi.R "${outDir}/salmon"

# Sleuth analysis for kallisto

echo "Sleuthing around kallisto..."

# make array with comparisons as each element
((n_elements=${#samples[@]}, max_index=n_elements - 1))

# i starts at one to avoid header
for ((i = 1; i <= max_index; i++)); do
  echo "Running sleuth for: " "${samples[i]}"
  Rscript sleuth_analysis.R "${outDir}/kallisto" "${outDir}/sleuth_kallisto" \
  "$cores" ${samples[i]}
done

# Sleuth analysis for Salmon

echo "Sleuthing around Salmon..."

for ((i = 1; i <= max_index; i++)); do
  echo "Running sleuth for: " "${samples[i]}"
  Rscript sleuth_analysis.R "${outDir}/kallisto" "${outDir}/sleuth_salmon" \
  "$cores" ${samples[i]}
done


echo "Finished at" $(date)
