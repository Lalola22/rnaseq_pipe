# /bin/bash
set -e
set -u
set -o pipefail

# usage: scriptname.sh "/path/to/fastq/dirs" "/path/to/output/dir"

rawDIR="$1"
outdir="2"
