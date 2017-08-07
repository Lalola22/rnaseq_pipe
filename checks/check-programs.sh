#! /bin/bash
set -e
set -u
set -o pipefail

# This is just a very basic check for software
# run from top of repository
# If one is not available will fail with an error
# pass the directory to trimmomatic as the first arguements

fastqc

java jar ${1}/trimmomatic-0.36.jar

salmon cite

which R

which python3

echo "All required software are installed :~)"
