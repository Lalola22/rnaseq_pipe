#! /bin/bash
set -e
set -u
set -o pipefail

# This is just a very basic check for software
# run from top of repository
# If one is not available will fail with an error


# get the path to trimmomatic
oldIFS="$IFS"
IFS=$'\n' extraPaths=($(<"aux_files/extra_paths.txt"))
IFS="$oldIFS"
trimmomaticPath=${extraPaths[0]}

fastqc -version

# java jar ${trimmomaticPath}/trimmomatic-0.36.jar

echo Salmon:
salmon -v

which R

which python3

echo "All required software are installed :~)"
