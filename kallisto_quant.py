"""
created 09/05/17
# Purpose

- Parallel execution of "Kallisto" for alignment and quant
    The bootstrapping is being run in parallel as opposed to each kallisto call
- reference transcriptome should be placed within the kallisto folder

basedir is the top level output dir
inputdirectory should contain all folders with .fastq.gz reads to be processed
processed is where the trimmed read files will go
log is self explanatory
program is directory for flexbar and the .idx file
max_threads is the maximum number of threads to allocate to the process



An example pair of files is:

    25uM_1_R1_trimmed_1P.fastq.gz
    25uM_1_R1_trimmed_2P.fastq.gz

Outputs kallisto files for each read pair and
associated log files in a nested directory
"""
# --- packages

import os
import sys
from subprocess import call


# --- variables using sys.argv

basedir = sys.argv[1]
inputdirectory = sys.argv[2]
max_threads = sys.argv[3]
processed = basedir + "kallisto/"


# --- functions

def kallisto_call(read1):
    """
    l is the lock object
    read1 is the forward read
    calls kallisto quant for the read pair specified by the arguements
    Rewrite this to be more specific for a single read pair
    so it can be parallelised
    also review how to actually do this... current way does not seem to.
    """
    dividing = read1.split(".")
    basename = dividing[0].replace("_1P", "")
    read2 = read1.replace("1P", "2P")
    call("kallisto quant -i aux_files/GRCh38transcriptome_kal.idx -t " +
        max_threads + " -o " + processed + basename + " -b 100 " +
        inputdirectory + read1 + " " + inputdirectory + read2, shell=True)

# --- __main__ call


if __name__ == "__main__":
    # --- check dirs and create if neccessary
    if not os.path.exists(processed):
        os.makedirs(processed)
    # --- create list of read1 pair file names
    read_list = []
    for fname in os.listdir(inputdirectory):
        if "1P" in fname:
            read_list.append(fname)
    # --- call kallisto_call on each read pair in parallel
    for read in read_list:
        kallisto_call(read)
