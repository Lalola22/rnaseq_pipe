"""
Created 03/06/17
For executation of the Salmon quantification step
To be run with three arguements
    * basedir - top level output directory
    * input directory - contains folders with .fastq.gz files
    * max_threads - how many threads to allocate to salmon
Returns salmon quantifications and associated log files to a directory
within the top level output dir.
"""

# --- packages

import time
import os
import sys
import subprocess


# --- variables using sys.argv

basedir = sys.argv[1]
inputdirectory = sys.argv[2]
max_threads = sys.argv[3]
processed = basedir + "salmon/processed/"
log = basedir + "salmon/log/"


# --- functions

def salmon_call(read1):
    dividing = read1.split(".")
    basename = dividing[0].replace("1P", "")
    read2 = read1.replace("1P", "2P")
    subprocess.call(" echo salmon quant -i aux_files/GRCh38transcriptome_sal.idx " +
        "-l A --numBootstraps 100 -1 " + inputdirectory + read1 + " -2 " +
        inputdirectory + read2 + " -p " + max_threads + " -o " + basename +
        " >" + log + time.strftime("%Y%m%d") + "_" + basename + ".txt" + " 2>&1",
        shell=True)

# --- main call

if __name__ == "__main__":
    # --- check dirs and create if neccessary
    if not os.path.exists(processed):
        os.makedirs(processed)
    if not os.path.exists(log):
        os.makedirs(log)
    # --- create list of read1 pair file names
    read_list = []
    for fname in os.listdir(inputdirectory):
        if "1P" in fname:
            read_list.append(fname)
    # --- call kallisto_call on each read pair in parallel
    for read in read_list:
        salmon_call(read)
