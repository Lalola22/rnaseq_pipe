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

Outputs kallisto files for each read pair and associated log files in a nested directory
"""
# --- packages

import time
import os
import sys
import subprocess
from multiprocessing import Process, Lock

# --- variables using sys.argv

basedir         = sys.argv[1]
inputdirectory  = sys.argv[2]
program         = sys.argv[3]
processed       = basedir + "kallisto/processed/"
log             = basedir + "kallisto/log/"
max_threads     = (os.cpu_count() * 2)

# --- functions

def kallisto_call(read1, threads):
    """
    l is the lock object
    read1 is the forward read
    calls kallisto quant for the read pair specified by the arguements
    Rewrite this to be more specific for a single read pair
    so it can be parallelised
    also review how to actually do this... current way does not seem to.
    """
    dividing = read1.split(".")
    read2 = read1.replace("1P", "2P")
    subprocess.call(program + "kallisto quant -i " + program + "GRCh38transcriptome.idx -t " + threads + " -o \\" +
    processed + dividing[0] + " -b 100 " + inputdirectory + read1 + " " + inputdirectory + read2 + " >" + log + time.strftime("%Y%m%d-%H%M%S") + "_" + dividing[0] + ".txt" + " 2>&1", shell=True)


    #
    #     ## need to work out what the direcotry structure will be before I can write this code.
    #     ## ie is it more than main dir --> subdirs (if so write function to extract subdirs)
    #     fulldir = inputdirectory + subdir + "/"
    #     ## add a for read1 loop to idnetify read1's if needed? based on dir structure
    #     if os.path.isfile(fulldir + read1 + fastq.gz):
    #         r1 = read1 + fastq.gz
    #         r2 = r1.replace("1","2") #FIXME check what actually will have to be replaced from the trimmed names
    #         output = processed + r1 #FIXME check this is correct with what docs want
    #         subprocess.call("echo " + program + "kallisto quant -i " + program + "GRCh38transcriptome.idx -o " output + " \ " +
    #         "-b 100 " + r1 + " " + r2 + " >" + log + target + "_kallisto_quant.txt" +" 2>&1", shell=True)
    # finally:
    #     l.release()

# --- __main__ call

## For parallel
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
        kallisto_call(read, str(max_threads))
