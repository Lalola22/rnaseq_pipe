#! /bin/python3

"""
# Purpose

To batch process raw fastq files for adaptor trimming and QC trimming
with trimmomatic.
basedir is the top level output dir
inputdirectory should contain all folders with .fastq.gz reads to be processed
processed is where the trimmed read files will go
log is self explanatory
trim is directory for trimmomatic jar file

An example pair of files is:

    Untreated_1_R1.fastq.gz
    Untreated_1_R2.fastq.gz

Calls trimmomatic with the following parameters
- PE
-threads automatic ## check this
-phred 33
- ILLUMINACLIP:TruSeq3-PE.fa:2:30:10
- LEADING:3
- TRAILING:3
- SLIDINGWINDOW:4:15
- MINLEN:35

"""
# --- packages
import os
import sys
from subprocess import call

from multiprocessing import Process, Lock
from joblib import Parallel, delayed

# --- variables using sys.argv

basedir = sys.argv[1]
inputdirectory = sys.argv[2]
trim = sys.argv[3]
max_threads = int(sys.argv[4])
processed = basedir + "batch_trim/"
aux = "aux_files/"

# --- code body


def call_trimmomatic_par(subdir):
    """
    call wrapper for trimmomatic for trimming of adaptors sequences
    designed to be able to be called for parallel processing
    call in a loop (or another function) to loop through each
    dir for different classes
    """
    fulldir = inputdirectory + subdir + "/"  # construct full path
    try:
        os.path.isdir(fulldir)
    except Exception as e:
        print("{} is not a directory".format(fulldir))
    for fileR1 in os.listdir(fulldir):
        dividing = fileR1.split(".")
        if "R1" in fileR1:  # select forward read files to loop through
            fileR2 = fileR1.replace("R1", "R2")
            if os.path.isfile(fulldir + fileR2):
                basename = dividing[0].replace("_R1", "")
                # This is the forward read sans extensions
                call("java -jar " + trim +
                     "trimmomatic-0.36.jar PE -phred33 " + fulldir + fileR1 +
                     " " + fulldir + fileR2 + " -baseout " + processed +
                     basename + ".fastq.gz ILLUMINACLIP:" + aux +
                     "TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 " +
                     "SLIDINGWINDOW:4:15 MINLEN:35",
                     shell=True)


# --- __main__ call


if __name__ == "__main__":
    lock = Lock()
    # --- check dirs and create if neccessary
    if not os.path.exists(processed):
        os.makedirs(processed)
    #  call the func in parallel, as many cores as indicated
    Parallel(n_jobs=max_threads)(delayed(call_trimmomatic_par)(sub) for sub in os.listdir(inputdirectory) if os.path.isdir(inputdirectory + sub))
