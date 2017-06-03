#! /bin/python3

"""
# Purpose

To batch process raw fastq files for adaptor trimming and QC trimming
with trimmomatic.
basedir is the top level output dir
inputdirectory should contain all folders with .fastq.gz reads to be processed
processed is where the trimmed read files will go
log is self explanatory
program is directory for trimmomatic jar file

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
import subprocess

from multiprocessing import Process, Lock

# --- variables using sys.argv

basedir = sys.argv[1]
inputdirectory = sys.argv[2]
program = sys.argv[3]
processed = basedir + "batch_trim/processed/"
log = basedir + "batch_trim/log/"
aux = "aux_files/"

# --- code body

def call_trimmomatic_par(l, subdir):
    """
    call wrapper for trimmomatic for trimming of adaptors sequences
    designed to be able to be called for parallel processing
    call in a loop (or another function) to loop through each dir for different classes
    """
    l.acquire()
    try:
        fulldir = inputdirectory + subdir + "/"  # construct full path of each class
        try:
            os.path.isdir(fulldir)
        except Exception as e:
            print("{} is not a directory".format(fulldir))
        for fileR1 in os.listdir(fulldir):
            dividing = fileR1.split(".")
            if "R1" in fileR1:  # select each of the forward read files to loop through
                fileR2 = fileR1.replace("R1", "R2")
                if os.path.isfile(fulldir + fileR2):
                    dividing1 = fileR2.split(".")
                    basename = dividing[0].replace("_R1", "")  # This is the forward read sans extensions
                    subprocess.call("echo java -jar " + program +
                    "trimmomatic-0.36.jar PE -phred33 " +
                    fulldir + fileR1 + " " + fulldir + fileR2 + " -baseout " +
                    processed + basename +"_trimmed.fastq.gz ILLUMINACLIP:" +
                    aux + "/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:35" +
                    # " > " + log + "_" + basename + "_trimmomatic.txt" +" 2>&1",
                    shell=True)
    finally:
        l.release()

# --- __main__ call


if __name__ == "__main__":
    lock = Lock()
    # --- check dors and create if neccessary
    if not os.path.exists(processed):
        os.makedirs(processed)
    if not os.path.exists(log):
        os.makedirs(log)
    # --- loop for call to trimmomatic
    for sub in os.listdir(inputdirectory):
        if os.path.isdir(inputdirectory + sub):
            Process(target=call_trimmomatic_par, args=(lock, sub)).start()
