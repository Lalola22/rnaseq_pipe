#! python3
"""
Created 04/05/17
code taken from https://www.biostars.org/p/237119/

# usage

`python3 batch_trim.py "/path/to/basedir" "path/to/inputdirectory" "path/to/program"`

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

--- ToDo

Check how to run this in parallel, with the loop??
Check number of threads to assign

--- In progress

Changing directories to use bash arguements...

"""
# --- packages
import os
import sys
import subprocess

from multiprocessing import Process, Lock

# --- variables - for testing
# basedir="/Users/slee/programming/outputs/"
# inputdirectory = "/Volumes/sbs_microbio_01/MMc_group/Sam_L/rna_seq_489/test_subdir/"
# processed = basedir + "batch_trim/processed/"
# log= basedir + "batch_trim/log/parallel/"
# program = "/Users/slee/Downloads/Trimmomatic-0.36/"

# --- variables using sys.argv

basedir         = sys.argv[1]
inputdirectory  = sys.argv[2]
program         = sys.argv[3]
processed       = basedir + "batch_trim/processed/"
log             = basedir + "batch_trim/log/"


# --- code body

# --- Function version

def call_trimmomatic(subdir):
    """
    call wrapper for trimmomatic for trimming of adaptors sequences
    designed to be able to be called for parallel processing
    call in a loop (or another function) to loop through each dir for different classes
    """
    fulldir = inputdirectory + subdir + "/" # construct full path of each class
    try:
        os.path.isdir(fulldir)
    except Exception as e:
        print("{} is not a directory".format(fulldir))
    for fileR1 in os.listdir(fulldir):
        dividing = fileR1.split(".")
        if "R1" in fileR1: # select each of the forward read files to loop through
            fileR2 = fileR1.replace("R1", "R2")
            if os.path.isfile(fulldir + fileR2):
                dividing1 = fileR2.split(".")
                basename = dividing[0] # This is the forward read sans extensions
                subprocess.call("echo java -jar " + program + "trimmomatic-0.36.jar PE -phred33 \ " +
                fulldir + fileR1 + " " + fulldir + fileR2 + " -baseout " + processed + basename +"_trimmed.fastq.gz \ " +
                "ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:35" +" >" + log + basename + "_trimmomatic.txt" + " 2>&1", shell=True)

def call_trimmomatic_par(l, subdir):
    """
    call wrapper for trimmomatic for trimming of adaptors sequences
    designed to be able to be called for parallel processing
    call in a loop (or another function) to loop through each dir for different classes
    """
    l.acquire()
    try:
        fulldir = inputdirectory + subdir + "/" # construct full path of each class
        try:
            os.path.isdir(fulldir)
        except Exception as e:
            print("{} is not a directory".format(fulldir))
        for fileR1 in os.listdir(fulldir):
            dividing = fileR1.split(".")
            if "R1" in fileR1: # select each of the forward read files to loop through
                fileR2 = fileR1.replace("R1", "R2")
                if os.path.isfile(fulldir + fileR2):
                    dividing1 = fileR2.split(".")
                    basename = dividing[0] # This is the forward read sans extensions
                    subprocess.call("java -jar " + program + "trimmomatic-0.36.jar PE -phred33 \\" +
                    fulldir + fileR1 + " " + fulldir + fileR2 + " -baseout \\" +
                    processed + basename +"_trimmed.fastq.gz ILLUMINACLIP:" + program + "adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:35" +" >" + log + sub + "_" + basename + "_trimmomatic.txt" + " 2>&1", shell=True)
    finally:
        l.release()
# #  Rewriting to allow the file list process to enter subdirectories
# #
# for filedict in os.walk(inputdirectory):
#     for fileR1 in filedict[2]:
#         dividing = fileR1.split(".")
#         if ("R1" in fileR1) :
#             # print(fileR1)
#             fileR2 = fileR1.replace('R1', 'R2')
#             if os.path.isfile(filedict[0] + "/" + fileR2) :
#                 print("\n", fileR1, "\n", fileR2)
#                 dividing1 = fileR2.split(".")
#                 log1 = dividing[0]
#                 output1 = dividing[0]
#                 output2 = dividing1[0]
#                 print("echo java -jar " + program + "trimmomatic-0.36.jar PE -threads 12 -phred33 \ "
#                 + inputdirectory + fileR1 + " " + inputdirectory + fileR2 + " -baseout " + processed + output1 +"_trimmed.fastq.gz \ "
#                 + "ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:35" +" >" + log + log1 + "_trimmomatic.txt" + " 2>&1")#, shell=True) FIXME


# --- __main__ call

## For parallel
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

# # for non-parallel
# if __name__ == "__main__":
#     # --- check dors and create if neccessary
#     if not os.path.exists(processed):
#         os.makedirs(processed)
#     if not os.path.exists(log):
#         os.makedirs(log)
#     # --- loop for call to trimmomatic
#     for sub in os.listdir(inputdirectory):
#         if os.path.isdir(inputdirectory + sub):
#             call_trimmomatic(sub)
