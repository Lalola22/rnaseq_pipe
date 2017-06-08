# Created 08/06/17
# Sam Lee
# For Deseq2 analysis of kallisto results
# takes five inputs
# location of kallisto files
# location of results to be saved
# number of cores to use.
# treatment condition
# control condition

# Setup -------------------------------------------------------------------

library(tidyverse)
library(tximport)
library(DESeq2)
library(ensembldb)
library(EnsDb.Hsapiens.v79)

args <- commandArgs(trailingOnly = TRUE)
# test if there is at least one argument: if not, return an error
if (length(args) < 6) {
     stop("Six arguments must be supplied.", call. = FALSE)
}

basedir <- args[1]
outdir <- args[2]
cores <- args[3]
treatment <- args[4]
control <- args[5]

# Test values

basedir <- "/home/slee/outputs/june_05/kallisto"

# whether the results are for kallisto or sleuth data
type <- basename(basedir)

# sample table construction -----------------------------------------------

sample_table <- as.data.frame(matrix("", nrow = replicates * 2, ncol = 3),
                                     stringsAsFactors = FALSE)
colnames(sample_table) <- c("sample", "condition", "path")
sample_table$condition <- (c(rep(treatment, replicates),
                            rep(control, replicates)))

sample_table <- sample_table %>%
  mutate(sample = c(1:replicates,
                    1:replicates)) %>%
  mutate(sample = paste(condition, sample, sep = "_")) %>%
  mutate(path = file.path(basedir, sample))



## print to display as a check while program runs
print("Sample table:", quote = FALSE)
print(sample_table)

## check that all dirs in the table exist and exit if not
dircheck <- all(dir.exists(sample_table$path))
if (! dircheck == TRUE){
  stop("All the quantification directories for the conditions do not exist.")
}

# tximport -------------------------------------------------------------------

## Create the transcript_id to gene mappings


txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
k <- keys(txdb, keytype = "GENEID")
df <- select(txdb, keys = k, keytype = "gene_id", columns = "tx_name")
tx2gene <- df[, 2:1]  # tx ID, then gene ID


## Import the quantification data

txi <- tximport(files, type = type, tx2gene = tx2gene)

ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ condition)


# DESeq2 -------------------------------------------------------------------
