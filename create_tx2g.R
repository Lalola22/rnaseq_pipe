# 09/06/2017
# Sam Lee
# Creation of the tx 2 gene table to DESeq2 mappings


# Setup -------------------------------------------------------------------

library(tidyverse, quietly = TRUE)
library(BiocParallel, quietly = TRUE)
library(ensembldb, quietly = TRUE)
library(EnsDb.Hsapiens.v79, quietly = TRUE)

args <- commandArgs(trailingOnly = TRUE)

outdir <- args[1]
cores <- args[2]

# Code body -------------------------------------------------------------------

txdb <- EnsDb.Hsapiens.v79
k <- keys(txdb, keytype = "GENENAME")
df <- select(txdb, keys = k, keytype = "GENENAME", columns = "TXID")
tx2gene <- df[, 2:1]  # tx ID, then gene ID
