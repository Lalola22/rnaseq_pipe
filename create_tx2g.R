# 09/06/2017
# Sam Lee
# Creation of the tx 2 gene table to DESeq2 mappings


# Setup -------------------------------------------------------------------

library(tidyverse, quietly = TRUE, verbose = FALSE, warn.conflicts = FALSE)
library(BiocParallel, quietly = TRUE, verbose = FALSE, warn.conflicts = FALSE)
library(ensembldb, quietly = TRUE, verbose = FALSE, warn.conflicts = FALSE)
library(EnsDb.Hsapiens.v86, quietly = TRUE,
    verbose = FALSE, warn.conflicts = FALSE)

args <- commandArgs(trailingOnly = TRUE)

outdir <- args[1]
cores <- args[2]

register(MulticoreParam(cores))

# Code body -------------------------------------------------------------------

print("Forming gene to transcript mapping...", quote = FALSE)

txdb <- EnsDb.Hsapiens.v86
k <- keys(txdb, keytype = "GENENAME")
df <- select(txdb, keys = k, keytype = "GENENAME", columns = "TXID")
tx2gene <- df[, 2:1]  # tx ID, then gene ID

print("Saving the gene to transcript mapping...", quote = FALSE)

write_csv(tx2gene, path = file.path(outdir, "tx2g_table.csv"))
