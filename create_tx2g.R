# 09/06/2017
# Sam Lee
# Creation of the tx 2 gene table to DESeq2 mappings
# Two arguments required
# 1. dir to save map file to
# 2. species to make the map for, there are three options
#       a. Human
#       b. Mouse
#       c. Rat

# Setup -------------------------------------------------------------------

library(tidyverse, quietly = TRUE, verbose = FALSE, warn.conflicts = FALSE)
library(ensembldb, quietly = TRUE, verbose = FALSE, warn.conflicts = FALSE)


args <- commandArgs(trailingOnly = TRUE)

outdir <- args[1]
species <- toLower(args[2])

# check that the species is correct
speciesCheck <- species %in% c("human", "mouse", "rat")
if (! dircheck == TRUE){
  stop('Species must be one of "human", "mouse", or "rat"')
}

# assign species to  correct Db
txdb <- if (species == "human") {
   library(EnsDb.Hsapiens.v86, quietly = TRUE, verbose = FALSE, warn.conflicts = FALSE)
} else {if (species == "mouse") {
   library(EnsDb.Mmusculus.v79, quietly = TRUE, verbose = FALSE, warn.conflicts = FALSE)
} else {
   library(EnsDb.Rnorvegicus.v79, quietly = TRUE, verbose = FALSE, warn.conflicts = FALSE)
}
   
}

# Code body -------------------------------------------------------------------

print("Forming gene to transcript mapping...", quote = FALSE)

txdb <- EnsDb.Hsapiens.v86
k <- keys(txdb, keytype = "GENENAME")
df <- select(txdb, keys = k, keytype = "GENENAME", columns = "TXID")
tx2gene <- df[, 2:1]  # tx ID, then gene ID

print(sprintf("Saving the gene to transcript mapping for %s", species), quote = FALSE)

write_csv(tx2gene, path = file.path(outdir, "tx2g_table.csv"))
