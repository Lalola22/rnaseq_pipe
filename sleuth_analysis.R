#!/usr/bin/env Rscript
# Base R script for sleuth analysis
# Sam Lee
# Created 04/06/17

# Usage

# sleuth_basic.R "path/to/kallisto/files/" "treatment" "condition 2" path/to/outputs/"
# assumes three samples per condition, only two conditions
# also outs diagnostic and analytic plots of the seluth data

# Setup -------------------------------------------------------------------

# source required packages
### need to go through all the effort of re sorting the packages
library(tidyverse)
library(sleuth)
library(stringr)
library(biomaRt)



## change this back to 16 later.
## commandline arguments
args <- commandArgs(trailingOnly = TRUE)
# test if there is at least one argument: if not, return an error
if (length(args) < 6) {
  stop("Six arguments must be supplied.", call. = FALSE)
}
basedir    <- args[1]
outdir     <- args[2]
max_cores  <- as.numeric(args[3])
treatment <- args[4]
control <- args[5]
replicates <- as.numeric(args[6])

#  Set number of cores available for sleuth
options(mc.cores = max_cores)


# complete output directory name ------------------------------------------

outdir <- file.path(outdir, paste(treatment, control, sep = "_"))


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


# Fetch gene names --------------------------------------------------------


mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                         dataset = "hsapiens_gene_ensembl",
                         host    = "ensembl.org")

t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                     "external_gene_name"), mart = mart)

t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
                     ens_gene = ensembl_gene_id, ext_gene = external_gene_name)

so <- sleuth_prep(sample_table, ~ condition, target_mapping = t2g)
## look into this properly
# create sleuth object ----------------------------------------------------


so <- sleuth_fit(so)
so <- sleuth_fit(so, ~1, "reduced")

so <- sleuth_wt(so, paste("condition", control, sep = ""))
print("Wald test completed...", quote = FALSE)

# flip the beta ----------------------------------------------------

results_table <- sleuth_results(so, paste("condition", control, sep = "")) %>%
    mutate(b = -1 * b)



# make and save results table ---------------------------------------------
dir.create(outdir, recursive = TRUE) # will show a warning if this exists...

write_csv(results_table, file.path(outdir, "sleuth_results.csv"))

# Sleuth plots --------------------------------------------------------
print("making graphs now...", quote = FALSE)


# MA plot
pdf(file.path(outdir, "MA_plot.pdf"), width = 7, height = 4)
plot_ma(so, paste("condition", control, sep = ""))
dev.off()

# qq plot
pdf(file.path(outdir, "qq_plot.pdf"), width = 7, height = 4)
plot_qq(so, paste("condition", control, sep = ""))
dev.off()

# PCA
pdf(file.path(outdir, "pca_plot.pdf"), width = 7, height = 4)
plot_pca(so, text_labels = TRUE)
dev.off()

# PC variance
pdf(file.path(outdir, "pc_variance_plot.pdf"), width = 7, height = 4)
plot_pc_variance(so)
dev.off()

# Mean variance
pdf(file.path(outdir, "mean_variance_plot.pdf"), width = 7, height = 4)
plot_mean_var(so)
dev.off()

# Sample heatmap
pdf(file.path(outdir, "sample_heatmap_plot.pdf"), width = 7, height = 4)
plot_sample_heatmap(so)
dev.off()

# volcano plot
pdf(file.path(outdir, "volcano_plot.pdf"), width = 7, height = 4)
plot_volcano(so, paste("condition", control, sep = ""))
dev.off()

# top-20 heatmap
# First order the top 20 transcripts based on abs(b) and then extract transcript names
tab <- results_table[order(abs((results_table[abs(results_table$b) > 1.5, ])$b),
                           decreasing = TRUE), ]
## FIXME This is fine for now, but check code actually gets the top expression
pdf(file.path(outdir, "transcript_heatmap.pdf"), width = 7, height = 10)
plot_transcript_heatmap(so, tab$target_id[1:20], trans = "log2")
dev.off()
