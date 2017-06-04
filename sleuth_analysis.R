#!/usr/bin/env Rscript
# Base R script for sleuth analysis
# Sam Lee
# 11/05/17

# Usage

# sleuth_basic.R "path/to/kallisto/files/" "condition1" "condition 2" path/to/outputs/"
# assumes three samples per condition, only two conditions
# also outs diagnostic and analytic plots of the seluth data

# Setup -------------------------------------------------------------------

# source required packages
### need to go through all the effort of re sorting the packages
source("https://bioconductor.org/biocLite.R")
library(sleuth)
library(stringr)


#  Set number of cores available for sleuth
options(mc.cores = 16L) # set number of cores sleuth will use
## change this back to 16 later.
## commandline arguments
args <- commandArgs(trailingOnly = TRUE)
# test if there is at least one argument: if not, return an error
if (length(args) < 3) {
  stop("At least three arguments must be supplied (input file).n", call. = FALSE)
} else if (length(args) == 3) {
  # default output dir
  args[4] <- "~/outputs/"
}
basedir    <- args[1]
condition1 <- args[2]
condition2 <- args[3]
outdir     <- args[4]

# test parameters
# stop("Test parameters being used!!")
# basedir    <- "/home/slee/outputs/kallisto/processed/"
# condition1 <- "DMSO"
# condition2 <- "Untreated"
# outdir     <- "/home/slee/outputs/sleuth/"


# complete output directory name ------------------------------------------

outdir <- file.path(outdir, paste(condition1, condition2, sep = "_"))

# create object that has sample paths ------------------
# kal_paths <- list.dirs(basedir, full.names = TRUE)
# kal_paths <- c(kal_paths[grep(x= kal_paths, pattern = condition2)] ,kal_paths[grep(x= kal_paths, pattern = condition1)])




# sample table construction -----------------------------------------------

sample_table <- as.data.frame(matrix("", nrow = 6, ncol = 3),
                                     stringsAsFactors = FALSE)
colnames(sample_table) <- c("sample", "condition", "path")
sample_table$condition <- (c(rep(condition1, 3), rep(condition2, 3)))

c <-  0
for (n in c(1:3, 1:3)){
  c <- c + 1
  sample_table$sample[c] <- paste(sample_table$condition[c], "_", n, sep = "" )
  sample_table$path[c]   <- file.path(basedir, sample_table$sample[c])
    }
rm(c)
## rename the sample names so they're more readable


## print to display as a check while program runs
print("Sample table:", quote = FALSE)
print(sample_table)

## check that all dirs in the table exist and exit if not
dircheck <- all(dir.exists(sample_table$path))
if (! dircheck == TRUE){
  stop("All the kallisto output directories for the conditions do not exist.")
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

so <- sleuth_wt(so, paste("condition", condition2, sep = ""))
print("Wald test completed...", quote = FALSE)

# make and save results table ---------------------------------------------
dir.create(outdir, recursive = TRUE) # will show a warning if this exists...
results_table <- sleuth_results(so, paste("condition", condition2, sep = ""))
write.csv(file = file.path(outdir, "sleuth_results.csv"), x = results_table)


# Sleuth plots --------------------------------------------------------
print("making graphs now...", quote = FALSE)


# MA plot
pdf(file.path(outdir, "MA_plot.pdf"), width = 7, height = 4)
plot_ma(so, paste("condition",condition2,sep = ""))
dev.off()

# qq plot
pdf(file.path(outdir, "qq_plot.pdf"), width = 7, height = 4)
plot_qq(so, paste("condition", condition2, sep = ""))
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
plot_volcano(so, paste("condition", condition2, sep = ""))
dev.off()

# top-20 heatmap
# First order the top 20 transcripts based on abs(b) and then extract transcript names
tab <- results_table[order(abs((results_table[abs(results_table$b) > 1.5, ])$b),
                           decreasing = TRUE), ]
## FIXME This is fine for now, but check code actually gets the top expression
pdf(file.path(outdir, "transcript_heatmap.pdf"), width = 7, height = 10)
plot_transcript_heatmap(so, tab$target_id[1:20], trans = "log2")
dev.off()
