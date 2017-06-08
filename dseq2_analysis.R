# Created 08/06/17
# Sam Lee
# For Deseq2 analysis of kallisto results
# takes six inputs
# location of kallisto files
# location of results to be saved
# number of cores to use.
# treatment condition
# control condition
# number of bioloigcal replicates in each condition

# Setup -------------------------------------------------------------------

library(tidyverse)
library(tximport)
library(DESeq2)
library(ensembldb)
library(EnsDb.Hsapiens.v79)
library(ReportingTools)

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
replicates <- as.numeric(args[6])

# Test values

basedir <- "/home/slee/outputs/bcl6_paper/salmon"
outdir <-  "/home/slee/outputs/test"
cores <- 8
treatment <- "25uM"
control <- "DMSO"
replicates <- 3

# whether the results are for kallisto or sleuth data
type <- basename(basedir)

if (treatment == "25uM"){
  condt <- "FX1"
} else {
  condt <- treatment
}


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

sample_table$condition <- (as.factor(c(rep(condt, replicates),
                                       rep(control, replicates))))

sample_table$condition <- relevel(sample_table$condition,
                                  paste(condt, replicates, sep = "_") )

rownames(sample_table) <- sample_table$sample

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

txdb <- EnsDb.Hsapiens.v79
k <- keys(txdb, keytype = "GENENAME")
df <- select(txdb, keys = k, keytype = "GENENAME", columns = "TXID")
tx2gene <- df[, 2:1]  # tx ID, then gene ID

# make the file paths
if (type == "kallisto"){
    files <- file.path(sample_table$path, "abundance.h5")
    names(files) <- sample_table$sample
} else {
    files <- file.path(sample_table$path, "quant.sf")
    names(files) <- sample_table$sample
}

## Import the quantification data

txi <- tximport(
    files,
    type = type,
    tx2gene = tx2gene,
    txOut = TRUE)
# TXOut true means it stays at transcript level

# [, "condition", drop = FALSE]

ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = sample_table,
                                   design = ~condition)


# DESeq2 -------------------------------------------------------------------

dds <- DESeq(ddsTxi)

res <- results(dds)

# order by padj

resOrdered <- res[order(res$padj),]

# ReportingTools -------------------------------------------------------------------

des_report <- HTMLReport(shortName = "RNAseq_analysis_with_DESeq2",
    title = "RNA-seq analysis of differential expression using DESeq",
    basePath = file.path(outdir, "report") )

# check the countTable is right
# check what annotations should be set to
publish(
    dds,
    des_report,
    pvalueCutoff = 0.05,
    annotation.db = txdb,
    factor = sample_table$condition,
    reportDir = outdir
    )


# sample Plots -------------------------------------------------------------------

pdf(file.path(outdir, "FX1_DMSO_MA_plot.pdf"), width = 7, height = 4)
plotMA(res, ylim = c(-2,2))
dev.off()

pdf(file.path(outdir, "FX1_DMSO_top_padj_counts.pdf"), width = 7, height = 4)
plotCounts(dds, gene=which.min(res$padj), intgroup="condition")
dev.off()

write_csv(as.data.frame(resOrdered), path = file.path(outdir, "FX1_DMSO_deseq2_results.csv"))