# Created 08/06/17
# Sam Lee
# For Deseq2 analysis of kallisto results
# takes six inputs
# location of quant files
# location of results to be saved
# number of cores to use.
# treatment condition
# control condition
# number of bioloigcal replicates in each condition

# Done in parallel for both transcript Level
# and for collapsed gene level
# at this point tx level only for salmon

# Setup -------------------------------------------------------------------

library(tidyverse, quietly = TRUE)
library(BiocParallel, quietly = TRUE)
library(tximport, quietly = TRUE)
library(DESeq2, quietly = TRUE)
library(ensembldb, quietly = TRUE)
library(EnsDb.Hsapiens.v79, quietly = TRUE)
library(ReportingTools, quietly = TRUE)

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
outdir <-  "/home/slee/outputs/test/june_09"
cores <- 6
treatment <- "25uM"
control <- "DMSO"
replicates <- 3

# whether the results are for kallisto or sleuth data
type <- basename(basedir)

if (type == "kallisto"){
    stop("Kallisto isn't working for scripting yet!!")
}s

# if 25uM is the treatment
# we want to call it FX1
if (treatment == "25uM"){
  condt <- "FX1"
} else {
  condt <- treatment
}

# Set number of cores for multi-processing

register(MulticoreParam(cores))


# sample table construction -----------------------------------------------

sample_table <- as.data.frame(matrix("", nrow = replicates * 2, ncol = 3),
                                     stringsAsFactors = FALSE)
colnames(sample_table) <- c("sample", "condition","path")


sample_table <- sample_table %>%
    mutate(
        sample = paste(c(rep(treatment, replicates),
                            rep(control, replicates)),
                        c(1:replicates),
                        sep = "_")
        ) %>%
    mutate(path = file.path(basedir, sample)) %>%
    mutate(condition = as.factor(
        c(rep(condt, replicates),
          rep(control, replicates)))) %>%
    mutate(condition = relevel(condition, control))

## check if this is needed...
# rownames(sample_table) <- sample_table$sample

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

# make the quant file paths
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
txi.gene <- tximport(
    files,
    type = type,
    tx2gene = tx2gene)

ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = sample_table,
                                   design = ~condition)

ddsTxi.gene <- DESeqDataSetFromTximport(txi.gene,
                                  colData = sample_table,
                                                  design = ~condition)

# DESeq2 -------------------------------------------------------------------

dds <- DESeq(ddsTxi, parallel = TRUE)
dds.gene <- DESeq(ddsTxi.gene, parallel = TRUE)

res <- results(dds)
res.gene <- results(dds.gene)

# order by padj

res_ordered <- res[order(res$padj),]
res_ordered.gene <- res.gene[order(res.gene$padj),]

print("Transcript Level Results", quote = FALSE)
res_ordered

print("Gene Level Results", quote = FALSE)
res_ordered.gene

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
plotCounts(dds, gene = which.min(res$padj), intgroup = "condition")
dev.off()

write.csv(
    as.data.frame(res_ordered),
    file = file.path(
        outdir,
        paste(condt, control, "gene_deseq2_results_02.csv", sep = "_")) )

# save session info --------------------------------------------------------


sink(file.path()"sessionInfo.txt")
sessionInfo()
sink()
