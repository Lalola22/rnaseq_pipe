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

# 1. top level results directory
# 2. quants to use (kallisto or salmon)
# 3. number of cores
# 4. treatment
# 5. control
# 6. replicates
# Done in parallel for both transcript Level
# and for collapsed gene level
# at this point tx level only for salmon

# Setup -----------------------------------------------

library(tidyverse, quietly = TRUE)
library(BiocParallel, quietly = TRUE)
library(tximport, quietly = TRUE)
library(DESeq2, quietly = TRUE)
library(ReportingTools, quietly = TRUE)
library(org.Hs.eg.db, quietly = TRUE)

args <- commandArgs(trailingOnly = TRUE)
# test if there is at least one argument:
# if not, return an error
if (length(args) < 6) {
     stop("Six arguments must be supplied.", call. = FALSE)
}

top_dir <- args[1]
type <- args[2]
cores <- as.numeric(args[3])
treatment <- args[4]
control <- args[5]
replicates <- as.numeric(args[6])

# Test values

# top_dir <- "/home/slee/outputs/test/june_09"
# type <- "salmon"
# cores <- 8
# treatment <- "25uM"
# control <- "DMSO"
# replicates <- 3

# make the various directory paths

quant_dir <- file.path(top_dir, type)
ref_dir <- file.path(top_dir, "reference_files")
res_dir <- file.path(
    top_dir,
    "deseq2",
    paste(treatment, control, sep = "_"))

dir.create(res_dir, showWarnings = FALSE)

#### REMOVE BEFORE USAGE -- TESTING ONLY

quant_dir <- "/home/slee/outputs/bcl6_paper/salmon"

if (type == "kallisto"){
    stop("Kallisto isn't working for scripting yet!!")
}

# if 25uM is the treatment
# we want to call it FX1
if (treatment == "25uM"){
  condt <- "FX1"
} else {
  condt <- treatment
}

# Set number of cores for multi-processing

register(MulticoreParam(cores))


# sample table construction --------------------------

sample_table <- as.data.frame(
    matrix("", nrow = replicates * 2, ncol = 3),
    stringsAsFactors = FALSE)

colnames(sample_table) <- c("sample", "condition", "path")

sample_table <- sample_table %>%
    mutate(
        sample = paste(c(rep(treatment, replicates),
                            rep(control, replicates)),
                        c(1:replicates),
                        sep = "_")
        ) %>%
    mutate(path = file.path(quant_dir, sample)) %>%
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

# tximport ---------------------------------------------

## read in the tx2g table

tx2gene <- read_csv(file.path(ref_dir, "tx2g_table.csv"))

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

# DESeq2 --------------------------------------------

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

des_report <- HTMLReport(
    shortName = "RNAseq_analysis_with_DESeq2",
    title = paste("Gene Level RNA-seq analysis", condt, control sep = "_"),
    basePath = file.path(res_dir, "report") )
publish(
    dds.gene,
    des_report,
    pvalueCutoff = 0.05,
    annotation.db = "org.Hs.eg.db",
    factor = sample_table$condition,
    reportDir = outdir
    )
finish(des_report)

# sample Plots -------------------------------------------

pdf(file.path(
    res_dir, "FX1_DMSO_MA_plot.pdf"),
    width = 7, height = 4)
plotMA(res, ylim = c(-2,2))
dev.off()

write.csv(
    as.data.frame(res_ordered),
    file = file.path(
        res_dir,
        "tx_deseq2_results.csv")
        )

write.csv(
    as.data.frame(res_ordered.gene),
    file = file.path(
        res_dir,
        "gene_deseq2_results.csv")
        )


# save session info --------------------------------------


sink(file =file.path(res_dir, "dseq2_sessionInfo.txt") )
sessionInfo()
sink()
