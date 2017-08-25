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
# 2. paths to the aux files
# 3. number of cores
#
# Done in parallel for both transcript Level
# and for collapsed gene level woth tximport
# at this point tx level only for salmon

# Setup -----------------------------------------------

suppressPackageStartupMessages(library(tidyverse, quietly = TRUE))
suppressPackageStartupMessages(library(BiocParallel, quietly = TRUE))
suppressPackageStartupMessages(library(tximport, quietly = TRUE))
suppressPackageStartupMessages(library(DESeq2, quietly = TRUE))


args <- commandArgs(trailingOnly = TRUE)
# test if there is at least one argument:
# if not, return an error
# if (length(args) < 6) {
#      stop("Six arguments must be supplied.", call. = FALSE)
# }

top.dir <- args[1]
aux.dir <- args[2]
cores <- as.numeric(args[3])


dte <- format(Sys.time(), "%Y-%m-%d")

# Set number of cores for multi-processing
register(MulticoreParam(cores))


# make the various directory paths

ref.dir <- file.path(top.dir, "reference_files")
res.dir <- file.path(top.dir, "deseq2_results")


# Functions ---------------------------------------------------------------

tx_gene_conv <- function(tx, tx2g){
    # requires tx2g be loaded
    # takes a vec of TX_IDS and returns a vec of HGNC_IDS
    vec <-  match(tx, tx2g$TXID)
    gene <- tx2g[vec, 2]
    return(gene)
}

deseq_test <- function(dds.x, treat, control, level){
  # Calculates Wald test results for specified treatment
  #   and control in the DeSeq2 object passed in
  # dds.x must be made with betaPrior = TRUE
  # Saves the results table to file
  # Saves the head of the raw results so test details etc
  #   can be tested
  # also saves an MA plot of the data
  # level is gene or transcript - just used for naming
  #   the results dir
  # ---
  # TEST VALS
  # dds.x <- dds.tx
  # treat <- "25uM"
  # control <- "DMSO"
  # level <- "transcript"
  # ---
  ## create results dir
  results.dir <- file.path(res.dir, level, paste(treat, control, sep ="_"))
  dir.create(results.dir, recursive = TRUE, showWarnings = F)
  res.x <- results(
    dds.x,
    contrast = c("condition", treat, control)
  )
  res.x.ordered <- res.x[order(res.x$padj),]
  ## Save the MA plot
  pdf(file.path(
    results.dir,
    paste(dte, treat, control, "ma-plot.pdf", sep="_"))
      )
  plotMA(res.x.ordered, ylim = c(-2, 2))
  dev.off()
  ## save the results file
  write.csv(
    as.data.frame(res.x.ordered),
    file.path(
      results.dir,
      paste(dte, treat, control, level, "deseq2_results.csv", sep="_")
      ),
    quote = FALSE
  )
  ## save info about test
  # sink(
  #   file.path(
  #     results.dir,
  #     paste(dte, treat, control, level, "deseq2_info.csv", sep="_")
  #   )
  # )
  # res.x@elementMetadata
  # summary(res.x.ordered)
  # sink()
  write.csv(
    as.data.frame(res.x.ordered@elementMetadata),
    file.path(
      results.dir,
      paste(dte, treat, control, level, "deseq2_info.csv", sep="_")
    )
  )
}



# read the gene - tx mappings  --------------------------------------------

tx2gene <- read_csv(file.path(ref.dir, "tx2g_table.csv"))

# sample table construction --------------------------

condition.table <- read_csv(file.path(aux.dir, "condition_table.txt"))

# build the list of sample names from the condition table
# they will be foo_1 .. foo_n etc
samples <- c()
for (i in 1:length(condition.table$condition)){
  reps <- condition.table$replicates[[i]]
  x <-  paste(rep(condition.table$condition[[i]], reps),
              rep(1:reps),
              sep = "_"
              )
  samples <- append(samples, x)
}

# create sample table with sample names, conditions, and paths to quant.sf's
sample_table <- cbind(samples,
                      gsub(pattern = "_[0-9]+", "", samples),
                      file.path(top.dir, "kallisto", samples, "abundance.h5"))
sample_table <- as.data.frame(sample_table, stringsAsFactors = FALSE)
colnames(sample_table) <- c("sample", "condition", "path")
# make conditions factors for DESeq2 contrast testing.
sample_table$condition <- factor(
  sample_table$condition,
  levels = condition.table$condition
)


## print to display as a check while program runs
print("Sample table:", quote = FALSE)
print(sample_table)

## check that all dirs in the table exist and exit if not
dircheck <- all(file.exists(sample_table$path))
if (! dircheck == TRUE){
  stop("All the quantification directories for the conditions do not exist.")
}

# tximport ---------------------------------------------

## Import the quantification data

print("Importing the kallisto quantifications...", quote = FALSE)

txi.tx <- tximport(
  sample_table$path,
  type = "kallisto",
  tx2gene = tx2gene,
  txOut = TRUE,
  countsFromAbundance = "lengthScaledTPM",
  ignoreTxVersion = T
  )
txi.gene <- tximport(
  sample_table$path,
  type = "kallisto",
  tx2gene = tx2gene,
  txOut = FALSE,
  countsFromAbundance = "lengthScaledTPM",
  ignoreTxVersion = T
  )

rownames(txi.tx$counts) <- as.character(rownames(txi.tx$counts))
rownames(txi.gene$counts) <- as.character(rownames(txi.gene$counts))

ddsTxi.tx <- DESeqDataSetFromTximport(txi.tx,
                                   colData = sample_table,
                                   design = ~condition)

ddsTxi.gene <- DESeqDataSetFromTximport(txi.gene,
                                    colData = sample_table,
                                    design = ~condition)

# DESeq2 --------------------------------------------

print("Creating transcript level DESeq2 object", quote = FALSE)
dds.tx <- DESeq(ddsTxi.tx, parallel = TRUE, betaPrior = TRUE)

print("Creating gene level DESeq2 object", quote = FALSE)
dds.gene <- DESeq(ddsTxi.gene, parallel = TRUE)

# shrink models with rld for visualisation
print("Regular log transform the transcript level data", quote = FALSE)
rld.tx <- rlog(dds.tx, blind = FALSE)
print("Regular log transform the gene level data", quote = FALSE)
rld.gene <- rlog(dds.gene, blind = FALSE)

# Make PCA plots of the sample distances

dds.pca.tx <-  plotPCA(rld.tx)
dds.pca.gene <-  plotPCA(rld.gene)




# Do the test specified in de_test_table ----------------------------------

de.tests.table <- read_csv(file.path(aux.dir, "de_tests_table.txt"))
colnames(de.tests.table) <- toupper(colnames(de.tests.table))

# Transcript level
print("Extracting transcript level comparisons", quote = FALSE)
mapply(deseq_test,
       treat = de.tests.table$TREATMENT,
       control = de.tests.table$CONTROL,
       MoreArgs = list(dds.x = dds.tx,
                       level = "transcript"))

# Gene level
print("Extracting gene level comparisons", quote = FALSE)
mapply(deseq_test,
       treat = de.tests.table$TREATMENT,
       control = de.tests.table$CONTROL,
       MoreArgs = list(dds.x = dds.gene,
                       level = "gene"))

# save the PCAs -----------------------------------------------------------

ggsave(
  file.path(res.dir, "transcript", paste(dte, "transcrpt_PCA.pdf", sep = "_")),
  dds.pca.tx
)

ggsave(
  file.path(res.dir, "gene", paste(dte, "gene_PCA.pdf", sep = "_")),
  dds.pca.gene
)

# save session info --------------------------------------

sink(file =file.path(res.dir, paste(dte, "dseq2_sessionInfo.txt", sep = "_")))
sessionInfo()
sink()


print("Remember to check all info files to ensure accurate results", quote = FALSE)
