# plot_figures.R
# Created 18/10/2017
# Sam Lee
# creation of plots based on (a) DE test and (b) level [gene, transcript]
#   Volcano plot
#
# To create custom figures read the docs for the functions for additional options
# to change colours / gene sets / cutoffs etc
#
# Takes XXX arguments
# 1. top.dir (top output dir for analysis)
# 2. comp (the de test to make e.g. condtA_condtB)
# 3. (level [gene,transcript])
# ---
# packages required
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(RColorBrewer))

# Arguments 
args <- commandArgs(trailingOnly = TRUE)

top.dir <- args[1]
comp <- args[2]
level <- as.numeric(args[3])

# parameters
l2fc <- 1.5
pva <- 0.05

# Functions
plot.volcano <- function(deseq.df, l2fc, pval, id_list = NULL, sig_col = "#ca0020"){
  # function to create a volcano plot of DE results data at specified
  #     cutoffs for significance.
  # Optionally can highlight certain features (genes / tx's)
  # Can change the color for sig features using sig_col (default red) 
  #     This should ideally be a hex code but is not strict
  # returns a ggplot object
  # ---
  # deseq.df is a deseq2 results file read in with the following columns kept
  #   ID, log2FoldChange, padj
  # l2fc / padj are obvious
  # id_list is a vector containing features id's to be highlighted.
  # ---
  # filter table for sig results + add label status row
  sig.tbl <- deseq.df %>%
    mutate(significant = abs(log2FoldChange) > l2fc & padj < pval) %>%
    filter(significant == TRUE | significant == FALSE) %>%
    mutate(lab = ID %in% id_list)
  
  # create the color mapping # sig_col (default red) - grey
  my.colours <- c(sig_col, "#bababa")
  names(my.colours) <- c(TRUE, FALSE)
  colScale <- scale_colour_manual(name = "significant", values = my.colours)
  
  # make the volcano plot object
  plot.volc <- ggplot() +
    geom_point(data = sig.tbl, aes(x = log2FoldChange, y = -log10(padj), col = significant), alpha = 0.5) +
    # geom_hline(aes(yintercept= -log10(pval))) + 
    colScale +
    geom_point(data = filter(sig.tbl, lab == TRUE), aes(x = log2FoldChange, y = -log10(padj)), col = "black", size = 0.5) +
    theme_pubr()
  
  # return the plot 
  return(plot.volc)
}

# read in the deseq2 results file
dir.path <- file.path(top.dir, "deseq2_results", level, comp, )
file.path <- list.files(dir.path, pattern = "results", full.names = T)

deseq.file <- read_csv(file.path)


# make plot object
plot.volc <- plot.volcano(dplyr::select(deseq.file, ID, log2FoldChange, padj), l2fc, pval, sig_col = dmso.col)


# save plot objects
ggsave(file.path(dir.path, paste(comp, "volcano-plot.pdf", sep = "-")), plot = plot.volc)