# Created 03/06/2017
# Script to run Wasabi on Sleuth samples
library(wasabi)


args <- commandArgs(trailingOnly = TRUE)

res_dir <- args[1]

sfdirs <- list.dirs(res_dir, recursive = FALSE)

prepare_fish_for_sleuth(sfdirs)
