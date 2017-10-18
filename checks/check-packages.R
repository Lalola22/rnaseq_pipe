# This is just a very basic check for packages
# If one is not available will fail with an error

print("Checking required packages...")

suppressPackageStartupMessages(library(BiocParallel))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(EnsDb.Hsapiens.v86))
suppressPackageStartupMessages(library(EnsDb.Mmusculus.v79))
suppressPackageStartupMessages(library(EnsDb.Rnorvegicus.v79))
suppressPackageStartupMessages(library(ensembldb))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(tximport))

print("All required packages are installed :~)", quote = F)
