library(tidyverse)
library(Biostrings)

setwd("./reverse_translation/")
source("./reverse_translation.R")
proteins <- readAAStringSet("./data/test.fasta")

head(enzymes)

str_split(paste0(enzymes, collapse = ""), pattern = "")[[1]] %>% unique()

protein_nt_seqs <- 
  lapply(proteins, reverse_translate, 
         codon_usage_tables = codon_usage_tbl_split, 
         exclude_sites = enzymes[c("BsaI", "NotI", "PacI", "SbfI", "SgfI", "AscI")]) %>% 
  DNAStringSet()
protein_nt_seqs



