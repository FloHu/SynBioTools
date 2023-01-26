library(Biostrings)
library(testthat)
source("./reverse_translation/reverse_translation.R")

enzymes <- 
  read_tsv("/Users/fhuber/PROJEKTE/SynBioTools/reverse_translation/data/enzymes_and_sequences.tsv", 
           col_names = c("enzyme", "cut_site")) %>% 
  deframe()

get_data <- function() {
  codon_usage_tbl_yeast <<- 
    read_tsv("/Users/fhuber/PROJEKTE/SynBioTools/reverse_translation/data/codon_usage_table_highly_expressed_genes.tsv", 
             locale = locale(decimal_mark = ",")) %>% 
    janitor::clean_names()
  # exclude everything below 0.1
  codon_usage_tbl_yeast$frequency_adj <<- ifelse(codon_usage_tbl_yeast$frequency < 0.1, 0, codon_usage_tbl_yeast$frequency)
  
  codon_usage_tbl_yeast_split <- split(codon_usage_tbl_yeast, codon_usage_tbl_yeast$aa_single)
  codon_usage_tbl_yeast_split <<- 
    map(codon_usage_tbl_yeast_split, function(tab) {
      tab$sampling_weight <- tab$frequency / sum(tab$frequency)
      tab$sampling_weight_adj <- tab$frequency_adj / sum(tab$frequency_adj)
      return(tab)
    })
}


test_seqs <- readAAStringSet("./reverse_translation/data/test.fasta")

reverse_translate(test_seqs[1], codon_usage_tables = codon_usage_tbl_yeast_split)

reverse_translated <- list()

for (i in seq_len(length(test_seqs))) {
  cat("Processing ", i, "\n")
  reverse_translated[[i]] <- reverse_translate(aa_seq = test_seqs[i], codon_usage_tables = codon_usage_tbl_yeast_split)
}

(reverse_translated <- DNAStringSet(reverse_translated))


