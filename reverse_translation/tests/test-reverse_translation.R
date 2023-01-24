library(Biostrings)
library(testthat)
source("./reverse_translation/reverse_translation.R")

test_seqs <- readAAStringSet("./reverse_translation/data/test.fasta")

reverse_translate(test_seqs[1], codon_usage_tables = codon_usage_tbl_yeast_split)

reverse_translated <- list()

for (i in seq_len(length(test_seqs))) {
  cat("Processing ", i, "\n")
  reverse_translated[[i]] <- reverse_translate(aa_seq = test_seqs[i], codon_usage_tables = codon_usage_tbl_yeast_split)
}

(reverse_translated <- DNAStringSet(reverse_translated))


