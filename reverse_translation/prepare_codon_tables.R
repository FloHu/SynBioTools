library(tidyverse)
library(Biostrings)
library(stringr)

codon_usage_tbl_yeast <- 
  read_tsv("./data/codon_usage_table_highly_expressed_genes.tsv", 
           locale = locale(decimal_mark = ",")) %>% 
  janitor::clean_names() %>% 
  dplyr::rename(codon_frequency = frequency)

codon_usage_tbl_yeast_rare_excluded <- 
  mutate(codon_usage_tbl_yeast, 
         codon_frequency = ifelse(codon_frequency < 0.1, 0, codon_frequency))

single_to_triple_lut <- 
  codon_usage_tbl_yeast[, c("aa_single", "aa_triple")] %>% 
  distinct() %>% 
  deframe()

process_codon_table <- function(path, exclude_rare = FALSE) {
  # expects a table as obtained, for example, from here:
  # https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=413997
  # that has been reformatted to contain the following fields:
  # triplet,frequency_per_thousand,number
  # exclude_rare: if TRUE, will set codon_frequencies < 0.1 to 0
  codon_table <- 
    suppressMessages(read_csv(path)) %>% 
    dplyr::rename("codon" = "triplet")
  
  codon_table$codon <- str_replace_all(codon_table$codon, pattern = "U", replacement = "T")
  
  codon_table$aa_single <- 
    as.character(Biostrings::translate(DNAStringSet(codon_table$codon), no.init.codon = TRUE))
  
  codon_table$aa_triple <- single_to_triple_lut[codon_table$aa_single]
  # sum(codon_usage_tbl_ecoli$frequency_per_thousand)
  
  codon_table <- 
    group_by(codon_table, aa_single) %>% 
    mutate(codon_frequency = round(frequency_per_thousand / sum(frequency_per_thousand), digits = 4)) %>% 
    dplyr::select(aa_single, aa_triple, codon, codon_frequency) %>% 
    arrange(aa_single)
  
  if (exclude_rare) {
    codon_table$codon_frequency <- ifelse(codon_table$codon_frequency < 0.1, 
                                          0, 
                                          codon_table$codon_frequency)
  }
  
  return(codon_table)
}

codon_usage_tbl_ecoli <- process_codon_table(path = "./data/codon_usage_Ecoli.csv") 
codon_usage_tbl_ecoli_rare_excluded <- process_codon_table(path = "./data/codon_usage_Ecoli.csv", exclude_rare = TRUE)

codon_usage_tbl_hsapiens <- process_codon_table(path = "./data/codon_usage_Hsapiens.csv")
codon_usage_tbl_hsapiens_rare_excluded <- process_codon_table(path = "./data/codon_usage_Hsapiens.csv", exclude_rare = TRUE)

codon_usage_tbl_hsapiens

codon_usage_tables <- list(
  "codon_usage_tbl_yeast" = codon_usage_tbl_yeast, 
  "codon_usage_tbl_yeast_rare_excluded" = codon_usage_tbl_yeast_rare_excluded, 
  "codon_usage_tbl_ecoli" = codon_usage_tbl_ecoli, 
  "codon_usage_tbl_ecoli_rare_excluded" = codon_usage_tbl_ecoli_rare_excluded, 
  "codon_usage_tbl_hsapiens" = codon_usage_tbl_hsapiens, 
  "codon_usage_tbl_hsapiens_rare_excluded" = codon_usage_tbl_hsapiens_rare_excluded
)

saveRDS(object = codon_usage_tables, file = "./data/codon_usage_tables.rds")

# split them up
codon_usage_tables_split <- 
  map(codon_usage_tables, function(x) {
    split(x, x$aa_single)
  }) -> foo

saveRDS(object = codon_usage_tables_split, file = "./data/codon_usage_tables_split.rds")


