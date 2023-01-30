reverse_translate <- function(aa_seq, codon_usage_tables, aa_seq_name = "", exclude_sites = NULL) {
  # a wrapper around reverse_translate which will exclude some restriction enzyme sites
  # aa_seq: a sequence of amino acids
  # exclude_sites: a character vector of restriction enzyme cutting sites that 
  # need to be excluded
  if (! "codon_frequency" %in% colnames(bind_rows(codon_usage_tables))) {
    stop("Column 'codon_frequency' missing in codon_usage_tables. Please provide a column with the frequencies at which codons should be picked.")
  }
  aa_seq <- toupper(as.character(aa_seq))
  allowed_letters <- c("G", "A", "L", "M", "F", "W", "K", "Q", "E", "S", 
                       "P", "V", "I", "C", "Y", "H", "R", "N", "D", "T")
  stopifnot(length(aa_seq) == 1)
  aa_seq_split <- strsplit(aa_seq, split = "")[[1]]
  aa_seq_split <- aa_seq_split[aa_seq_split %in% allowed_letters]
  stopifnot(all(aa_seq_split %in% allowed_letters))
  aa_seq <- paste0(aa_seq_split, collapse = "")
  
  nt_seq = character(length = length(aa_seq_split))
  cur_ind <- 1
  last_ind <- length(aa_seq_split)
  n_backtracks <- 0
  while (cur_ind <= last_ind) {
    #test_passed <- FALSE
    n_repeats <- 0
    weight_var <- "codon_frequency"
    repeat {
      aa <- aa_seq_split[cur_ind]
      
      if (n_repeats > 20) weight_var <- NULL
      
      codon <- aa_to_codon_optimised(aa, 
                                     codon_usage_table = codon_usage_tables[[aa]], 
                                     weight_var = weight_var)
      nt_seq[cur_ind] <- codon
      n_repeats <- n_repeats + 1
      if (is.null(exclude_sites) || (!contains_exclude_sites(paste0(nt_seq[1:cur_ind], collapse = ""), exclude_sites))) {
        break
      } else if (n_repeats > 50) {
        # if there are too many tries, backtracking occurs
        cur_ind <- pmax(cur_ind - 5, 1)
        n_backtracks <- n_backtracks + 1
        break
      }
    }
    ## TO DO: this should only happen if backtracking did not occur
    cur_ind <- cur_ind + 1
    if (n_backtracks > length(aa_seq_split)) {
      warning("It seems to be impossible to perform reverse translation under the given constraints, returning empty DNA string instead.")
      return(DNAString())
    } 
  }
  nt_seq_collapsed <- DNAString(paste0(nt_seq, collapse = ""))
  # make the test:
  ## no.init.codon = TRUE: important! otherwise if and only if TTG and CTG are the first 
  # codon they will be translated as M instead of L, see here:
  # https://github.com/Bioconductor/Biostrings/issues/15
  backtranslated <- as.character(Biostrings::translate(nt_seq_collapsed, 
                                                       no.init.codon = TRUE))
  stopifnot(identical(unname(aa_seq), unname(backtranslated)))
  return(nt_seq_collapsed)
}


aa_to_codon_optimised <- function(aa, codon_usage_table, 
                                  weight_var = "codon_frequency") {
  # reverse translates a single amino acid using a codon usage table
  if (!is.null(weight_var)) {
    probs = codon_usage_table[[weight_var]]
  } else {
    probs = NULL
  }
  codon <- sample(x = codon_usage_table[["codon"]], 
                  size = 1, 
                  prob = probs)
  return(codon)
}


contains_exclude_sites <- function(nt_seq, exclude_sites) {
  # nt_seq: a DNA sequence
  # exclude_sites: a character vector with sequences to screen for
  nt_seq <- DNAString(nt_seq)
  stopifnot({
    is.character(exclude_sites)
  })
  tests_fwd <- map_dbl(exclude_sites, function(pattern) {
    countPattern(pattern = pattern, subject = nt_seq, fixed = FALSE)
  })
  tests_rev <- map_dbl(exclude_sites, function(pattern) {
    countPattern(pattern = pattern, subject = reverseComplement(nt_seq), fixed = FALSE)
  })
  tests <- c(tests_fwd, tests_rev)
  return(any(tests > 0))
}


get_gc_content <- function(nt_seq) {
  nt_seq_split <- strsplit(nt_seq, split = "")[[1]]
  gc_content <- round(mean(nt_seq_split %in% c("G", "C")), digits = 3)
  return(gc_content)
}

enzymes <- 
  read_tsv("./data/enzymes_and_sequences.tsv", 
           col_names = c("enzyme", "cut_site")) %>% 
  deframe()

codon_usage_tables <- readRDS("./data/codon_usage_tables.rds")


# codon_usage_tables <- list(
#   "sc" = codon_usage_tbl_yeast, 
#   "ec" = NULL, 
#   "hs" = NULL
# )

codon_usage_tables_split <- 
  lapply(codon_usage_tables, function(tab) {
    if (is.null(tab)) return(NULL)
    split(tab, tab$aa_single)
  })


  

