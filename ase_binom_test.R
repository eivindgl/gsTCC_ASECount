#!/usr/bin/env Rscript
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  tidyverse,
  stringr,
  purrr,
  docopt,
  assertthat
)


parse_numeric_arg <- function(s, default_value) {
  if (is_null(s)) {
    default_value
  } else {
    assert_that(!is_null(as.numeric(s)))
    as.numeric(s)  
  }
}

'usage: ase_binom_test.R [-l <min_reads> -p <padj_cutoff> -s <sample_name> -d <donor_name>] <input_path> <output_path>
options:
  -l <min_readsl>     Drop SNPs with fewer reads
  -p <padj_cutoff>    Drop SNPs above FDR cutoff
  -s <sample_name>    Add sample name to output files
  -d <donor>          Add donor name
' -> doc

# load the docopt library
library(docopt)
#a <- docopt(doc, 'out/filtered_ase_counts/batch3__TCC.06.1__t10.csv out/sigfiltered/batch3__TCC.06.1__t10.csv')
a <- docopt(doc)
input_path <- a$input_path
output_path <- a$output_path
min_total_count <- parse_numeric_arg(a$l, default = 20)
padj_lim <- parse_numeric_arg(a$p, default = 1e-3)

dir.create(dirname(output_path), showWarnings = TRUE)
df <- read_csv(input_path)
if (!is_null(a$s)) {
  df <- df %>% 
    mutate(sample = a$s)
}

if (!is_null(a$d)) {
  df <- df %>% 
    mutate(donor = a$d)
}
df %>% 
  filter(totalCount >= min_total_count) %>% 
  mutate(pval = map2_dbl(refCount, totalCount, ~ binom.test(.x, .y)$p.value),
         padj = p.adjust(pval)) %>% 
  dplyr::select(padj, pval, everything()) %>% 
  arrange(pval) %>% 
  filter(padj < 1e-2) %>% 
  write_csv(output_path)
