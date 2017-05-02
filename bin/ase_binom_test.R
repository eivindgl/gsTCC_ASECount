#!/usr/bin/env Rscript
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  tidyverse,
  stringr,
  purrr,
  docopt,
  assertthat
)

'usage: ase_binom_test.R [-l <min_reads> -p <padj_cutoff> -s <sample_name> -d <donor_name> -t <timepoint>] <input_path> <output_path>
options:
  -s <sample_name>    Add sample name to output files
  -d <donor>          Add donor name
  -t <timepoint>      Minutes after stimultation
' -> doc

#a <- docopt(doc, '-t t10 out/filtered_ase_counts/batch3__TCC.06.1__t10.csv out/sigfiltered/batch3__TCC.06.1__t10.csv')
a <- docopt(doc)
input_path <- a$input_path
output_path <- a$output_path
assert_that(!is_null(a$s))
assert_that(!is_null(a$d))
assert_that(!is_null(a$t))
dir.create(dirname(output_path), showWarnings = FALSE)
read_csv(input_path) %>% 
  filter(totalCount > 0) %>% 
  mutate(sample = a[['-s']],
         donor = a[['-d']],
         timepoint = a[['-t']]) %>%
  mutate(pval = map2_dbl(refCount, totalCount, ~ binom.test(.x, .y)$p.value),
         padj = p.adjust(pval)) %>%
  arrange(pval) %>%
   write_csv(output_path)
