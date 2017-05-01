#!/usr/bin/env Rscript
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  tidyverse,
  stringr,
  assertthat
)

args = commandArgs(trailingOnly=TRUE)
assert_that(length(args) > 2)
out_path <- args[1]
snp_gene_path <- args[2]
input_csv_files <- args[3:length(args)]
assert_that(is.count(length(input_csv_files)))
assert_that(all(map_lgl(input_csv_files, file.exists)))

assert_that(file.exists(snp_gene_path))
snp_gene_map <- read_tsv(snp_gene_path) %>% 
  dplyr::rename(variantID = snp)

input_csv_files %>% 
  map(read_csv) %>% 
  bind_rows() %>% 
  inner_join(snp_gene_map) %>% 
  write_csv(out_path)
