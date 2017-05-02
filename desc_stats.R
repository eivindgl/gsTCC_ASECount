#!/usr/bin/env Rscript
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  tidyverse,
  stringr,
  assertthat
)

df <- read_csv('out/sig_ASE_effects.csv') %>% 
  mutate(timepoint = str_extract(sample, 't\\d+$'),
         timepoint = factor(timepoint, levels = c('t0', 't10', 't30', 't180')))

df %>% 
  group_by(variantID, gene_id, donor) %>% 
  summarize(donor_effect = n(), gene_name = unique(gene_name)) %>% 
  arrange(desc(donor_effect)) %>% 
  print(n = 30)

# we have 87 samples with significant results (of 88 possible)
df %>% 
  distinct(sample) %>% 
  nrow()
  

df %>% 
  distinct(donor, sample) %>% 
  group_by(donor) %>% 
  count() %>% 
  arrange(desc(n))

df %>% 
  filter(gene_name == 'CCR5') %>% 
  distinct(donor, sample) %>% 
  group_by(donor) %>% 
  count() %>% 
  arrange(desc(n))

df %>% 
  filter(gene_name == 'CCR5' & donor == 'TCC-01') %>% 
  dplyr::select(donor, sample, gene_name, variantID)


##
## Identify timepoint specifc ASE patterns
##
gdf <- df %>% 
  distinct(gene_id, donor, sample, timepoint) %>% 
  group_by(gene_id, donor, timepoint) %>% 
  summarize(ase = n() > 0) %>% 
  ungroup() %>% 
  group_by(gene_id, timepoint) %>% 
  summarize(n = n()) %>% 
  spread(timepoint, n, fill = 0)

gdf <- gdf %>% 
  gather(timepoint, n, -gene_id) %>% 
  group_by(gene_id) %>% 
  summarise(total = sum(n)) %>% 
  inner_join(gdf) %>% 
  arrange(desc(total))
gdf %>% 
  filter(total > 15)

##
## Join with GWAS risc SNP genotyp info
##
gtdf <- read_csv('tmp_data/tag_snp_genotype.csv')%>%
  gather(donor, genotype, -SNP) %>% 
  dplyr::rename(tag_snp = SNP) %>% 
  mutate(category = if_else(genotype %in% c('A|A', 'C|C', 'G|G', 'T|T'),
                            'homozygous', 'heterozygous'))
gene_tag_snp_map <- read_csv('tmp_data/gene_tag_snp_map.csv')

top_N_genes <- gdf$gene_id %>% head(n = 40)
df %>% 
  filter(gene_id %in% top_N_genes) %>% 
  dplyr::select(SNP=variantID, donor, sample, gene_id, gene_name, timepoint) %>% 
  inner_join(gene_tag_snp_map) %>% 
  inner_join(gtdf) %>% 
  distinct(gene_id,gene_name, donor, category) %>% 
  group_by(gene_id, gene_name) %>% 
  summarize(het = sum(category == 'heterozygous'),
            hom = sum(category == 'homozygous'),
            het_frac = het / (het + hom),
            hom_frac = hom / (het + hom)) %>% 
  arrange(desc(het + hom)) %>% 
  ggplot() +
  geom_jitter(aes(het, hom))

