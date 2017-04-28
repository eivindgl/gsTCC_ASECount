# Writes coding snps to tmp_data/unambigous_coding_snps.vcf.bgz
# Writes snp/gene map to tmp_data/snp_gene_map.tsv
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  tidyverse,
  stringr,
  rtracklayer,
  assertthat,
  GenomicFeatures,
  VariantAnnotation
)

gene_annotation <- import('../tmp_data/gencode.v26lift37.annotation.gff3.gz')
tag_snps <- import('../input_data/CeD_tag_SNPs.bed')
hits <- tag_snps %>% 
  flank(5e5, both = TRUE) %>% 
  findOverlaps(gene_annotation)
# a verbose way of finding overlapping genes
overlapping_genes <- gene_annotation[to(hits)] %>%
  mcols %>% 
  as.data.frame() %>% 
  filter(type == 'gene') %>% 
  distinct(ID) %>% 
  `[[`('ID')
cat('We got', length(overlapping_genes), 'genes within proximity of CeD tag SNPs.\n')


# I don't know how to read a subset of the gencode annotation
# I therefore compute all the exons and then select a subset.
gdb <- makeTxDbFromGFF('../tmp_data/gencode.v26lift37.annotation.gff3.gz', organism = 'Homo sapiens')
gexons <- exonsBy(gdb, by = 'gene')

# all overlapping genes must also be present in gexons list
assertthat::are_equal(
  intersect(overlapping_genes, names(gexons)),
  overlapping_genes)
uxons <- reduce(unlist(gexons[overlapping_genes]), ignore.strand = TRUE) # remove overlap

rna_snps <- readVcf('../input_data/gsTCC.RNA.SNPs.vcf.gz')
snps <- rowRanges(rna_snps)
seqlevelsStyle(snps) <- 'UCSC'
snpExonHits <- findOverlaps(snps, uxons)
  
coding_snps <- rna_snps[from(snpExonHits), ]

coding_snps %>% 
  writeVcf('../tmp_data/coding_snps.gz', index = TRUE)

snps <- rowRanges(coding_snps)
seqlevelsStyle(snps) <- 'UCSC'
snps <- keepSeqlevels(snps, paste0('chr', 1:22))

# filter low mappability snps
# check simulation filter from review
# some mapscores span multiple bases and may therefore cover multiple SNPs
mapscore <- import.bw('../tmp_data/wgEncodeCrgMapabilityAlign50mer.bigWig', which = snps)
num_SNPs_overlapping_multiple_mapscores <- countOverlaps(snps, mapscore) %>% 
  purrr::keep(~ . > 1) %>% 
  length()
assert_that(are_equal(num_SNPs_overlapping_multiple_mapscores, 0))
idx <- findOverlaps(snps, mapscore)
assert_that(are_equal(length(snps), max(from(idx)))) # from(idx) index snps
assert_that(are_equal(length(mapscore), max(to(idx)))) # to(idx) index mapscore
assert_that(are_equal(from(idx), sort(from(idx)))) # snps are the reference
snp_score <- unlist(mapscore[to(idx)])$score
filtered_coding_snps <- coding_snps[which(snp_score >= 1), ]
filtered_coding_snps %>% 
  writeVcf('../tmp_data/unambigous_coding_snps.vcf', index = TRUE)
snps <- snps[which(snp_score >= 1), ]

snp_genes <- findOverlaps(gene_annotation, snps)
x <- unlist(gene_annotation[from(snp_genes), ])
df <- tibble(snp=names(snps[to(snp_genes), ]),
             gene=x$ID,
             type=x$type) %>% 
  filter(type == 'gene') %>% 
  dplyr::select(-type)
df %>%
  write_tsv('../tmp_data/snp_gene_map.tsv')

## read SNPs filtered by Panousis et al 2014
## all SNPs in this list are associated with ambigous mapping
## COMMENTED OUT BECAUSE THIS IS A SUBSET OF ENCODE mappability < 1 data
# df <- read_tsv('../tmp_data/EUR01_50bp_result_stats_05bias.txt',
#               col_types = cols(CHR = col_character())) %>% 
#   filter(str_detect(VARIANT_ID, 'snp')) # only include SNPs (not indels etc) 
# sim_snps <- GRanges(seqnames = df$CHR, 
#         ranges = IRanges(start = df$POSITION, width = 1))
# seqlevelsStyle(sim_snps) <- 'UCSC'
# sim_idx <- findOverlaps(snps, sim_snps)
  
  
# filter(abs(XTU_CORRECT_RATIO_SCALED - 0.5) > 0.05) %>% 