#!/usr/bin/env nextflow
import java.nio.file.Paths
import java.nio.file.Files

assert Files.exists(Paths.get(params.count_dir))
odir = new File(params.out_dir)
if (!odir.exists()) {
  odir.mkdirs()
}
 
sample_names = file('input_data/yang_sample_names.txt')	// list of samples used in project
bamfile_names = file('input_data/bamfilenames.txt')	// all possible bamfilenames
coding_snps = file('tmp_data/unambigous_coding_snps.vcf.bgz')
coding_snps_idx = file('tmp_data/unambigous_coding_snps.vcf.bgz.tbi')
assert Files.exists(sample_names)
assert Files.exists(bamfile_names)
assert Files.exists(coding_snps)
assert Files.exists(coding_snps_idx)

process generate_experiment_metadata {
  storeDir "$params.out_dir"
  input:
    file sample_names
    file bamfile_names
  output:
    file 'experiment_metadata.csv' into expmeta
  """
    generate_metadata_file.py $sample_names $bamfile_names -o experiment_metadata.csv
  """
}

expmeta
  .splitCsv(header: true)
  .map { 
    [it.sample_name, it.genotype_name, it.timepoint, Paths.get(params.count_dir, it.sample_name + '.csv') ] }
  .set{ ase_counts_ch}


process filter_ASE_counts {
  input:
    file vcf from coding_snps
    file vcf_idx from coding_snps_idx
    set sample_name, genotype_name, timepoint, 'counts.csv' from ase_counts_ch
  output:
    set sample_name, genotype_name, timepoint, "${sample_name}.csv" into filtered_ase_counts_ch
    """
    filter_ase_counts.py --lim-reads 30 --vcf $vcf --counts 'counts.csv' -n $genotype_name -o "${sample_name}.csv"
    """
}

process test_ASE_counts {
  publishDir "$params.out_dir/ASE_test", mode: 'copy'
  input:
    set sample_name, genotype_name, timepoint, filtered_counts from filtered_ase_counts_ch
  output:
    file "${sample_name}.csv" into tested_ase_counts_ch
    """
    ase_binom_test.R -p 1e-2 -l 20 $filtered_counts -t $timepoint -s $sample_name -d $genotype_name "${sample_name}.csv"
    """
}

snp_gene_map = file('tmp_data/snp_gene_map.tsv')

process combine_significant_ASE {
  publishDir "$params.out_dir", mode: 'copy'
  input:
    file csv_files from tested_ase_counts_ch.toSortedList() 
    file snp_gene_map
  output:
    file 'ASE_effects_full.csv' into ASE_ch
    """
    merge_and_annotate_ASE_results.R ASE_effects_full.csv $snp_gene_map $csv_files 
    """
}

process combine_significant_ASE {
  publishDir "$params.out_dir", mode: 'copy'
  input:
    file ASE_ch
  output:
    file 'ASE_effects.csv' into sigASE_ch
    """
    #!/usr/bin/env Rscript
    if (!require("pacman")) install.packages("pacman")
    pacman::p_load(tidyverse)
    df <- read_csv('$ASE_ch')
    df %>% 
      filter(padj < 1e-2) %>% 
      dplyr::select(SNP=variantID, gene_id, gene_name, chrom=contig, pos=position, sample, donor, timepoint, padj) %>% 
      write_csv('ASE_effects.csv')
    """
}
