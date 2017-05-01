#!/usr/bin/env nextflow
import java.nio.file.Paths

Channel
  .fromPath('tmp_data/run_meta.csv')
  .splitCsv(header: true)
  .map { 
    [it.sample_name, it.genotype_name, Paths.get(params.bam_base_dir, it.bam) ] }
  .set{ sample_meta }

coding_snps = file('tmp_data/unambigous_coding_snps.vcf.bgz')
coding_snps_idx = file('tmp_data/unambigous_coding_snps.vcf.bgz.tbi')

process ASECount {
  storeDir 'out/ase_counts'
  module 'GATK'
  //  executor 'slurm'
  cpus 2
  memory '6 GB'
  time '4 h'
  input:
    file vcf from coding_snps
    file vcf_idx from coding_snps_idx
    set sample_name, genotype_name, 'rna.bam' from sample_meta
  output:
    set sample_name, genotype_name, "${sample_name}.ASE.csv" into ase_counts_ch
    """
    java \
     -XX:ParallelGCThreads=1  -Xmx4g \
     -jar \${EBROOTGATK}/GenomeAnalysisTK.jar \
     -R "$params.genome_ref" \
     -T ASEReadCounter \
     -o "${sample_name}.ASE.csv" \
     -I rna.bam \
     -sites "$vcf" \
     -L "$vcf" \
     -U ALLOW_N_CIGAR_READS \
     --minMappingQuality 10
    """
}

process filter_ASE_counts {
  publishDir "out/filtered_ase_counts", mode: 'copy'
  input:
    file vcf from coding_snps
    file vcf_idx from coding_snps_idx
    set sample_name, genotype_name, counts from ase_counts_ch
  output:
    set sample_name, genotype_name, "${sample_name}.csv" into filtered_ase_counts_ch
    """
    filter_ase_counts.py --vcf $vcf --counts $counts -n $genotype_name -o "${sample_name}.csv"
    """
}

process test_ASE_counts {
  publishDir "out/sig_ase", mode: 'copy'
  input:
    set sample_name, genotype_name, filtered_counts from filtered_ase_counts_ch
  output:
    set sample_name, genotype_name, "${sample_name}.csv" into sig_ase_counts_ch
    """
    ase_binom_test.R -p 1e-2 -l 20 $filtered_counts -s $sample_name -d $genotype_name "${sample_name}.csv"
    """
}

snp_gene_map = file('tmp_data/snp_gene_map.tsv')
process combine_significant_ASE {
  publishDir "out", mode: 'copy'
  input:
    file csv_files from sig_ase_counts_ch.map{ it[2] }.toSortedList() 
    file snp_gene_map
  output:
    file 'sig_ASE_effects.csv' into sig_ASE_ch
    """
    merge_and_annotate_ASE_results.R sig_ASE_effects.csv $snp_gene_map $csv_files 
    """
}

sig_ASE = sig_ASE_ch.first()



  
