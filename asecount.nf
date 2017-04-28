#!/usr/bin/env nextflow
import java.nio.file.Paths

Channel
  .fromPath('tmp_data/run_meta.csv')
  .splitCsv(header: true)
  .map { 
    [it.sample_name, Paths.get(params.bam_base_dir, it.bam) ] }
  .set{ sample_meta }

coding_snps = file('tmp_data/unambigous_coding_snps.vcf.bgz')

process ASECount {
  storeDir 'out/ase_counts'
  module 'GATK'
  executor 'slurm'
  cpus 2
  memory '6 GB'
  time '4 h'
  input:
    file vcf from coding_snps
    set sample_name, 'rna.bam' from sample_meta
  output:
    set sample_name, "${sample_name}.ASE.csv" into ase_counts 
    """
    java \
     -XX:ParallelGCThreads=1  -Xmx4g \
     -jar \${EBROOTGATK}/GenomeAnalysisTK.jar \
     -R "$params.genome_ref" \
     -T ASEReadCounter \
     -o "${sample_name}.ASE.csv" \
     -I "$bam" \
     -sites "$vcf" \
     -L "$vcf" \
     -U ALLOW_N_CIGAR_READS \
     --minMappingQuality 10
    """
}

  
