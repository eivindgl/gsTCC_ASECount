#!/usr/bin/env nextflow
import java.nio.file.Paths

# TODO add all calculon specific code here

rna_snps_all = file('tmp_input_data/gsTCC.RNA.all.vcf.gz')

process filter_SNPs_from_RNA_genotypes {
  module  'GATK'
  input:
    file rna_snps_all
  output:
    file 'gsTCC.RNA.SNPs.vcf.gz' into rna_snps_ch
    """
    java -Xmx2g -jar ${EBROOTGATK}/GenomeAnalysisTK.jar -T SelectVariants \
	-R "$params.genome_ref" \
	-V "$rna_snps_all" \
	-selectType SNP \
	-o gsTCC.RNA.SNPs.vcf.gz

    """
}


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



  
