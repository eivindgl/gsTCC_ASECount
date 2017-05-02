#!/usr/bin/env nextflow
import java.nio.file.Paths

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

