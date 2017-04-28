module load GATK
GENOME_BUILD='/groups/umcg-wijmenga/tmp04/resources/b37/indices/human_g1k_v37.fasta'
vcf_file='input_data/gsTCC.RNA.all.vcf.gz'
java -Xmx2g -jar ${EBROOTGATK}/GenomeAnalysisTK.jar -T SelectVariants \
	-R $GENOME_BUILD \
	-V $vcf_file \
	-selectType SNP \
	-o gsTCC.RNA.SNPs.vcf.gz

