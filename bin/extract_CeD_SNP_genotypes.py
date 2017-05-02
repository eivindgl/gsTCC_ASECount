import collections
import vcf
import pandas as pd

with open('input_data/CeD_tag_SNPs.bed') as f:
    tag_snp_IDs = {x.strip().split()[-1] for x in f}

vx = vcf.Reader(filename='input_data/gsTCC.IChip.all.vcf.gz', compressed=True)

sample_names = sorted(vx.samples)
tag_snps = [x for x in vx if x.ID in tag_snp_IDs]
tag_snps_not_ichip = tag_snp_IDs.difference(x.ID for x in tag_snps)
print('The following CeD tag SNPs are not covered on immunochip:', ', '.join(tag_snps_not_ichip))

d = collections.defaultdict(list)
for snp_row in tag_snps:
    for sample in sample_names:
        variant = snp_row.genotype(sample).gt_bases
        d[sample].append(variant)

df = pd.DataFrame(d, index=[x.ID for x in tag_snps])
df.to_csv('tmp_data/tag_snp_genotype.csv', index_label='SNP')
