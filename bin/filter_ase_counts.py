#!/usr/bin/env python3
'''
This scripts filters all heterozygous coding SNPs for a given donor.
Input counts.csv is a list of counts for all relevant coding SNPs.
However, for each donor, only heterozygous coding SNPs are relevant when measuring ASE effects.
'''
import sys
import vcf
import argparse
import pandas as pd

def parse_args(argv):
    p = argparse.ArgumentParser(description='This scripts filters all heterozygous coding SNPs for a given donor.')
    p.add_argument('--vcf', required = True)
    p.add_argument('--counts', required = True)
    p.add_argument('-n', '--genotype-name', required = True)
    p.add_argument('-o', '--out-csv', required = True)
    p.add_argument('--lim-reads', type=int, default = 30, help='Minimum required reads per SNP')
    return p.parse_args(argv)

def test_argparse():
    s = '--vcf tmp_data/unambigous_coding_snps.vcf.bgz --counts out/ase_counts/batch2__TCC.03.1__t0.csv -n TCC-03 -o /tmp/test.csv'.split()
    args = parse_args(s)
    return args

def parse_vcf(vcf_path, is_compressed):
    return vcf.Reader(filename=vcf_path, compressed=is_compressed)

def get_het_snps(vcf, sample):
    return [snp.ID for snp in vcf if snp.genotype(sample).is_het]

def parse_counts(tsv_path):
    df = pd.read_table(tsv_path)
    return df

def main(args):
    vcf = parse_vcf(args.vcf, is_compressed=args.vcf.endswith('gz'))
    snps = list(vcf)

    het_snps = get_het_snps(snps, args.genotype_name)
    counts = parse_counts(args.counts)
    is_heterozygous = counts.variantID.isin(het_snps)
    has_reads = counts.totalCount >= args.lim_reads
    het_counts = counts[is_heterozygous & has_reads]
    het_counts.to_csv(args.out_csv, index=False)

#vcf_path = 'tmp_data/unambigous_coding_snps.vcf.bgz'
if __name__ == '__main__':

    args = parse_args(sys.argv[1:])
    main(args)