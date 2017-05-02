#!/usr/bin/env python3
import sys
import vcf
import argparse
import pandas as pd

def parse_args(argv):
    p = argparse.ArgumentParser()
    p.add_argument('--vcf', required = True)
    p.add_argument('--counts', required = True)
    p.add_argument('-n', '--genotype-name', required = True)
    p.add_argument('-o', '--out-csv', required = True)
    return p.parse_args(argv)

def test_argparse():
    s = '--vcf tmp_data/unambigous_coding_snps.vcf.bgz --counts out/ase_counts/batch2__TCC.03.1__t0.ASE.csv -n TCC-03 -o /tmp/test.csv'.split()
    args = parse_args(s)

def parse_vcf(vcf_path, is_compressed):
    return vcf.Reader(filename=vcf_path, compressed=is_compressed)

def get_het_snps(vcf, sample):
    return [snp.ID for snp in vcf if snp.genotype(sample).is_het]

def parse_counts(tsv_path):
    df = pd.read_table(tsv_path)
    return df

def main(args):
    vcf = parse_vcf(args.vcf, is_compressed=args.vcf.endswith('gz'))
    genotype_names = vcf.samples
    snps = list(vcf)

    het_snps = get_het_snps(snps, args.genotype_name)
    counts = parse_counts(args.counts)
    het_counts = counts[counts.variantID.isin(het_snps)]
    het_counts.to_csv(args.out_csv, index=False)

#vcf_path = 'tmp_data/unambigous_coding_snps.vcf.bgz'
if __name__ == '__main__':

    args = parse_args(sys.argv[1:])
    main(args)