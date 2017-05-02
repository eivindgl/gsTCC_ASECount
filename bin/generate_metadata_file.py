#!/usr/bin/env python3
'''
Parses a set a samplenames and bamfiles with specific format conventions and generates a metadata csv file.
So in essence, makes conventional information explicit.
'''


import argparse
import collections
import csv

import regex

class Sample:

    def __init__(self, batch, donor_id, replicate, timepoint, sample_name):
        self.batch = batch
        self.donor_id = donor_id
        self.replicate = replicate
        self.timepoint = timepoint
        self.sample_name = sample_name
        self.idstring = f'{batch}_TCC-{donor_id}-{replicate}_{timepoint}'

    @property
    def donor(self):
        return f'TCC-{self.donor_id}'

    def __hash__(self):
        return hash(str(self))

    def __str__(self):
        return self.idstring

    def __eq__(self, other):
        return str(self) == str(other)

    @classmethod
    def parse_string(cls, sample_name):
        x = regex.match(r'(.*)__TCC\.(\d+)\.(\d+)__(.+)', sample_name).groups()
        return cls(*x, sample_name = sample_name)


def parse_sample_to_bamfile_map(filenames):
    p = regex.compile(r'.*_t\d+_')
    sample_names = [p.search(filename).captures()[0][:-1] for filename in filenames]
    return dict(zip(sample_names, filenames))

def main(args):
    samples = [Sample.parse_string(line.strip()) for line in args.sample_names]
    bamfile_map = parse_sample_to_bamfile_map([x.strip() for x in args.bamfiles.readlines()])

    sample_list = [collections.OrderedDict(
            sample_name = sample.sample_name,
            bam = bamfile_map[sample],
            genotype_name = sample.donor,
            timepoint = sample.timepoint)
        for sample in samples]
    writer = csv.DictWriter(args.out, fieldnames=sample_list[0].keys())
    writer.writeheader()
    writer.writerows(sample_list)

if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('sample_names', type = argparse.FileType('r'))
    p.add_argument('bamfiles', type=argparse.FileType('r'))
    p.add_argument('-o', '--out', type=argparse.FileType('w'), default = '-')
    main(p.parse_args())



