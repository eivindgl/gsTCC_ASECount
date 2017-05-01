import collections
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


samples = list()
with open('input_data/yang_sample_names.txt') as f:
    for line in f:
        x = line.strip().split()[-1].strip('"')
        y = Sample.parse_string(x)
        samples.append(y)

with open('input_data/bamfilenames.txt') as f:
    xs = [x.strip() for x in f]
    p = regex.compile(r'.*_t\d+_')
    ys = [p.search(x).captures()[0][:-1] for x in xs]
    bamfiles = dict(zip(ys, xs))
    name = ys[0]

with open('tmp_data/genotype_names.txt') as f:
    genotype_names = [x.strip() for x in f]

import csv
with open('tmp_data/run_meta.csv', 'w') as f:

    xs = [collections.OrderedDict(
            sample_name = sample.sample_name,
            bam = bamfiles[sample],
            genotype_name = sample.donor)
    for sample in samples]
    writer = csv.DictWriter(f, xs[0].keys())
    writer.writeheader()
    writer.writerows(xs)


