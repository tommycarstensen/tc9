#!/bin/env python3

## Tommy Carstensen
## Wellcome Trust Sanger Institute
## September 2015

import argparse
import pysam
import re
import gzip
import sys
import urllib.request
import natsort


def main():

    args = parse_args()

    d_samples = parse_samples(args.samples)

    for vcf in natsort.natsorted(args.vcf):

        ## Read tabix indexed VCF file.
        tbx = pysam.TabixFile(args.vcf)

        ## Loop metadata lines and header line.
        for line in tbx.header:
            pass
        ## Parse header.
        header = line.decode().split('\t')
        ## Get sample indices for each population.
        for pop in d_samples.keys():
            ## Remove samples not in VCF file from dict to avoid ValueError.
            for sampleID in list(d_samples[pop]):
                if not sampleID in header:
                    d_samples[pop].remove(sampleID)
            ## Get index for each sample in the population.
            d_samples[pop] = set([
                header.index(sampleID) for sampleID in d_samples[pop]])
            continue

    loop_records(args, d_samples, tbx)

    return


def loop_records(args, d_samples, tbx):

    ## Print to header.
    print('#Formatting is allele count of REF,ALT1,...,ALTn,nonmiss,miss alleles per population.')

    ## Sort populations alphabetically.
    l_pops = list(sorted(d_samples.keys()))
    ## Print populations to header.
    print('#CHROM', 'POS', 'ID', 'REF', 'ALT', sep='\t', end='')
    for pop in l_pops:
        print('\t', pop, end='')
    print('')

    ## Compile re pattern.
    pattern = re.compile('[/|]')
    ## Loop over records.
    for record in tbx.fetch(args.chrom, args.start, args.end):
        l = record.split()
        ## Parse REF and ALT from the split line.
        CHROM = l[0]
        POS = l[1]
        ID = l[2]
        REF = l[3]
        ALT = l[4]
        ## Count of allele types.
        n = len(ALT.split(',')) + 1
        print(CHROM, POS, ID, REF, ALT, sep='\t', end='')
        ## Loop over sorted populations.
        for pop in l_pops:
            print('\t', end='')
            d_GT = {}
            nonmiss = 0
            miss = 0
            for index in d_samples[pop]:
                GT = re.split(pattern, l[index].split(':')[0])
                if GT == ['.', '.']:
                    miss += 2
                    continue
                nonmiss += 2
                for _ in map(int, GT):
                    try:
                        d_GT[_] += 1
                    except KeyError:
                        d_GT[_] = 1
                    continue
                continue
            for i in range(n):
                try:
                    print(d_GT[i], end=',')
                except KeyError:
                    print(0, end=',')
            print(nonmiss, miss, end='', sep=',')
        print('')

    return


def parse_samples(uri):

    ## Build dict with k=pop and v=set of samples.
    d_samples = {}
    f = open_uri(uri)
    for line in f:
        ## Skip blank lines.
        if not line.rstrip():
            continue
        sampleID, pop = line.rstrip().split()[:2]
        try:
            d_samples[pop].add(sampleID)
        except KeyError:
            d_samples[pop] = set([sampleID])

    return d_samples


def open_uri(uri):

    if uri.startswith('ftp://'):
        req = urllib.request.Request(uri)
        with urllib.request.urlopen(req) as response:
            f = response.read().decode('utf-8').split('\n')
            return f
    else:
        if uri.endswith('.gz'):
            with gzip.open(uri, 'rt') as f:
                return f
        else:
            with open(uri) as f:
                return f


def parse_args():

    parser = argparse.ArgumentParser()

    parser.add_argument('--vcf', required=True, nargs='+')
    parser.add_argument(
        '--samples', required=True,
        help='columns: 1=sampleID, 2=pop; e.g. ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel')
    parser.add_argument('--chrom', type=str, default=None)
    parser.add_argument('--start', type=int, default=None)
    parser.add_argument('--end', type=int, default=None)
##    parser.add_argument('--out', required=False, default=sys.stdout)

    args = parser.parse_args()

    return args


if __name__ == '__main__':
    main()
