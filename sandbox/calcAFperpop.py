#!/bin/env python3

## Tommy Carstensen
## Wellcome Trust Sanger Institute
## September 2015, September 2016

import argparse
import pysam
import re
import gzip
import sys
import urllib.request
import natsort
import contextlib


def main():

    args = parse_args()

    d_samples = parse_samples(args.samples)

    for vcf in natsort.natsorted(args.vcf):

        ## Read tabix indexed VCF file.
        tbx = pysam.TabixFile(vcf)

        ## Loop metadata lines and header line.
        for line in tbx.header:
            ## Print metadata lines.
            if line[:2] == b'##':
                print(line.decode())

        ## Parse header after looping over metadata lines.
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
        n_samples = len(header)-9

        ## Append new metadata lines.
        for abbreviation, text in (
            ('AF', 'Allele Frequency'),
            ('AC', 'Allele Count'),
            ):
            print('##INFO=<ID={},Number=A,Type=Float,'.format(
                abbreviation), end='')
            print('Description="{}">'.format(text))
            for pop in sorted(d_samples.keys()):
                print('##INFO=<ID={}_{},Number=A,Type=Float,'.format(
                    abbreviation, pop), end='')
                print('Description="{}'.format(text), end=' ')
                print('in the {} population with {} samples">'.format(
                    pop, len(d_samples[pop])))

        ## Print header line with sample IDs.
        if args.bool_keep_GTs:
            print(line.decode())
        ## Print header line without sample IDs.
        else:
            print('\t'.join(line.decode().split('\t')[:9]))

    loop_records(args, d_samples, tbx, n_samples)

    return


def loop_records(args, d_samples, tbx, n_samples):

    assert 'all' not in d_samples.keys()

    ## Sort populations alphabetically.
    l_pops = ['all']+list(sorted(d_samples.keys()))

    d_samples['all'] = range(9, n_samples+9)

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
        QUAL = l[5]
        FILTER = l[6]
        INFO = l[7]
        FORMAT = l[8]
        ## Count of allele types.
        n = len(ALT.split(',')) + 1
        print(
            CHROM, POS, ID, REF, ALT, QUAL, FILTER, sep='\t', end='\t')
        d_AFs = {}
        d_ACs = {}
        ## Loop over sorted populations.
        for pop in l_pops:
##            print('\t', end='')
            d_GT = {i: 0 for i in range(n)}
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
            d_AFs[pop] = ','.join(['0' if d_GT[i] == 0 else str(round(
                d_GT[i] / (
                    2 * len(d_samples[pop])), 5)) for i in range(1, n)])
            d_ACs[pop] = ','.join([str(d_GT[i]) for i in range(1, n)])
            ## Continue loop over pops.
            continue

        ## Overwrite or append to INFO string.
        if args.bool_overwrite_INFO:
            INFO = ''
        else:
            INFO += ';'
        for pop in l_pops:
            if pop == 'all':
                k = ''
            else:
                k = '_{}'.format(pop)
                INFO += ';'
            INFO += 'AF{}={}'.format(k, d_AFs[pop])
            INFO += ';AC{}={}'.format(k, d_ACs[pop])
        print(INFO, FORMAT, sep='\t', end='')
        ## Append GTs to the line.
        if args.bool_keep_GTs:
            print('\t', '\t'.join(l[9:]), sep='', end='')
        print('')

    return


def parse_annotation(vcf, args, chrom, start):

    tbx = pysam.TabixFile(vcf)
    for record in tbx.fetch(chrom, start-1, None):
        l = record.split('\t')
        INFO = l[7]
        POS = l[1]
        yield INFO, POS


def parse_samples(uri):

    ## Build dict with k=pop and v=set of samples.
    d_samples = {}
    with open_uri(uri) as f:
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


@contextlib.contextmanager
def open_uri(uri):

    if uri.startswith('ftp://'):
        req = urllib.request.Request(uri)
        with urllib.request.urlopen(req) as response:
            f = response.read().decode('utf-8').split('\n')
            yield f
    else:
        if uri.endswith('.gz'):
            with gzip.open(uri, 'rt') as f:
                yield f
        else:
            with open(uri) as f:
                yield f


def parse_args():

    parser = argparse.ArgumentParser()

    parser.add_argument('--vcf', required=True, nargs='+')
    parser.add_argument(
        '--samples', required=True,
        help='columns: 1=sampleID, 2=pop; e.g. \
ftp://ftp.1000genomes.ebi.ac.uk/vol1/\
ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel')
    parser.add_argument('--chrom', type=str, default=None)
    parser.add_argument('--start', type=int, default=None)
    parser.add_argument('--end', type=int, default=None)
##    parser.add_argument('--out', required=False, default=sys.stdout)
    parser.add_argument('--bool_keep_GTs', action='store_true')
    parser.add_argument(
        '--bool_overwrite_INFO', action='store_true',
        help='Overwrite existing INFO field.')

    args = parser.parse_args()

    return args


if __name__ == '__main__':
    main()
