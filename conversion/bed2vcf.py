#!/bin/python3

## Tommy Carstensen, Wellcome Trust Sanger Institute
## December 2013, January 2015

import argparse
import os
import math
import re
import datetime
import sys
import gzip
import Bio
from Bio.bgzf import BgzfWriter


bpb = n_bits_per_byte = 8


def main():

    args = argparser()

    l_fam_sampleIDs = parse_fam(args.bfile + '.fam')
    l_vcf_sampleIDs = shorten_sampleIDs(l_fam_sampleIDs)

    if args.keep:
        l_keep_sampleIDs = parse_keep(path_keep, l_fam_sampleIDs)
        l_keep_index = index_keep(l_keep_sampleIDs, l_fam_sampleIDs)
    else:
        l_keep_index = list(range(len(l_fam_sampleIDs)))

    assert len(l_vcf_sampleIDs) == len(l_fam_sampleIDs)

    ## Get count of samples.
    n_samples = len(l_fam_sampleIDs)
    ## Get count of SNPs. Only needed to keep track of completion percent.
    with open(args.bfile + '.bim', 'r') as bim:
        n_SNPs = len(bim.readlines())

    d_fai = read_fai(args.ref + '.fai')

##         gzip.open(args.vcf, 'wt') as vcf, \
    with open(args.bfile + '.bed', 'rb') as bed, \
         BgzfWriter(args.vcf, 'wb') as vcf, \
         open(args.bfile + '.bim', 'r') as bim, \
         open(args.ref, 'r') as ref:
        convert(
            args, bed, vcf, bim, ref, d_fai,
            l_vcf_sampleIDs, l_keep_index, n_samples, n_SNPs)

    return

def convert(
    args, bed, vcf, bim, ref, d_fai,
    l_vcf_sampleIDs, l_keep_index, n_samples, n_SNPs):

    n_bytes_per_SNP = math.ceil(n_samples/4)

    ## Prepare list of genotypes instead of appending to empty list.
    ## The latter is slow.
    l_GT = [None]*len(l_keep_index)

    write_metadata(args, vcf)
    write_header(vcf, l_vcf_sampleIDs)

    ## Write first 3 bytes of bed file.
    magic_number = bytearray([108,27])
    mode = bytearray([1])
    bed.read(len(magic_number)+len(mode))

    QUAL = '.'
    FILTER = 'PASS'
    FORMAT = 'GT'

    ## data lines
    for i_SNP, line_bim in enumerate(line_bim):
        i_SNP += 1
        ## By default, the minor allele is coded A1
        ## and the major allele is coded A2
        CHROM, ID, _, POS, A1, A2 = line_bim.rstrip().split()
        if args.chrom and CHROM != args.chrom:
            continue
        if i_SNP % 1000 == 0:
            print('SNP {} of {} SNPs. CHROM={} POS={} time={}'.format(
                i_SNP, n_SNPs, CHROM, POS,
                datetime.datetime.now().strftime("%H:%M:%S")))

        ## Parse reference allele from reference sequence.
        REF = parse_ref(d_fai, ref, CHROM, int(POS))

        ## Parse alleles and genotypes from bed file.
        cnt_allele_nonmissing, cnt_allele_A1, genotype_fields = parse_bed(
            bed, l_keep_index, REF, A2, n_bytes_per_SNP, l_GT)
##        NCHROBS = cnt_allele_nonmissing

        ## Calculate minor allele frequency.
        AF_A1 = cnt_allele_A1/cnt_allele_nonmissing

        ## ALT = major
        ## Convert MAF to allele frequency for the ALT allele.
        if REF == A1:
            ALT = A2
            AF = 1-AF_A1
        ## ALT = minor
        elif REF == A2:
            ALT = A1
            AF = AF_A1
        elif A1 == '0':
            ## 
            if not AF_A1 == 0:
                print('WARNING: REF={} CHROM={} POS={} ID={} A1={} A2={} AF={:.4f}'.format(
                    REF, CHROM, POS, ID, A1, A2, AF_A1), file=sys.stderr)
                ALT = A2
                AF = 1-AF_A1
                continue
            ## monomorphic
            else:
                ALT = '.'
                AF = 0
        ## Neither A1 nor A2 is identical to the reference allele.
        ## Should not happen. Print an error and continue.
        else:
            print('ERROR: REF={} CHROM={} POS={} ID={} A1={} A2={} AF={:.4f}'.format(
                REF, CHROM, POS, ID, A1, A2, AF_A1), file=sys.stderr)
            continue

        ## Join fixed fields.
        INFO = 'AF={0:.4f}'.format(AF)
        fixed_fields = '\t'.join((
            CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT))

        ## Combine fixed fields and genotype fields.
        line_vcf = fixed_fields+'\t'+genotype_fields+'\n'

        ## Write line to file.
        vcf.write(line_vcf)

    return


def parse_bed(bed, l_keep_index, REF, A2, n_bytes_per_SNP, l_GT):

    '''http://pngu.mgh.harvard.edu/~purcell/plink/binary.shtml'''

    bytesSNP = bed.read(n_bytes_per_SNP)
    ## Convert to integer from bytes.
##    int_bytesSNP = int.from_bytes(bytesSNP, 'big')
    int_bytesSNP2 = int.from_bytes(bytesSNP, 'little')
    cnt_allele_A1 = 0
    cnt_allele_nonmissing = 0
    for i, i_keep in enumerate(l_keep_index):
        ## Calculate size of bit shift for sample.
        shift2 = 2*i_keep
        ## Query status of bit with bit masking
        bits = int_bytesSNP2 >> shift2 & 0b11
        ## query status of bit
        if bits == 0:
            ## hom
            s2b = '00'
            cnt_allele_A1 += 2
            cnt_allele_nonmissing += 2
            if REF == A2:
                GT = '1/1'
            else:
                GT = '0/0'
        elif bits == 2:
            ## het
            s2b = '10'
            cnt_allele_A1 += 1
            cnt_allele_nonmissing += 2
            GT = '0/1'
        elif bits == 3:
            ## hom
            s2b = '11'
            cnt_allele_nonmissing += 2
            if REF == A2:
                GT = '0/0'
            else:
                GT = '1/1'
        elif bits == 1:
            ## missing
            s2b = '01'
            GT = './.'
        else:
            print(s2, x)
            stop2
        l_GT[i] = GT

    genotype_fields = '\t'.join(l_GT)

    return cnt_allele_nonmissing, cnt_allele_A1, genotype_fields


def write_metadata(args, vcf):

    ## metadata
    vcf.write('##fileformat=VCFv4.3\n')
    vcf.write('##fileDate={}\n'.format(datetime.date.today()))
    vcf.write('##source={}.bed\n'.format(args.bfile))
    vcf.write('##source={}\n'.format(sys.argv[0]))
    vcf.write('##reference={}\n'.format(args.ref))
    vcf.write('##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">\n')
    vcf.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')

    return

def write_header(vcf, l_vcf_sampleIDs):

    ## header line
    vcf.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT')
    for sampleID in l_vcf_sampleIDs:
        vcf.write('\t{}'.format(sampleID))
    vcf.write('\n')

    return

def read_fai(path_fai):

    assert os.path.isfile(path_fai)

    d = {}
    with open(path_fai) as file_fai:
        for line_fai in file_fai:
            l = line_fai.rstrip().split()
            chrom = l[0]
            byte_length = int(l[1])
            byte_start = int(l[2])
            bytes_per_line_excl_line_break = int(l[3])
            bytes_per_line_incl_line_break = int(l[4])
            d[chrom] = {'length': byte_length, 'start': byte_start}

    return d


def parse_ref(d_fai, fd_ref, CHROM, POS, size=1):

    ## parse ref
    cnt = cnt_chars_per_line_excluding_newline = 60
    row1 = (POS - 1) // cnt
    row2 = (POS - 1 + size) // cnt
    size += row2 - row1
    col = (POS - 1) % cnt
    byte_init = d_fai[CHROM]['start']
    offset = byte_init + (cnt + 1) * row1 + col
    fd_ref.seek(offset)
    read = fd_ref.read(size).replace('\n', '')

    return read


def index_keep(l_keep_sampleIDs,l_fam_sampleIDs):

    if not l_keep_sampleIDs:
        l_keep_sampleIDs = l_fam_sampleIDs

    if len(set(l_keep_sampleIDs)-set(l_fam_sampleIDs)) > 0:
        print('keep is not a subset of fam')
        sys.exit()

##    l_keep_index = []
##    for i in range(len(l_fam_sampleIDs)):
##        sampleID = l_fam_sampleIDs[i]
##        if not sampleID in l_keep_sampleIDs:
##            continue
##        l_keep_index += [i]
    l_keep_index = list(sorted(
        [l_fam_sampleIDs.index(sampleID) for sampleID in l_keep_sampleIDs]))

    return l_keep_index
    

def parse_fam(fp_fam,):
    
    l_keep_index = []

    with open(fp_fam, 'r') as fam:
        l_sampleIDs = [line.rstrip().split()[0] for line in fam]
    with open(fp_fam, 'r') as fam:
        line = fam.readline()
        if line.split()[0] != line.rstrip().split()[1]:
            print('havent written code for when FID and IID are different')
            sys.exit()

    if len(l_sampleIDs) != len(set(l_sampleIDs)):
        print('duplicate sample IDs')
        print(l_sampleIDs)
        stop

    return l_sampleIDs


def shorten_sampleIDs(l_sampleIDs_long):

##    if fp_update_ids:
##        with open(fp_update_ids) as f:
##            d_update = {
##                line.strip().split()[0]:line.strip().split()[2] for line in f}
##    else:
##        d_update = {sampleID:sampleID for sampleID in l_sampleIDs_long}

    keyword = re.compile(r'(\d\d\d\d\d\d_[A-H]\d\d_)(.+\d\d\d\d\d\d\d)')

    l_sampleIDs_short = []

    for sampleID_long in l_sampleIDs_long:
##        sampleID_long = d_update[sampleID_long]
        match = result = keyword.search(sampleID_long)
        if match:
            sampleID_short = match.group(2)
        else:
            sampleID_short = sampleID_long
        l_sampleIDs_short += [sampleID_short]

    return l_sampleIDs_short


def parse_keep(fp_keep,l_fam_sampleIDs):

    if fp_keep:
        with open(fp_keep) as keep:
            l_keep_sampleIDs_unsorted = [line.rstrip().split()[0] for line in keep]
##        for sampleID in l_fam_sampleIDs:
##            if not sampleID in l_keep_sampleIDs_unsorted:
##                continue
        l_keep_sampleIDs = [
            sampleID for sampleID in l_fam_sampleIDs
            if sampleID in l_keep_sampleIDs_unsorted]
    else:
        l_keep_sampleIDs = l_fam_sampleIDs

    return l_keep_sampleIDs


def argparser():

    parser = argparse.ArgumentParser()

    ## Required.
    parser.add_argument('--vcf', '--out', required = True)
    parser.add_argument('--bfile', '--in', required = True)
    s_help = 'The A2 allele is saved as the reference by PLINK2, '
    s_help = 'but a reference sequence is needed to determine REF and ALT.'
    parser.add_argument(
        '--ref', required = True, default=None,
        help=s_help,
        )

    ## Optional. Good to have instead of creating intermediate bed files.
    parser.add_argument('--keep', required = False, default=None)
##    parser.add_argument('--update-ids', required = False, default=None)
    parser.add_argument('--chrom', required = False)

    args = namespace_args = parser.parse_args()

    ## Assert that input exists.
    assert all([
        os.path.isfile('{}.{}'.format(args.bfile, suffix))
        for suffix in ('bed','bim','fam')])
    assert os.path.isfile(args.ref)

    ## Do not allow non-compressed output vcf file.
    assert os.path.splitext(args.vcf)[1] == '.gz'

    ## Do not allow compressed reference sequence for now.
    assert os.path.splitext(args.ref)[1] == '.fa'

    return args


if __name__ == '__main__':
    main()
