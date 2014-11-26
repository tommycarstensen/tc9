#!/bin/python3

## Tommy Carstensen, Wellcome Trust Sanger Institute, September 2014

import argparse
import gzip
import contextlib
import os
import collections
import re


def main():

    d_args = argparser()

##    with contextlib.ExitStack() as stack:
##        vcf = stack.enter_context(open_file(d_args['vcf']))
##        fam = stack.enter_context(open('{}.fam'.format(d_args['bed'])))
##        for line in vcf:
##            if line[0] != '#':
##                break
##        print(line)
##        stop1

    ## sample removal
    if d_args['remove']:
        with open(d_args['remove']) as f:
            set_samples_remove = set(f.read().splitlines())

    with contextlib.ExitStack() as stack:
        bed = stack.enter_context(open('{}.bed'.format(d_args['bed']), 'wb'))
        bim = stack.enter_context(open('{}.bim'.format(d_args['bed']), 'w'))
        fam = stack.enter_context(open('{}.fam'.format(d_args['bed']), 'w'))
        if d_args['extract']:
            extract = stack.enter_context(open(d_args['extract']))

        magic_number = bytearray([108,27])
        mode = bytearray([1])
        bed.write(magic_number+mode)

        if d_args['extract']:
            chrom_e, pos_e, ref_e, alt_e, ID_e = next(parse_extract(extract))
        ## Loop over sorted VCF files.
        for i_vcf, filename_vcf in enumerate(sort_nicely(d_args['vcf'])):
            vcf = stack.enter_context(open_file(filename_vcf))
            ## Determine which columns to parse.
            for line in vcf:
                if line[:2] == '##':
                    continue
                if i_vcf == 0:
                    l_columns = []
                    for i_sample, sample in enumerate(
                        line.rstrip().split('\t')[9:], 9):
                        if sample in set_samples_remove:
                            continue
                        l_columns.append(i_sample)
                        fam.write('{ID} {ID} 0 0 -9 -9\n'.format(ID=sample))
                        ## Do not parse genotypes yet.
                break
            ## Parse first variant.
            try:
                l, chrom, pos, ref, alt = next(parse_vcf(vcf))
            ## Skip empty VCF file (e.g. near centromere)
            except StopIteration:
                continue
            while True:
                if d_args['extract']:
                    if int(chrom) < int(chrom_e):
                        try:
                            l, chrom, pos, ref, alt = next(parse_vcf(vcf))
                        except StopIteration:
                            break
                        continue
                    elif int(chrom) > int(chrom_e):
                        chrom_e, pos_e, ref_e, alt_e, ID_e = next(parse_extract(extract))
                        continue
                    if pos < pos_e:
                        try:
                            l, chrom, pos, ref, alt = next(parse_vcf(vcf))
                        except StopIteration:
                            break
                        continue
                    elif pos > pos_e:
                        chrom_e, pos_e, ref_e, alt_e, ID_e = next(parse_extract(extract))
                        continue
                    else:
                        bool_continue = False
                        ## ref=major
                        if ref == ref_e and alt == alt_e:
                            pass
                        ## ref=minor
                        elif ref == alt_e and ref_e == alt:
                            pass
                        ## monomophic on chip
                        elif ref == alt_e and ref_e == '0' and len(alt) == 1:
                            bool_continue = True
                        ## monomophic on chip
                        elif alt == alt_e and ref_e == '0' and len(ref) == 1:
                            bool_continue = True
                        ## INDEL
                        elif any(len(GT) > 1 for GT in ref.split(',')):
                            bool_continue = True
                        ## INDEL
                        elif any(len(GT) > 1 for GT in alt.split(',')):
                            bool_continue = True
                        ## triallelic
                        elif ref == ref_e and alt_e in alt.split(',') and all(len(GT) == 1 for GT in alt.split(',')):
                            bool_continue = True
                        ## triallelic
                        elif ref == ref_e and alt_e not in alt.split(',') and all(len(GT) == 1 for GT in alt.split(',')):
                            bool_continue = True
                        ## triallelic
                        elif ref == alt_e and ref_e in alt.split(',') and all(len(GT) == 1 for GT in alt.split(',')):
                            bool_continue = True
                        ## triallelic
                        elif ref == alt_e and ref_e not in alt.split(',') and all(len(GT) == 1 for GT in alt.split(',')):
                            bool_continue = True
                        ## triallelic
                        elif ref == alt_e and ref_e == '0' and all(len(GT) == 1 for GT in alt.split(',')):
                            bool_continue = True
                        ## triallelic
                        elif ref_e == '0' and alt_e in alt.split(',') and all(len(GT) == 1 for GT in alt.split(',')):
                            bool_continue = True
                        ## allele difference
                        elif ref == ref_e and alt != alt_e and len(alt) == 1:
                            bool_continue = True
                        ## allele difference
                        elif ref == alt_e and alt != ref_e and len(alt) == 1:
                            bool_continue = True
                        ## allele difference
                        elif ref != ref_e and alt != alt_e and ref_e == '0':
                            bool_continue = True
                        else:
                            print(chrom, pos, ref, alt)
                            print(chrom_e, pos_e, ref_e, alt_e)
                            bool_continue = True
                        if bool_continue:
                            print('skip', chrom, pos, ref, alt, ref_e, alt_e)
                            l, chrom, pos, ref, alt = next(
                                parse_vcf(vcf))
                            try:
                                chrom_e, pos_e, ref_e, alt_e, ID_e = next(
                                    parse_extract(extract))
                            except StopIteration:
                                break
                            continue
                        pass
                    pass
                ## write bim
                bim.write(
                    '{chrom}\t{ID}\t0\t{pos}\t{major}\t{minor}\n'.format(
                        chrom=chrom, ID=ID_e, pos=pos, major=ref, minor=alt))
                ## write bed
                b = ''
                barray = bytearray()
                for column in l_columns:
                    GT = l[column].split(':',1)[0]
                    ## HOMREF
                    if GT == '0/0':
                        b += '00'
                    ## HOMALT
                    elif GT == '1/1':
                        b += '11'
                    ## HET
                    elif GT == '0/1':
                        b += '01'
                    ## unknown
                    elif GT == './.':
                        b += '10'
                    else:
                        print(chrom, pos, ref, alt, ref_e, alt_e)
                        print(l[column])
                        stopGT
                    if len(b) == 8:
                        barray.append(int(b[::-1],2))
                        b = ''
                    ## Continue loop over samples.
                    continue
                if b:
                    barray.append(int(b[::-1].zfill(8),2))
                ## write bed
                bed.write(bytes(barray))
                ## Read next line of vcf.
                l, chrom, pos, ref, alt = next(parse_vcf(vcf))
                chrom_e, pos_e, ref_e, alt_e, ID_e = next(parse_extract(extract))
                if pos % 100 == 0:
                    print(chrom, pos)
                ## Continue while loop.
                continue
            ## Continue loop over vcf files.
            continue
        ## Exit contextlib stack.
        pass

    check_output(d_args)

    return


def check_output(d_args):

    with open(d_args['bed']+'.fam') as f:
        sample = len(f.readlines())
    with open(d_args['bed']+'.bim') as f:
        for SNP, line in enumerate(f, 1):
            continue
    print('size', os.path.getsize(d_args['bed']+'.bed'))
    print('SNP', SNP, 'sample', sample, (sample+sample%4)/4)
    assert os.path.getsize(d_args['bed']+'.bed') == SNP*(sample+sample%4)/4+3

    return


def alphanum_key(s):
    ## http://dave.st.germa.in/blog/2007/12/11/exception-handling-slow/
    NUM_RE = re.compile('([0-9]+)')
    return [ int(c) if c.isdigit() else c for c in NUM_RE.split(s) ]


def sort_nicely(l):
    ## http://nedbatchelder.com/blog/200712/human_sorting.html
    """ Sort the given list in the way that humans expect.
    """
    l.sort(key=alphanum_key)
    return l


def parse_extract(extract):

    for line in extract:
        l = line.rstrip().split()
        chrom = l[0]
        ID = l[1]
        pos = int(l[3])
        ref = l[4]
        alt = l[5]
        yield chrom, pos, ref, alt, ID


def parse_vcf(vcf):

    for line in vcf:
        l = line.rstrip().split()
        chrom = l[0]
        pos = int(l[1])
        ref = l[3]
        alt = l[4]
        yield l, chrom, pos, ref, alt


def open_file(file_name):

    if os.path.splitext(file_name)[-1] == '.gz':
        fd = gzip.open(file_name, 'rt')
    else:
        fd = open(file_name)

    return fd


def argparser():

    parser = argparse.ArgumentParser()

    parser.add_argument('--vcf', required=True, help='vcf input file', nargs='+')
    parser.add_argument(
        '--bed', '--bfile', required=True,
        help='bed output file')
    parser.add_argument(
        '--extract', required=False,
        help='SNPs to extract (<CHROM>\t<POS>\t<REF>\t<ALT>)')
    parser.add_argument(
        '--remove', required=False,
        help='samples to remove (<ID>)')
##    parser.add_argument(
##        '--merge', required=False,
##        help='bed file to merge vcf file with')

    d_args = vars(parser.parse_args())

    return d_args


if __name__ == '__main__':
    main()


### format in (or out)
##http://pngu.mgh.harvard.edu/~purcell/plink/binary.shtml
##
### format out
##http://stackoverflow.com/questions/2872381/how-to-read-a-file-byte-by-byte-in-python-and-how-to-print-a-bytelist-as-a-binar
##
### chunk
##http://stackoverflow.com/questions/1035340/reading-binary-file-in-python
##
### buffer size
##http://stackoverflow.com/questions/236861/how-do-you-determine-the-ideal-buffer-size-when-using-fileinputstream
