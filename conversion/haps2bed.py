#!/bin/python3

## Tommy Carstensen
## Wellcome Trust Sanger Institute
## April 2014, August 2016

## http://pngu.mgh.harvard.edu/~purcell/plink/binary.shtml

import argparse
import contextlib
import fileinput
import sys
import os
import gzip

def main():

    d_args = argparser()

    magic_number = bytearray([108,27])
    mode = bytearray([1])

    with contextlib.ExitStack() as stack:
        if d_args['haps'] == [sys.stdin]:
            haps = stack.enter_context(d_args['haps'][0])
        else:
            for path in d_args['haps']:
                if not os.path.isfile(path):
                    print(path, 'does not exist')
                    exit()
            haps = stack.enter_context(fileinput.FileInput(
                files=d_args['haps'], openhook=hook_compressed_text))
        bed = stack.enter_context(open(d_args['bfile']+'.bed', 'wb'))
        bim = stack.enter_context(open(d_args['bfile']+'.bim', 'w'))
        bed.write(magic_number+mode)
        line = haps.readline()
        n_columns = len(line.split())
        n_samples = int((n_columns-5)/2)
        print(n_samples,'samples')
        if n_samples % 4:
            print('possibly need to write more code')
            sys.exit()
        while line:
            barray, line_bim = line2barray(line, n_columns)
            line = haps.readline()
            bed.write(bytes(barray))
            bim.write(line_bim)

    with open(d_args['sample']) as f:
        with open(d_args['bfile']+'.fam', 'w') as fam:
            for i in range(2):
                f.readline()
            for line in f:
                l = line.rstrip().split()
                fam.write(l[0], l[1], 0, 0, -9, -9, sep='\t')

    return


def hook_compressed_text(filename, mode):

    ##http://stackoverflow.com/questions/21529163/python-gzipped-fileinput-returns-binary-string-instead-of-text-string/21529243

    print(filename)
    ext = os.path.splitext(filename)[1]
    if ext == '.gz':
        f = gzip.open(filename, mode + 't')
    else:
        f = open(filename, mode)

    return f


def line2barray(line, n_columns):

    l = line.split()

    try:
        chrom = l[0].split(':')[0]
    except:
        chrom = l[0]
        print(line)
        stop
    rsID = l[1]
    pos = l[2]
    alleleA = l[3]
    alleleB = l[4]
    line_bim = '%s\t%s\t%s\t%s\t%s\t%s\n' %(
#        chrom, rsID, '0', pos, alleleA, alleleB)
        chrom, chrom+':'+pos+':'+alleleA+":"+alleleB, '0', pos,
        alleleA, alleleB)

    b = ''
    barray = bytearray()
    for i in range(5, n_columns, 2):
        if l[i:i+2] == ['0','0']:
            b += '00'
        elif l[i:i+2] == ['1','1']:
            b += '11'
        elif l[i:i+2] == ['0','1']:
            b += '01'
        elif l[i:i+2] == ['1','0']:
            b += '01'
        else:
            print(l[i:i+2])
            stop
        if len(b) == 8:
            print('bbb', i, b)
            barray.append(int(b[::-1],2))
            b = ''
##    if len(b) > 0:
    if b:
        barray.append(int(b[::-1].zfill(8),2))

    return barray, line_bim


def argparser():

    parser = argparse.ArgumentParser()

    parser.add_argument('--haps', '--in', nargs='+', default=[sys.stdin])
    parser.add_argument('--sample', required=True)
    parser.add_argument('--bfile', '--out', default='haps2bed')

    namespace_args = parser.parse_args()

    d_args = vars(namespace_args)

    return d_args


if __name__ == '__main__':
    main()
