#!/bin/python3

## Tommy Carstensen (tc9)
## Wellcome Trust Sanger Institute, June 2013

## built-ins
import glob
import re
import fileinput
import argparse
##import math
####sys.path.append('/nfs/users/nfs_t/tc9/github/sandbox')
####import gnuplot

'''This is a tool for extracting/excluding chromosomal positions
from one file (file 2)
from another file (file 1)
and converting to another file format (file3),
when the files are sorted by chromosomal position (human sort instead of machine sort).'''

def main():

    d = parse_args()

##    column2 = 1
##    column1 = 1
##    separator1 = ' '
##    separator2 = ' '
##    bool_header2 = True
##    bool_header1 = False
##    header_char1 = '#'
##    header_length1 = 1
##    format1 = 'tped'

##    format1 = 'gprobs'
##    format2 = 'prune.in'
##    format3 = 'EIGENSOFT'
##    mode = 'extract'
##    file_chrom_seq = None
##    files1 = ['../pipeline/uganda4x/out_BEAGLE/*.gprobs']
####    files1 = ['tped/*.tped']
##    files2 = ['prune/*.prune.in']
##    affix3 = 'EIGENSOFT'

    l_chrom_seq = [str(chrom) for chrom in range(1,23)]
    affix3 = d['affix3']
    format1 = d['format1']
    format2 = d['format2']
    format3 = d['format3']
    files1 = d['files1']
    files2 = d['files2']

    files1 = file_str2list(files1)
    files2 = file_str2list(files2)
    if not (files1 and files2):
        sys.exit()

    loop(
        files1,files2,affix3,format1,format2,format3,l_chrom_seq,
        )
    
    return


def loop(files1,files2,affix3,format1,format2,format3,l_chrom_seq,):

    d_read = {
        'tped':generate_line_tped,
        'prune.in':generate_line_prune,
        'bgl':generate_line_BEAGLE,
        'gprobs':generate_line_BEAGLE,
        }
    d_write = {
        'tped':write_line_tped,
        'EIGENSOFT':write_line_EIGENSOFT,
        }
    d_convert = {
        'gprobs':{'EIGENSOFT':convert_BEAGLE2EIGENSOFT},
        }
    read1 = d_read[format1]
    read2 = d_read[format2]
    if format3 == format1:
        convert = return_line
    else:
        convert = d_convert[format1][format3]
    write = d_write[format3]

    if format3 == 'EIGENSOFT':
        file3 = '%s.eigenstratgeno' %(affix3)
        file3b = '%s.snp' %(affix3)
    else:
        print(format3)
        sotp

    if format3 == 'EIGENSOFT':
        if format1 == 'gprobs':
            with fileinput.input(files=files1) as f1, open('%s.ind' %(affix3),'w') as f3:
                for sampleID in f1.readline().strip().split()[3:-1:3]:
                    f3.write('\t%s\tU\tNA\n' %(sampleID))
        else:
            print(format1)
            stop

    n_intersection = 0

    with fileinput.input(files=files1) as f1, fileinput.input(files=files2) as f2, open(file3,'w') as f3, open(file3b,'w') as f3b:

        chrom1,pos1,l1,line1 = next(read1(f1))
        chrom2,pos2,l2,line2 = next(read2(f2))
        print(chrom1,chrom2,pos1,pos2)

        while True:

            if chrom1 != chrom2:
                if l_chrom_seq.index(chrom2) > l_chrom_seq.index(chrom1):
                    try:
                        chrom1,pos1,l1,line1 = next(read1(f1))
                    except StopIteration:
                        break
                else:
                    try:
                        chrom2,pos2,l2,line2 = next(read2(f2))
                    except StopIteration:
                        break
            else:
                if pos1 == pos2:
                    if format3 == format1:
                        write_line(f3,line1)
                    elif format3 == format2:
                        write_line(f3,line2)
                    elif format3 == 'EIGENSOFT':
                        write_line_EIGENSOFT(f3,convert,l1)
                        write_line_EIGENSOFTsnp(f3b,chrom1,pos1,l1,line1,)
                    else:
                        print(format3)
                        stop
                    n_intersection += 1
                    try:
                        chrom1,pos1,l1,line1 = next(read1(f1))
                        chrom2,pos2,l2,line2 = next(read2(f2))
                    except StopIteration:
                        break
                elif pos2 > pos1:
                    try:
                        chrom1,pos1,l1,line1 = next(read1(f1))
                    except StopIteration:
                        break
                else:
                    try:
                        chrom2,pos2,l2,line2 = next(read2(f2))
                    except StopIteration:
                        break

    print('positions at intersection:',n_intersection)

    return


def write_line_EIGENSOFTsnp(f,chrom,pos,l,line,):

    ## /nfs/team149/Software/EIG4.2/CONVERTF/example.snp

    rsID = l[0]
    alleleA = l[1]
    alleleB = l[2]
##    chrom = chrom.replace('X','23')
    cM = 0
    line = '\t%s\t%s\t%i\t%i\t%s %s\n' %(rsID,chrom,cM,pos,alleleA,alleleB,)
    f.write(line)

    return


def write_line(f,line,):

    f.write(line)

    return


def convert_BEAGLE2EIGENSOFT(l):

    line = ''
    alleleA = l[1]
    alleleB = l[2]
    x = 0
    for i in range(3,len(l),3):
        x += 1
        l_probs = [float(l[i+j]) for j in range(3)]
        line += '%i' %(l_probs.index(max(l_probs)))
    line += '\n'

    return line


def write_line_EIGENSOFT(f,convert,l,):

    line = convert(l)
    f.write(line)

    return


def write_line_tped():

    line_tped = '%s %s:%s 0 %s' %(chrom,chrom,pos,pos)

    return


def return_line(line):

    return line


def file_str2list(l_fp):

    if len(l_fp) == 1 and '*' in l_fp[0]:
        l_fp = glob.glob(l_fp[0])
    l_fp = sort_nicely(l_fp)

    return l_fp


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


def parse_samples_tfam(fp_samples,):

    ## convert "plate_well_sample" format to "sample" format
    ## i.e. remove info about plate and well
    keyword = re.compile(r'(\d\d\d\d\d\d_[A-H]\d\d_)(.+\d\d\d\d\d\d\d)')
    l_samples = []
    with open(fp_samples) as lines:
        for line in lines:
            sampleID = line.split()[1]
            match = result = keyword.search(sampleID)
            if match:
                l_samples += [match.group(2)]
            else:
                l_samples += [sampleID]

    return l_samples


def generate_line_BEAGLE(f):

    set_nt = set(['A','C','G','T',])

    for line in f:
        l = line.strip().split(' ')
        if not l[1] in set_nt: continue
        if not l[2] in set_nt: continue
        chrom,pos = l[0].split(':')
        pos = int(pos)

        yield chrom,pos,l,line


def generate_line_tped(f):

    for line in f:
        l = line.strip().split(' ')
        chrom = l[0]
        pos = int(l[3])

        yield chrom,pos,l,line

def generate_line_prune(f):

    for line in f:
        l = line.strip().split(':')
        chrom = l[0]
        pos = int(l[1])
        yield chrom,pos,l,line


def parse_args():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        '--files1',
        dest='files1',nargs='+',
        metavar='FILE',default=None,
        required = True,
        )

    parser.add_argument(
        '--files2',
        dest='files2',nargs='+',
        metavar='FILE',default=None,
        required = True,
        )

    parser.add_argument(
        '--affix3',
        dest='affix3',
        metavar='FILE',default=None,
        required = True,
        )

    parser.add_argument(
        '--format1',
        dest='format1',
        choices=['gprobs'],
        metavar='STRING',default=None,
        required = True,
        )

    parser.add_argument(
        '--format2',
        dest='format2',
        choices=['prune.in'],
        metavar='STRING',default=None,
        required = True,
        )

    parser.add_argument(
        '--format3',
        dest='format3',
        choices=['EIGENSOFT'],
        metavar='STRING',default=None,
        required = True,
        )

##    ## http://docs.python.org/2/library/functions.html#vars
##    d_args = {}
##    print(parser.parse_args())
##    print(dict(parser.parse_args()))
##    for k,v in vars(parser.parse_args()).items():
##        d_args[k] = v
##        continue

    d = dict(vars(parser.parse_args()).items())

    return d


if __name__ == '__main__':

    main()
