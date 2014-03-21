#!/bin/python3

## Tommy Carstensen, Wellcome Trust Sanger Institute, March 2014

## http://pngu.mgh.harvard.edu/~purcell/plink/binary.shtml

import argparse
import os
##import re
##import datetime
import sys
import contextlib
import fileinput
import gzip
import bz2
import math


def main():

    d_args = argparser()

##    l_keep_sampleIDs = parse_keep(fp_keep,l_fam_sampleIDs)
##
##    l_keep_index = index_keep(l_keep_sampleIDs,l_fam_sampleIDs)
##
##    l_sampleIDs_vcf = shorten_sampleIDs(l_keep_sampleIDs,fp_update_ids)

##    with contextlib.ExitStack() as stack:
##        n_samples1 =
##        n_samples2 = 

    ## in the future don't use contextlib and fileinput,
    ## but rather loop over files and then lines within a generator
    with contextlib.ExitStack() as stack:
        within_stack(stack, d_args)

    return


def within_stack(stack, d_args):

    ext_out = os.path.splitext(d_args['out'])[1]

    d_kwargs, d_ext = open_files(stack, d_args)
    d_out, print_out = write_header(stack,d_args)

    ext1 = d_ext[1]
    ext2 = d_ext[2]

    d_parse_pos = {
        '.vcf':parse_pos_vcf, '.bed':parse_pos_bed, '.gen':parse_pos_gen}
    parse_pos1 = d_parse_pos[ext1](ext_out, **d_kwargs[1])
    parse_pos2 = d_parse_pos[ext2](ext_out, **d_kwargs[2])

    d_parse_samples = {
        '.vcf':parse_samples_vcf, '.bed':parse_samples_bed,
        '.gen':parse_samples_gen}
    l_samples1 = d_parse_samples[ext1](d_args['in1'][0])
    n_samples1 = len(l_samples1)
    l_samples2 = d_parse_samples[ext2](d_args['in2'][0])
    n_samples2 = len(l_samples2)

    if ext_out == '.bed':
        ## lose sex information... should probably change that...
        open(d_args['out'][:-4]+'.fam','w').writelines(
            ['%s %s 0 0 -9 -9\n' %(sampleID, sampleID) for sampleID in l_samples1+l_samples2])
    else:
        stop_write_more_code

    d_parse_gt = {
        '.vcf':parse_gt_vcf, '.bed':parse_gt_bed, '.gen':parse_gt_gen}
    parse_gt1 = d_parse_gt[ext1]
    parse_gt2 = d_parse_gt[ext2]

    chrom1, pos1, rsID1, alleleA1, alleleB1, gt1 = next(parse_pos1)
    chrom2, pos2, rsID2, alleleA2, alleleB2, gt2 = next(parse_pos2)

    kwargs_gt1 = get_kwargs_gt(ext1, n_samples1)
    kwargs_gt2 = get_kwargs_gt(ext2, n_samples2)

    while True:

        if chrom1 != chrom2:
            index1 = d_args['chromosomes'].index(chrom1)
            index2 = d_args['chromosomes'].index(chrom2)
            if index1 > index2:
                chrom2, pos2, rsID2, alleleA2, alleleB2, gt2 = next(parse_pos2)
            else:
                chrom1, pos1, rsID1, alleleA1, alleleB1, gt1 = next(parse_pos2)
        elif pos1 == pos2:
            if alleleA1 == alleleA2 and alleleB1 == alleleB2:
                bool_rev = False
            elif alleleA1 == alleleB2 and alleleB1 == alleleA2:
                bool_rev = True
            elif ext2 == '.bed' and alleleB1 == alleleB2 and alleleA2 == '0':
                bool_rev = False
            elif ext2 == '.bed' and alleleA1 == alleleB2 and alleleA2 == '0':
                bool_rev = True
            else:
                print('warning: allele mismatch')
                print(chrom1, pos1, alleleA1, alleleB1)
                print(chrom2, pos2, alleleA2, alleleB2)
                chrom1, pos1, rsID1, alleleA1, alleleB1, gt1 = next(parse_pos1)
                chrom2, pos2, rsID2, alleleA2, alleleB2, gt2 = next(parse_pos2)
                continue
            gt1 = parse_gt1(gt1, ext_out, **kwargs_gt1)
            gt2 = parse_gt2(gt2, ext_out, bool_rev = bool_rev, **kwargs_gt2)
            print_out(
                ext_out, d_out, gt1, gt2,
                chrom1, pos1, rsID1, alleleA1, alleleB1)
            chrom1, pos1, rsID1, alleleA1, alleleB1, gt1 = next(parse_pos1)
            chrom2, pos2, rsID2, alleleA2, alleleB2, gt2 = next(parse_pos2)
        elif pos1 > pos2:
            chrom2, pos2, rsID2, alleleA2, alleleB2, gt2 = next(parse_pos2)
##            elif pos2 > pos1:
        else:
            try:
                chrom1, pos1, rsID1, alleleA1, alleleB1, gt1 = next(parse_pos1)
            except StopIteration:
                break

    return


def write_header(stack,d_args):

    ext_out = os.path.splitext(d_args['out'])[1]

    if ext_out == '.bed':
        fd_out = stack.enter_context(open(d_args['out'],'wb'))
        magic_number = bytearray([108,27])
        mode = bytearray([1])
        fd_out.write(magic_number+mode)
        bim_out = stack.enter_context(open(d_args['out'][:-4]+'.bim','w'))
        d_out = {ext_out:fd_out,'.bim':bim_out}
        print_out = print_out_bed
    else:
        fd_out = stack.enter_context(open(d_args['out'],'w'))
        d_out = {ext_out:fd_out}
        stop_write_more_code

    return d_out, print_out


def print_out_vcf(
    ext_out, d_out, gt1, gt2, chrom, pos, rsID, alleleA, alleleB):

    stop_write_more_code

    return


def print_out_gen(
    ext_out, d_out, gt1, gt2, chrom, pos, rsID, alleleA, alleleB):

    stop_write_more_code

    return


def print_out_bed(
    ext_out, d_out, gt1, gt2, chrom, pos, rsID, alleleA, alleleB):

    fd_out = d_out[ext_out]

    bits_per_byte = 8
    gt = gt1+gt2+'0'*int(len(gt1+gt2)%bits_per_byte)
    barray = bytearray()
    for i in range(0,len(gt),bits_per_byte):
        barray.append(int(gt[i:i+8][::-1],2))
    fd_out.write(bytes(barray))

    line_bim = '%s\t%s\t0\t%s\t%s\t%s\n' %(chrom,rsID,pos,alleleA,alleleB)
    d_out['.bim'].write(line_bim)

    return


def get_kwargs_gt(ext, n_samples):

    kwargs = {
        '.bed':{'n_samples':n_samples},
        '.vcf':{},
        }[ext]

    return kwargs


def parse_gt_bed(gt, ext_out, n_samples, bool_rev = False):

    if ext_out != '.bed':
        stop_write_more_code

    if bool_rev:
        d_gt = {
            '.vcf':{'00':'1/1','10':'./.','11':'0/0','01':'0/1',},
            '.bed':{'00':'11','10':'10','11':'00','01':'01',},
            }
    else:
        d_gt = {
            '.vcf':{'00':'0/0','10':'./.','11':'1/1','01':'0/1',},
            '.bed':{'00':'00','10':'10','11':'11','01':'01',},
            }

    d_sep = {'.bed':''}

    gt_out = ''
    cnt_samples = 0
    for byte in gt:
##                        s8 = str(bin(byte)[2:]).zfill(8)
        s8 = '{:08b}'.format(byte)#[::-1]
        for i in range(8-1,1-1,-2):
            cnt_samples += 1
            if cnt_samples > n_samples:
                break
            GT = s8[i-1:i+1][::-1]
            gt_out += d_sep[ext_out]+d_gt[ext_out][GT]

    return gt_out


def parse_pos_bed(ext_out, fileinput=None, bim=None, fam=None):

    files = fileinput

    n_samples = len(fam.readlines())
    n_bytes_per_SNP = math.ceil(n_samples/4)

##    n_bytes_per_SNP = math.ceil(n_samples/4)
##
##    d_gt = {
##        '.vcf':{'00':'0/0','10':'./.','11':'1/1','01':'0/1',},
##        '.bed':{'00':'00','10':'10','11':'11','01':'01',},
##        }
##    d_sep = {'.vcf':'\t','.bed':''}
##    ext_out = '.vcf' ## tmp!!!

    for file in files:
        with open(file,'rb') as bed:
            magic_number = bytearray([108,27])
            mode = bytearray([1])
            header = bed.read(len(magic_number)+len(mode))
            if header != magic_number+mode:
                print(header)
                print(magic_number, mode)
                stop
            while True:
                l = bim.readline().rstrip().split()
                chrom = l[0]
                rsID = l[1]
                pos = int(l[3])
                alleleA = l[4]
                alleleB = l[5]
                bytesSNP = bed.read(n_bytes_per_SNP)
                ## skip mirror alleles
                if alleleA == 'C' and alleleB == 'G':
                    continue
                if alleleA == 'A' and alleleB == 'T':
                    continue
                if alleleA == 'G' and alleleB == 'C':
                    continue
                if alleleA == 'T' and alleleB == 'A':
                    continue
                ## skip monomorphic alleles
                if alleleA == '0' or alleleB == '0':
                    continue

##                if False: ## slow
##                    genotype_fields = ''
##                    for i_sample in range(n_samples):
##                        i_byte = i_sample//4
##                        byte = bytesSNP[i_byte]
##        ##                s8 = str(bin(byte)[2:]).zfill(8)
##                        s8 = '{:08b}'.format(byte)#[::-1]
##                        i = i_sample%4
##                        s2 = s8[8-2*i-2:8-2*i][::-1]
##                        if s2 == '10':
##                            genotype_fields += '\t./.'
##                        else:
##                            genotype_fields += '\t%s' %('/'.join(s2))
##                if True: ## fast
####                    genotype_fields2 = ''
##                    gt = ''
##                    cnt_samples = 0
##                    for i_byte in range(n_bytes_per_SNP):
##                        byte = bytesSNP[i_byte]
##    ##                        s8 = str(bin(byte)[2:]).zfill(8)
##                        s8 = '{:08b}'.format(byte)#[::-1]
##                        for i in range(8-1,1-1,-2):
##                            cnt_samples += 1
##                            if cnt_samples > n_samples:
##                                break
##                            s2 = s8[i-1:i+1][::-1]
##                            gt += s2
####                            genotype_fields2 += '\t%s' %(d_gt[ext_out][s2])
####                print('a',genotype_fields[:40])
####                print('b',genotype_fields2[:40])
####                print(genotype_fields2==genotype_fields)
####                stop

                yield chrom, pos, rsID, alleleA, alleleB, bytesSNP
            
    return


def junk():
    STOPSTOP
                
    with open(fp_bim,'r') as bim:
        n_SNPs = len(bim.readlines())

    with open(fp_bed,'rb') as bed, open(fp_vcf,'w') as vcf, open(fp_bim,'r') as bim:
        magic_number = bytearray([108,27])
        mode = bytearray([1])
        header = bed.read(len(magic_number)+len(mode))
        if header != magic_number+mode:
            print(header)
            print(magic_number, mode)
            stop
        vcf.write('##fileformat=VCFv4.0\n')
        vcf.write('##fileDate=%s\n' %(datetime.date.today()))
        vcf.write('##source=%s\n' %(fp_bed))
        vcf.write('##source=%s\n' %(sys.argv[0]))
        vcf.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT')
        for sampleID in l_sampleIDs_vcf:
            vcf.write('\t%s' %(sampleID))
        vcf.write('\n')
        i_SNP = 0
        while True:
            i_SNP += 1
            if i_SNP%1000 == 0:
                print('SNP %i of %i SNPs' %(i_SNP,n_SNPs))
            line_bim = bim.readline()
            if not line_bim: break
            l_bim = line_bim.rstrip().split()
            CHROM = l_bim[0]
            POS = l_bim[3]
            ID = l_bim[1]
            REF = l_bim[4]
            ALT = l_bim[5]
            ## monomorphic
            if ALT == '0':
                ALT = '.'
            QUAL = '.'
            FILTER = '.'
            FORMAT = 'GT'

            bytesSNP = bed.read(n_bytes_per_SNP)
            cnt_allele_ALT = 0
            cnt_alleles = 0
            genotype_fields = ''
            for i_keep in l_keep_index:
                i_byte = int(i_keep/4)
                byte = bytesSNP[i_byte]
##                s8 = str(bin(byte)[2:]).zfill(8)
                s8 = '{:08b}'.format(byte)#[::-1]
                i = i_keep%4
                s2 = s8[8-2*i-2:8-2*i]
                if s2 == '01':
                    genotype_fields += '\t./.'
                else:
                    cnt_alleles += 2
                    cnt_allele_ALT += s2.count('1')
                    genotype_fields += '\t%s' %('/'.join(s2))

            AF_ALT = cnt_allele_ALT/cnt_alleles
            INFO = 'AF=%f' %(AF_ALT)
            fixed_fields = '\t'.join((
                CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT))

            line_vcf = fixed_fields+genotype_fields+'\n'
            vcf.write(line_vcf)

    return


def parse_samples_vcf(vcf):

    with hook_compressed_text(vcf,'r') as f:
        for line in f:
            if line[0] != '#':
                break
            header = line
    l_samples = header.rstrip().split('\t')[9:]

    return l_samples


def parse_samples_gen():

    stop1

    return


def parse_samples_bed(bed):

    with open(bed[:-4]+'.fam') as f:
        l_samples = [line.split()[1] for line in f.readlines()]

    return l_samples


def parse_pos_gen():

    stop

    return


def parse_gt_vcf(gt, ext_out, bool_rev=False):

    if bool_rev:
        d_gt = {
            '.bed':{'./.':'10','0/0':'11','0/1':'01','1/1':'00',},#'1/0':'01',},
            }
    else:
        d_gt = {
            '.bed':{'./.':'10','0/0':'00','0/1':'01','1/1':'11',},#'1/0':'01',},
            }

    d_sep = {'.bed':''}

    gt_out = ''
    for s in gt.rstrip().split('\t'):
        GT = s.split(':')[0]
        gt_out += d_sep[ext_out]+d_gt[ext_out][GT]

    return gt_out


def parse_gt_gen(gt):

    stop

    return


def parse_pos_vcf(ext_out, fileinput=None):

    for line in fileinput:
        if line[0] == '#':
            continue
        l = line.split('\t',9)
        ## skip INDELs and triallelic SNPs
        ALT = l[4]
        if not ALT in ['A','C','G','T','N','.']:
            if len(ALT) == 1:
                print(ALT)
                stop
            continue
        REF = l[3]
        if not REF in ['A','C','G','T','N','.']:
            if len(REF) == 1:
                print(REF)
                stop
            continue
        CHROM = l[0]
        POS = int(l[1])
        ID = l[2]
        gt = l[9]

        yield CHROM, POS, ID, REF, ALT, gt

    return


def open_files(stack, d_args):

    ## rename this function to get_extension (or add --ext1 --ext2 arguments) - yes, do the latter! restrict to input choices immediately!

    d_kwargs = {}
    d_ext = {}
    for i in (1,2):
        d_kwargs[i] = {}
        files = d_args['in%i' %(i)]
        root,ext = os.path.splitext(files[0])
        if ext == '.gz':
            root,ext = os.path.splitext(root)
        if not ext in ['.bed','.vcf']:
            print(files[0])
            stop
        d_ext[i] = ext
        if ext == '.bed':
            if len(files) > 1:
                print(files)
                stop
            bim = stack.enter_context(fileinput.input(
                files=[file[:-4]+'.bim' for file in files]))
            fam = stack.enter_context(open(files[0][:-4]+'.fam'))
            for file in files:
                n_samples = len(open(file[:-4]+'.bim').readlines())
                n_SNPs = len(open(file[:-4]+'.fam').readlines())
                if os.path.getsize(file)-3 != n_samples*n_SNPs/4:
                    print(os.path.getsize(file))
                    print(n_samples)
                    print(n_SNPs)
                    print(n_samples*n_SNPs/4)
                    stop
            d_kwargs[i]['bim'] = bim
            d_kwargs[i]['fam'] = fam
##            d_kwargs[i]['fileinput'] = stack.enter_context(open(files[0], 'rb'))
            d_kwargs[i]['fileinput'] = files
        else:
            d_kwargs[i]['fileinput'] = stack.enter_context(fileinput.FileInput(
                files=files,openhook=hook_compressed_text))

    return d_kwargs, d_ext


def alphanum_key(s):
    import re
    ## http://dave.st.germa.in/blog/2007/12/11/exception-handling-slow/
    NUM_RE = re.compile('([0-9]+)')
    return [ int(c) if c.isdigit() else c for c in NUM_RE.split(s) ]


def sort_nicely(l):
    ## http://nedbatchelder.com/blog/200712/human_sorting.html
    """ Sort the given list in the way that humans expect.
    """
    l.sort(key=alphanum_key)
    return l


def hook_compressed_text(filename, mode):

    ##http://stackoverflow.com/questions/21529163/python-gzipped-fileinput-returns-binary-string-instead-of-text-string/21529243

    ext = os.path.splitext(filename)[1]
    if ext == '.gz':
        f = gzip.open(filename, mode+'t')
    elif ext == '.bz2':
        f = bz2.open(filename, mode + 't')
    else:
        f = open(filename, mode)

    return f


def line2barray(line,n_columns,chrom):

    l = line.split()

##    rsID = l[0]
    rsID = l[1]
    pos = l[2]
    alleleA = l[3]
    alleleB = l[4]
    line_bim = '%s\t%s\t%s\t%s\t%s\t%s\n' %(chrom,rsID,'0',pos,alleleA,alleleB)
    
    b = ''
    barray = bytearray()
    for i in range(5,n_columns,3):
        l_triple = l[i:i+3]
        if l_triple == 3*[l_triple[0]]:
            b += '10'
        else:
            pmax = max(l_triple)
            index = l_triple.index(pmax)
            if index == 0:
                b += '00'
            elif index == 1:
                b += '01'
            else:
                b += '11'
        if len(b) == 8:
##            print('bbb',b)
            barray.append(int(b[::-1],2))
            b = ''
##    if len(b) > 0:
    if b:
        barray.append(int(b[::-1].zfill(8),2))

    return barray, line_bim


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
    

def shorten_sampleIDs(l_sampleIDs_long,fp_update_ids):

    if fp_update_ids:
        with open(fp_update_ids) as f:
            d_update = {
                line.strip().split()[0]:line.strip().split()[2] for line in f}
    else:
        d_update = {sampleID:sampleID for sampleID in l_sampleIDs_long}

    keyword = re.compile(r'(\d\d\d\d\d\d_[A-H]\d\d_)(.+\d\d\d\d\d\d\d)')

    l_sampleIDs_short = []

    for sampleID_long in l_sampleIDs_long:
        sampleID_long = d_update[sampleID_long]
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


def isfile(str_):

    if not os.path.isfile(str_) and not os.path.islink(str_):
        msg = '%s is neither a readable file nor a symbolic link' % str_
        raise argparse.ArgumentTypeError(msg)

    return str_


def argparser():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        '--mode', choices=['intersection','union'], default = 'intersection')
    parser.add_argument('--out', required=True, default=sys.stdout)
    parser.add_argument('--in1', type=isfile, nargs='+')
    parser.add_argument('--in2', type=isfile, nargs='+')
    parser.add_argument(
        '--chromosomes', nargs='+',
        default=[str(i) for i in range(1,23)]+['X','Y'])

    d_args = vars(parser.parse_args())

    d_args['in1'] = sort_nicely(d_args['in1'])
    d_args['in2'] = sort_nicely(d_args['in2'])

    return d_args


if __name__ == '__main__':
    main()
