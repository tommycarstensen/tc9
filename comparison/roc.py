#!/bin/python3

## Tommy Carstensen, Wellcome Trust Sanger Institute
## February-March 2013, October-November 2014

import argparse
import fileinput
import os
import gzip
import contextlib
import bz2
import itertools
import operator
import heapq

## Global dictionary for converting chromosomes to sortable integers.
d_chroms = {str(i):i-1 for i in range(1, 23)}
d_chroms['X'] = 22
d_chroms['Y'] = 23
d_chroms['MT'] = 24

def main():

    args = argparser()

    if not args.input_sorted:
        args.calls = sort_nicely(args.calls)

    ## Put file objects in a contextlib to ensure safe file object closure.
    with \
        fileinput.FileInput(
            files=args.calls, openhook=hook_compressed_text) as file_calls, \
        open_file(args.baseline, 'r') as file_vcf_baseline, \
        open_file(args.bed) as file_bed, \
        open_file(args.recal, file_type = 'recal') as file_recal, \
        open_file(args.out, 'w') as file_out, \
        open_file(args.out+'.tp.gz', 'w') as file_out_tp, \
        open_file(args.out+'.fp.gz', 'w') as file_out_fp, \
        open_file(args.out+'.fn.gz', 'w') as file_out_fn, \
        open_file(args.ref) as file_ref:

        initiate_loop_lines(
            args,
            file_calls, file_vcf_baseline, file_out, file_bed,
            file_recal, file_ref,
            file_out_tp, file_out_fp, file_out_fn,
            )

    return


def key_generator_recal(file):

    ## This function necessary prior to Python 3.5.

    for line in file:
        if line[0] == '#':
            continue
        chrom, pos, _ = line.split('\t',2)
        yield (d_chroms[chrom], int(pos)), line

    return

def heapq_merge(l_it):

    for key, line in heapq.merge(*l_it):
        yield line

@contextlib.contextmanager
def open_file(path, mode='r', file_type=None):

    if not path:
        yield None
    elif path == '-':
        assert mode in ('r', 'w')
        if mode == 'w':
            return sys.stdout
        elif mode == 'r':
            return sys.stdin
    elif type(path) == list:
        if file_type == 'recal':
            l_it = []
            for p in path:
                with open_file(p) as file:
                    it = key_generator_recal(file)
                    l_it.append(it)
            yield heapq_merge(l_it)
        else:
            print(path)
            stop
    elif type(path) == str:
        if os.path.splitext(path)[-1] == '.gz':
            yield gzip.open(path, mode + 't')
        else:
            yield open(path, mode)
    else:
        stop

    return

def parse_split_line(l, col, file_ref, d_fai):

    GT = l[col].split(':')[0].replace('|', '/')
    chrom = l[0]
    pos = int(l[1])
    ref = l[3]
    alt = l[4]

    ref_norm, gt_str, pos_norm = GT_int2str(
        file_ref, d_fai, chrom, pos, ref, alt, GT)

    d = {
        'chrom_str':chrom, 'pos':pos, 'l':l, 'GT':GT,
        'pos_norm':pos_norm, 'gt_str':gt_str, 'ref':ref, 'alt':alt,
        'chrom_int':d_chroms[chrom],
        }

    return d


def initiate_loop_lines(
    args,
    file_calls, file_baseline, file_out, file_bed,
    file_recal, file_ref,
    file_out_tp, file_out_fp, file_out_fn,
    ):

    col1 = col_calls = sampleID2col(args.sampleID_calls, file_calls)
    col2 = col_baseline = sampleID2col(args.sampleID_baseline, file_baseline)

    d_fai = read_fai(args.ref + '.fai')

    line2, l2, chrom2, pos2 = next(split_line_vcf(file_baseline))
    d_baseline = parse_split_line(l2, col2, file_ref, d_fai)

    if file_bed:
        chrom_bed, pos1_bed, pos2_bed = next(split_line_bed(file_bed))

    nFalse = 0
    nTrue = 0

    if args.type == 'INDEL':
        func_check_type = is_indel
    elif args.type == 'SNP':
        func_check_type = is_snp
    else:
        print(args.type)
        stop
        func_check_type = is_any

    if args.sort == 'VQSLOD':
        if file_recal:
            func_return_score = parse_VQSLOD_from_recal
        else:
            func_return_score = parse_VQSLOD_from_self
    elif args.sort == 'QUAL':
        func_return_score = operator.itemgetter(5)
    else:
        stop

    d_bed = {'chrom_str':chrom_bed, 'pos1':pos1_bed, 'pos2':pos2_bed}

    ## Also print these to file.
    d_type_errors = {'fp':[], 'tp':[], 'fn':0}

    for line1, l1, chrom1, pos1 in split_line_vcf(file_calls):

        ref1 = l1[3]
        alt1 = l1[4]

        process_line(
            line1, l1, chrom1, pos1, col1, ref1, alt1,
            col2,
            args,
            file_recal,
            file_bed, d_bed,
            file_ref, d_fai,
            file_baseline, d_baseline,
            file_out_tp, file_out_fp, file_out_fn,
            d_type_errors,
            func_check_type, func_return_score)
##    while True:
##        try:
##            l1, chrom1, pos1 = next(split_line_vcf(file_calls))
##        except StopIteration:
##            break

    print('False', nFalse)
    print('True', nTrue)
    print('tp', tp)  # baseline & calls
    print('fp', tp)  # calls - baseline
    print('fn', fn)  # baseline - calls
    compare_identity_to_vcfeval
    write_tp_fp_tpr_fdr_to_file_and_associated_pos
    write_fdr_to_file_by_counting_bed_range_or_maybe_not

    return

def process_line(
    line1, l1, chrom1, pos1, col1, ref1, alt1,
    col2,
    args,
    file_recal,
    file_bed, d_bed,
    file_ref, d_fai,
    file_baseline, d_baseline,
    file_out_tp, file_out_fp, file_out_fn,
    d_type_errors,
    func_check_type, func_return_score):

    ## skip SNP/indel
    if len(ref1) > 1 and all(
        [len(alt_i) == len(ref1) for alt_i in alt1.split(',')]):
        print('Please convert MNPs to SNPs before using this script.')
        print('MNP', chrom1, pos1, ref1, alt1)
        sys.exit()

    ## Parse GT.
    GT1 = l1[col1].split(':')[0]

    ## Keep recal file in sync if provided.
    score = func_return_score(
        file_recal=file_recal, chrom=chrom1, pos=pos1, l=l1)

    ## skip non-call and REFREF-call
    if GT1 == './.' or GT1 == '0/0':
        return

    ## Skip if not correct variant type.
    bool_type = func_check_type(ref1, alt1)
    if not bool_type:
        return

    ## Skip positions not in bed file.
    if file_bed:
        while d_chroms[chrom1] > d_chroms[d_bed['chrom_str']]:
            d_bed['chrom_str'], d_bed['pos1'], d_bed['pos2'] = next(
                split_line_bed(file_bed))
        while pos1 > d_bed['pos2']:
            d_bed['chrom_str'], d_bed['pos1'], d_bed['pos2'] = next(
                split_line_bed(file_bed))
        if d_bed['pos1'] > pos1:
            return
        pass

    QUAL = l1[5]

    ## Left align and trim. And convert GT_int to GT_str.
    print(pos1, ref1, alt1, GT1, col1, len(l1))
    ref1_norm, gt1_str, pos1_norm = GT_int2str(
        file_ref, d_fai, chrom1, pos1, ref1, alt1, GT1)

    if pos1_norm != pos1:
        print(chrom1, pos1, ref1, alt1, GT1)
        print(ref1_norm, gt1_str, pos1_norm)
        stop

    bool_TF = None
    bool_break = False
    bool_fn = False
    while True:
        if chrom1 != d_baseline['chrom_str']:
            if d_chroms[chrom1] > d_chroms[d_baseline['chrom_str']]:
                bool_fn = True
                pass
            else:
                bool_TF = False
                break
            pass
        elif pos1 == d_baseline['pos']:
            if bool_TF != None and l2[6] != 'PASS':
                stop1
                pass
            tmp = [ref1] + alt1.split(',')
            if d_baseline['pos'] != d_baseline['pos_norm']:
                stoptmp3
            ## Identical genotypes.
            if gt1_str == d_baseline['gt_str']:
                bool_TF = True
            ## Identical genotypes, but phasing different or missing.
            elif gt1_str == list(reversed(d_baseline['gt_str'])):
                bool_TF = True
            elif [tmp[int(i)] for i in GT1.split('/')] == d_baseline['gt_str']:
                stoptmp8
            ## different genotypes
            else:
                bool_TF = False
                ## identical REF and ALT (bi- and multiallelic)
                if ref1 == d_baseline['ref'] and len(
                    set(alt1.split(',')) & set(d_baseline['alt'].split(','))) > 0:
                    pass
                ## different REF and ALT
                else:
                    print(chrom1, pos1, 'GT1', GT1, 'GT2', d_baseline['GT'])
                    print('ref1', ref1, 'alt1', alt1)
                    print('ref2', d_baseline['ref'], 'alt2', d_baseline['alt'])
                    print(l1[:5], l1[col1])
                    print(d_baseline['l'][:5], d_baseline['l'][col2])
                    import time
                    time.sleep(.5)
            break
##                ## compare genotypes
####                bool_TF = compareGT(GT1,ref1,alt1,GT2,ref2,alt2,bool_TF)
##                ## do not compare genotypes/haplotypes,
##                ## but only positions/sites
##                bool_TF = True
##                pass
        elif pos1_norm == d_baseline['pos']:
            stoptmp4a
        elif pos1_norm == d_baseline['pos_norm']:
            stoptmp4b
        elif pos1 == d_baseline['pos']:
            stoptmp4c
        elif pos1 > d_baseline['pos']:
            if (
                pos1 - d_baseline['pos'] < 5 and
                len(ref1) > 1 and len(alt1) > 1):
                print(pos1, d_baseline['pos'])
                print(l1[:9], l1[-1])
                print(l2)
                stoptmp1
            bool_fn = True
            pass
        ## elif pos2 > pos1:
        else:
            if bool_TF == None:
                bool_TF = False  # FP or TN
            else:
                stop2
            break
        if bool_fn:
            d_type_errors['fn'] += 1
            print('fn', d_type_errors['fn'], d_baseline['chrom_str'], d_baseline['pos'])
            try:
                while True:
                    (
                        line2, d_baseline['l'], d_baseline['chrom_str'],
                        d_baseline['pos']) = next(split_line_vcf(
                            file_baseline))
                    for k, v in parse_split_line(
                        d_baseline['l'], col2, file_ref, d_fai).items():
                        d_baseline[k] = v
                    if d_baseline['ref'] == d_baseline['alt']:
                        print(d_baseline)
                        stoptmp
                        continue
                    print(
                        d_baseline['pos'], d_baseline['pos_norm'],
                        d_baseline['GT'], d_baseline['gt_str'])
##                    assert d_baseline['pos'] == d_baseline['pos_norm']  # tmp!!!
                    break
            except StopIteration:
                bool_break = True
            bool_fn = None
            pass
        if bool_break:
            break
        continue

    ## Last variant of baseline.
    if bool_break:
        stop_handle_this

    assert bool_TF != None
    ## Assert that we didn't read ahead in the baseline file.
    assert pos1 <= d_baseline['pos']

    QUAL = l1[5]

    if len(gt1_str[0]) == len(ref1) and len(gt1_str[1]) == len(ref1):
        bool_indel = False
        type_variant = 'SNP'
        if not len(ref1) == 1 and not len(_ref1) == 1:  # assert not MNP
            print('ref alt GT', ref1, alt1, GT1)
            print(gt1, _ref1)
            print('trim1', trim_right(ref1, gt1_str[0]))
            print('trim2', trim_right(ref1, gt1_str[1]))
            print(GT_int2str(GT1, ref1, alt1, pos1, file_ref))
            print(ref2, alt2)
            print(gt2)
            stop
    else:
        bool_indel = True
        type_variant = 'indel'

##    ## todo: write vcf file with tp/fp annotations
##    line = '{chrom}\t{pos}\t{bool_TF}\t'.format(
##        chrom=chrom1, pos=pos1, bool_TF=bool_TF)
##    line += '{}\t{}\t{}\t'.format(QUAL, score, GT1)
##    line += '{}\t{}\t'.format(ref1, alt1)
##    line += '\n'
##    file_out.write(line)
    print(pos1, bool_TF)
    if bool_TF == False:
        d_type_errors['fp'].append(score)
        file_out_fp.write(line1)
    else:
        d_type_errors['tp'].append(score)
        file_out_tp.write(line1)

    return


def sampleID2col(sampleID, file):

    for line in file:
        if line[:2] == '##':
            continue
        break
    if sampleID:
        col = line.rstrip().split('\t').index(sampleID)
    else:
        assert len(line.rstrip().split('\t')) == 10
        col = 9

    return col


def split_line_bed(file_bed):

    for line in file_bed:
        chrom, pos1, pos2 = file_bed.readline().split('\t')
        pos1 = int(pos1)
        pos2 = int(pos2)
        yield chrom, pos1, pos2

    return


def read_fai(path_fai):

    if not os.path.exists(path_fai):
        sys.stderr.write('Does not exist: {}'.format(path_fai))
        sys.exit()

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


def is_any():

    return True


def GT_int2str(file_ref, d_fai, chrom, pos, ref, alt, GT):

    l = [ref]+alt.split(',')
    gt = [l[int(i)] for i in GT.split('/')]
##    print('len(ref)', len(ref))
    ## Left align if variant.
    if gt[0] != gt[1]:
        ref, gt[0], gt[1], pos = trim_and_left_align(
            file_ref, d_fai, chrom, pos, ref, gt[0], gt[1])
##    if len(ref) > 1 and len(gt[0]) > 1 and len(gt[1]) > 1:
####        index = min(len(ref), len(gt[0]), len(gt[1]))
####        ## e.g. ug2g 1:98999 TTTTA TTTTATTTA,T 1/1
####        if ref[:index] == gt[0][:index] and ref[:index] == gt[1][:index]:
######            print(ref, gt)
######            print(ref[index-1:], gt[0][index-1:], gt[1][index-1:])
######            print(pos, pos+index-1)
######            print(GT)
####            print('xxxxxxxxxxxxxxxxxxxxxxxxxx')
####            return ref[index-1:], [gt[0][index-1:], gt[1][index-1:]], pos+index-1
######        ## e.g. ug2g 1:691544 TA AA,T 1/1
####        if len(ref) == len(gt[0]) == len(gt[1]):
####            print(pos, ref, alt, GT, gt)
####            import time
####            time.sleep(1)
####            for i in range(2):
####                gt[i] = trim_right(ref, gt[i])[1]
####        else:
##        if True:
####                ref, gt[0], gt[1], pos = trim3(ref, gt[0], gt[1], pos)
##            ref, gt[0], gt[1], pos = trim_and_left_align(
##                file_ref, d_fai, chrom, pos, ref, gt[0], gt[1])
####            return ref, gt, pos
####            print(pos, ref, alt, GT, gt)
####            print(pos, GT)
####            print(ref)
####            print(alt)
####            print(gt)
####            stoptmp7
####            for i in range(2):
####    ##            if ref[:index] == gt[i][:index]:
####    ##            if len(ref) == len(gt[i]) and ref[-1] == gt[i][-1]:  # removed Nov2
####    ##            print(id(gt_trimmed), id(gt[i]), gt[i])
####                try:
####                    gt[i] = trim_right(ref, gt[i])[1]
####                except:
####                    print(i, pos, GT)
####                    print(ref)
####                    print(alt)
####                    print(gt)
####                    stop
####    ##            print(id(gt_trimmed), id(gt[i]), id(ref))
    assert len(gt[0]) > 0 and len(gt[1]) > 0

    return ref, gt, pos


def trim_and_left_align(file_ref, d_fai, chrom, pos, ref, gt1, gt2):

    assert len(ref) > 0 and len(gt1) > 0 and len(gt2) > 0
    ## right trim and left align
    while ref[-1] == gt1[-1] and ref[-1] == gt2[-1]:
        ## right trim
        ref = ref[:-1]
        gt1 = gt1[:-1]
        gt2 = gt2[:-1]
        ## left align
        if len(ref) == 0 or len(gt1) == 0 or len(gt2) == 0:
            pos -= 1
            nt = parse_ref(d_fai, file_ref, chrom, pos)
            ref = nt+ref
            gt1 = nt+ref
            gt2 = nt+ref

    ## left trim
    while ref[0] == gt1[0] and ref[0] == gt2[0]:
        if len(ref) == 1 or len(gt1) == 1 or len(gt2) == 1:
            break
        pos += 1
        ref = ref[1:]
        gt1 = gt1[1:]
        gt2 = gt2[1:]

    return ref, gt1, gt2, pos


def parse_ref(d_fai, file_ref, CHROM, POS, size=1):

    ## parse ref
    cnt = cnt_chars_per_line_excluding_newline = 60
    row1 = (POS-1)//cnt
    row2 = (POS-1+size)//cnt
    size += row2-row1
    col = (POS-1)%cnt
    byte_init = d_fai[CHROM]['start']
    offset = byte_init+(cnt+1)*row1+col
    file_ref.seek(offset)
    read = file_ref.read(size).replace('\n','')

    return read


def trim3(s1, s2, s3, pos):

##    len_min = min(len(s1), len(s2), len(s3))

    assert s1 != s2 or s1 != s3
##    else: return s1,s2,s3

##    print()
##    print(s1)
##    print(s2)
##    print(s3)

    ## right trim
    while s1[-1] == s2[-1] and s1[-1] == s3[-1]:
        s1 = s1[:-1]
        s2 = s2[:-1]
        s3 = s3[:-1]
##    for j, (a, b, c) in enumerate(zip(reversed(s1), reversed(s2), reversed(s3))):
##        if a != b or a != c:
##            break
##    print()
##    print('j', j)
##    print(min(len(s1), len(s2), len(s3)))
##    print(s1[:-j])
##    print(s2[:-j])
##    print(s3[:-j])
##    if j > 0:
##        s1 = s1[:-j]
##        s2 = s2[:-j]
##        s3 = s3[:-j]

    assert len(s1) > 0 and len(s2) > 0 and len(s3) > 0

##    ## left trim unnecessary
##    if len_min-j == 1:
##        ## assert left alignment not necessary
##        assert s1[0] != s2[0] or s1[0] != s3[0]
##        return s1, s2, s3, pos

    ## left trim
    while s1[0] == s2[0] and s1[0] == s3[0]:
        s1 = s1[1:]
        s2 = s2[1:]
        s3 = s3[1:]
        pos += 1
##    for i, (a, b, c) in enumerate(zip(s1, s2, s3)):
##        if a != b or a != c:
##            break
##    print()
##    print('i', i)
##    print(min(len(s1), len(s2), len(s3)))
##    print(s1[i:])
##    print(s2[i:])
##    print(s3[i:])
##    s1 = s1[i:]
##    s2 = s2[i:]
##    s3 = s3[i:]
    assert len(s1) > 0 and len(s2) > 0 and len(s3) > 0

##    pos += i

##    print(s1, s2, s3, pos)
##    stop

    return s1, s2, s3, pos


def trim_left(s1, s2):

    return


def trim_right2(s1, s2):

    for i, (a, b) in enumerate(zip(s1, s2)):
        if a != b:
            break
    print('s1', s1, 's2', s2, 'i', i)
    print(s1[:-i], s2[:-i])
    stop
    if i == 0:
        return s1, s2
    else:
        return s1[:-i], s2[:-i]


def trim_right(s1, s2):

    assert len(s1) == len(s2)  # removed Nov2
    assert len(s1) >= len(s2)
    for i, (a, b) in enumerate(zip(reversed(s1), reversed(s2))):
        if a != b:
            break
    print('s1', s1, 's2', s2, 'i', i)
    if i == 0:
        return s1, s2
    else:
        return s1[:-i], s2[:-i]


def parse_VQSLOD_from_recal(file_recal=None, chrom=None, pos=None, **kwargs):

    line, l_recal, chrom_recal, pos_recal = next(split_line_vcf(file_recal))
    if chrom_recal != chrom:
        stop1
    if pos_recal != pos:
        print(chrom_recal, pos_recal, pos)
        stop2
    for s in l_recal[7].split(';'):
        l = s.split('=')
        if l[0] == 'VQSLOD':
            VQSLOD = l[1]
            break

    return VQSLOD


def parse_VQSLOD_from_self(l=None, **kwargs):

    for s in l[7].split(';'):
        l = s.split('=')
        if l[0] == 'VQSLOD':
            VQSLOD = l[1]
            break

    return VQSLOD


def compareGT(GT1,ref1,alt1,GT2,ref2,alt2,bool_TF):

    d = {0:[], 1:[]}
    for i, (GT, alt, ref) in enumerate(
        [(GT1,alt1,ref1),(GT2,alt2,ref2)]):
##                    print(pos1,pos2,i,GT,alt,ref)
        l = [ref]+alt.split(',')
        try:
            d[i] = [l[int(j)] for j in GT.split('/')]
        except ValueError:
            print(s)
            print(l2)
            stop5
    if d[0] == d[1] or d[0] == list(reversed(d[1])):
        if bool_TF == False:
            print('\n')
            print(pos1,pos2)
            print(GT1,ref1,alt1)
            print(GT2,ref2,alt2)
            print(l1[col1])
            print(l2[col2])
            stop4
        bool_TF = True  # TP and FN
    else:
        if set(d[0]) != set(d[1]):
            if bool_TF == True:
                print('\n')
                print(pos1,pos2)
                print(GT1,ref1,alt1)
                print(GT2,ref2,alt2)
                print(l1[col1])
                print(l2[col2])
                stop5
            bool_TF = False
        else:
            print('\n')
            print(d)
            print(pos1,ref1,alt1,GT)
            print(l2[:9])
            print(l2[col2])
            print(l1[col1])
            print(list(reversed(d[1])))
            stop1
##                break
    return bool_TF


def is_snp(ref, alt):

    return not is_indel(ref, alt)


def is_indel(REF, ALT):

    ## 1) biallelic SNP
    if len(REF) == 1 and len(ALT) == 1:
        bool_indel = False
    else:
        l_ALT = ALT.split(',')
        if all([len(ALT_i) == len(REF) for ALT_i in l_ALT]):
            ## 2) multiallelic SNP
            if len(REF) == 1:
                bool_indel = False
            ## 3) mnp
            else:
                bool_indel = False
                l = [
                    [x == y for x, y in zip(REF, l_ALT[i])].count(False) == 1
                    for i in range(len(l_ALT))]
                if all(l):
                    type_variant = 'SNP'
                elif any(l):
                    type_variant = 'SNP,MNP'
                else:
                    type_variant = 'MNP'
                stop_write_more_code
        ## skip deletions
        elif not REF in ['A', 'C', 'G', 'T',]:
            if REF in iter(','.join(tup) for tup in itertools.permutations('ACGT',2)):
                stoptmp
            bool_indel = True
        ## skip insertions
        else:
            bool_indel = True

    return bool_indel


def split_line_bgl(file_bgl):

    for line in file_bgl:
        if line[:6] == 'marker':
            continue
        l = line.rstrip().split()
        chrom, pos = l[0].split(':')
        pos = int(pos)
##        ref = l[3]
##        alt = l[4]
        yield l, chrom, pos

    return


def split_line_vcf(file_vcf):

    for line in file_vcf:
        if line == '':
            return
        if line[0] == '#':
            continue
        l = line.rstrip().split('\t')
        chrom = l[0]
        pos = int(l[1])
##        ref = l[3]
##        alt = l[4]
        yield line, l, chrom, pos

    return


def hook_compressed_text(filename, mode):

    ##http://stackoverflow.com/questions/21529163/python-gzipped-fileinput-returns-binary-string-instead-of-text-string/21529243

    ext = os.path.splitext(filename)[1]
    if ext == '.gz':
        f = gzip.open(filename, mode+'t')
    elif ext == '.bz2':
        f = bz2.open(filename, mode+'t')
    else:
        f = open(filename, mode)

    return f


def alphanum_key(s):
    import re
    ## http://dave.st.germa.in/blog/2007/12/11/exception-handling-slow/
    NUM_RE = re.compile('([0-9]+)')
    return [int(c) if c.isdigit() else c for c in NUM_RE.split(s)]


def sort_nicely(l):
    ## http://nedbatchelder.com/blog/200712/human_sorting.html
    """ Sort the given list in the way that humans expect.
    """
    l.sort(key=alphanum_key)
    return l


def argparser():

    parser = argparse.ArgumentParser()

    ## Baseline; e.g. GiaB.
    parser.add_argument(
        '--baseline', '--gold', '--truth', '--path_b', required=True)

    ## Variant calls.
    parser.add_argument(
        '--calls', '--predict', '--path_c', required=True, nargs='+',
        default='-')

    parser.add_argument('--input_sorted', action='store_true', default=False)

    ## For left alignment of indels.
    parser.add_argument('--ref', '--path_r', required=True)

    ## For VQSLOD values.
    parser.add_argument('--recal', '--path_recal', required=False, nargs='*')

    ## For annotating ROC curve.
    parser.add_argument('--tranches', '--path_tranches', required=False)

    parser.add_argument('--sampleID_calls', required=False)

    parser.add_argument(
        '--sampleID_baseline', '--sampleID_truth', required=False)

    parser.add_argument(
        '--sort', '--vcf-score-field', '-f', '--sortval', default='QUAL',
        choices=['QUAL', 'VQSLOD'])

    parser.add_argument('--out', '--o', '-o', default='-')

    ## For filtering of calls and baseline.
    parser.add_argument('--bed', required=False)

    ## For filtering of variant types.
    parser.add_argument(
        '--type', '--mode', default='SNP', choices=['SNP', 'INDEL'])

    args = parser.parse_args()

    return args


if __name__ == '__main__':
    main()
