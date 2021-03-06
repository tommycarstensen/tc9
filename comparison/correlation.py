#!/bin/python3

## Tommy Carstensen, Wellcome Trust Sanger Institute, Oct2013, Apr2015

## THIS SCRIPT IS JUST PLAIN UGLY!!! FIX IT WHEN TIME!!!

## http://pngu.mgh.harvard.edu/~purcell/plink/binary.shtml

## todo: make it possible to read multiple bed files

## todo: give the user an option to submit a file that links chromosomes and
## IMPUTE2 output files --- files should then be sorted by chromosome

## gen,bed,bgl,vcf,hap vs gen,bed,bgl,vcf,hap
##

import argparse
import os
import math
import contextlib
import fileinput
import re
import sys
import gzip
import itertools
import pysam


def main():

    d_args = argparser()

    d_samples, d_indexes = index_samples(d_args)

    os.makedirs(os.path.dirname(d_args['affix']), exist_ok=True)

    with contextlib.ExitStack() as stack:

        d = open_files(stack, d_args)

        if d_args['extract']:
            f_extract = stack.enter_context(open(d_args['extract']))
        else:
            f_extract = None

        if d_args['discordant']:
            f_discordant = stack.enter_context(open(d_args['discordant'], 'w'))
        else:
            f_discordant = None

        loop_main(
            d_args['affix'], d_args['chrom'], f_extract, f_discordant,
            d_samples, d_indexes, d_args, **d)

    return


def open_files(stack, d_args):

    ## tidy up this function and delete duplicate code!!!

    ## http://docs.python.org/3/tutorial/controlflow.html
    d = {}

    for i in 1, 2:
        files = d_args['files%i' %(i)]
        print(files)
        d['format%i' %(i)] = d_args['format%i' %(i)]
#        markers = [markers1, markers2,][i-1]

        if d_args['format%i' %(i)] == 'bed':
            if len(files) != 1:
                print(files)
                stoptmp
            f_bed = stack.enter_context(open('%s.bed' %(files[0]), 'rb'))
            f_bim = stack.enter_context(open('%s.bim' %(files[0])))
            ## skip magic number and mode (and check that they are correct)
            magic_number = [108, 27]
            mode = [1]
            for j in range(len(magic_number)+len(mode)):
                byte = f_bed.read(1)
                if ord(byte) != list(magic_number+mode)[j]:
                    stop
            d['file%i' %(i)] = f_bed
            d['bim%i' %(i)] = f_bim

        ## Oxford/IMPUTE2
        elif d_args['format%i' %(i)] == 'gen':
            d['fileinput%i' %(i)] = fileinput.FileInput(
                files=files, openhook=hook_compressed_text)
            d['file%i' %(i)] = stack.enter_context(d['fileinput%i' %(i)])

        ## Beagle3
        elif d_args['format%i' %(i)] == 'bgl':
            d['fileinput%i' %(i)] = fileinput.FileInput(
                files=files, openhook=hook_compressed_text)
            d['file%i' %(i)] = stack.enter_context(d['fileinput%i' %(i)])
            if d_args['markers%i' %(i)]:
                d['markers%i' %(i)] = fileinput.FileInput(files=markers)
            else:
                d['markers%i' %(i)] = None

        ## VCF
        elif d_args['format%i' %(i)] == 'vcf':
            d['fileinput%i' %(i)] = fileinput.FileInput(
                files=files, openhook=hook_compressed_text)
            d['fileinput%i' %(i)] = itertools.chain.from_iterable(
#            d['fileinput%i' %(i)] = fileinput.input(
                [
                    ## fetch records in a region using 0-based indexing
                    pysam.TabixFile(path).fetch(
                        d_args['chrom'],
                        d_args['start']-1,
                        d_args['end']-1,
                        parser=pysam.asTuple(),
                        ) for path in files]
#                openhook=hook_compressed_text,
                )
#            d['file%i' %(i)] = stack.enter_context(d['fileinput%i' %(i)])
            d['file%i' %(i)] = d['fileinput%i' %(i)]
#            d['paths%i' %(i)] = files

#            d['file%i' %(i)] = pysam.TabixFile(files).fetch(
#                d_args['chrom'],
#                d_args['start'],
#                d_args['end'],
#                parser=pysam.asTuple(),
#                )

##            ## skip header
##            for line in d['file%i' %(i)]:
##                if line[:6] == '#CHROM':
##                    break

        ## hap
        elif d_args['format%i' %(i)] == 'hap':
            d['fileinput%i' %(i)] = fileinput.FileInput(
                files=files, openhook=hook_compressed_text)
            d['file%i' %(i)] = stack.enter_context(d['fileinput%i' %(i)])
            d['legend%i' %(i)] = fileinput.FileInput(
                files=d_args['legend%i' %(i)], openhook=hook_compressed_text)

        else:
            print(d['format%i' %(i)])
            stop

    return d


def hook_compressed_text(filename, mode, encoding='utf8'):

    ext = os.path.splitext(filename)[1]
    if ext == '.gz':
        f = gzip.open(filename, mode + 't', encoding=encoding)
    else:
        f = open(filename, mode, encoding=encoding)

    return f


def loop_main(
    affix, chrom, f_extract, f_discordant, d_samples, d_indexes, d_args,
##    file1, file2, format1, format2,
    file1=None, file2=None, format1=None, format2=None,
    ## get chromosome number from gen/bgl file name
    fileinput1=None, fileinput2=None,
    ## get chromosome and position from bim file
    bim1=None, bim2=None, markers1=None, markers2=None,
    legend1=None, legend2=None,
    ):

    ## todo: split this loooooong function into sub-processes

    d_func = {
        'gen': parse_gen, 'bgl': parse_bgl, 'bed': parse_bed,
        'vcf': parse_vcf,
        'hap': parse_hap}
    func1 = d_func[format1]
    func2 = d_func[format2]
    d_func_dosage = {
        'gen': parse_dosage_gen, 'bgl': parse_dosage_bgl,
        'bed': parse_dosage_bed, 'vcf': parse_dosage_vcf,
        'hap': parse_dosage_hap}
    func_dosage1 = d_func_dosage[format1]
    func_dosage2 = d_func_dosage[format2]

    assert len(d_indexes[1]) == len(d_indexes[2])

    ## how many bed bytes to be read for each bim line?
    n_bytes1 = math.ceil(len(d_samples[1])/4)
    n_bytes2 = math.ceil(len(d_samples[2])/4)

    kwargs1 = {
        'file': file1, 'fileinput': fileinput1,
        'bim': bim1, 'n_bytes': n_bytes1, 'markers': markers1,
        'legend': legend1, 'chrom': chrom, 'format': format1,
        'start': d_args['start'], 'end': d_args['end'],
        }
    kwargs2 = {
        'file': file2, 'fileinput': fileinput2,
        'bim': bim2, 'n_bytes': n_bytes2, 'markers': markers2,
        'legend': legend2, 'chrom': chrom, 'format': format2,
        'start': d_args['start'], 'end': d_args['end'],
        }

##    d = {'00':'1 0 0','01':'0 1 0','11':'0 0 1','10':'0.3333 0.3333 0.3333'}
    d_correl = {
        'xx': {}, 'xy': {}, 'yy': {}, 'x': {}, 'y': {},
        'n': {}, 'r2': {}, 'cnt': {}}
    for MAF in range(0, 51):
        for k in d_correl.keys():
            d_correl[k][MAF] = 0

    with open('%s.per_SNP.r2.txt' %(affix), 'w') as fd_out:
        fd_out.write('#CHROM POS MAF R2 AC\n')
        d_correl, d_stats = loop_sub(
            kwargs1, kwargs2, fd_out, func1, func2, f_extract,
            d_indexes,
            func_dosage1, func_dosage2, f_discordant, d_correl,
            d_samples,
            )

    with open('{}.per_MAF.r2.txt'.format(affix), 'w') as fd_out:
        for MAF in range(0, 51):
##            nom = (
##                d_correl['xy'][MAF] - (
##                    d_correl['x'][MAF]*d_correl['y'][MAF] /
##                    d_correl['n'][MAF]))
##            f1 = d_correl['xx'][MAF]-d_correl['x'][MAF]**2/d_correl['n'][MAF]
##            f2 = d_correl['yy'][MAF]-d_correl['y'][MAF]**2/d_correl['n'][MAF]
##            den_sq = abs(f1*f2)
##            r2 = nom**2/den_sq
##            print(MAF)
            if d_correl['cnt'][MAF] == 0:
                continue
            r2 = d_correl['r2'][MAF]/d_correl['cnt'][MAF]
##            print('\b',r2)
            fd_out.write('%s %s\n' %(MAF, r2))

    with open('{}.per_sample.r2.txt'.format(affix), 'w') as fd_out:
        for i, d in d_stats.items():
            nom = (d['sumxy'] - d['sumx'] * d['sumy'] / d['n'])
            den_sq = (d['sumxx'] - d['sumx']**2 / d['n']) * (d['sumyy'] - d['sumy']**2 / d['n'])
            try:
                r2 = nom**2 / abs(den_sq)
            except ZeroDivisionError:
                r2 = 0
            nonref_conc1 = d['nonref_disc1']/(d['nonref_conc1']+d['nonref_disc1'])
            nonref_conc2 = d['nonref_disc2']/(d['nonref_conc2']+d['nonref_disc2'])
            fd_out.write('{} {} {} {} {}\n'.format(
                d_samples[1][d_indexes[1][i]], d_samples[2][d_indexes[2][i]], r2,
                nonref_conc1, nonref_conc2,
                ))

    return


def loop_sub(
    kwargs1, kwargs2, fd_out, func1, func2, f_extract,
    d_indexes,
    func_dosage1, func_dosage2, f_discordant, d_correl,
    d_samples,
    ):

    ## parse first line to be extracted
    if f_extract:
        extract_chrom, extract_pos = next(parse_extract(f_extract))

    ## parse first line
    chrom1, pos1, A1, B1, genotypes1 = next(func1(**kwargs1))
    chrom2, pos2, A2, B2, genotypes2 = next(func2(**kwargs2))

    d_stats = {
        k: {
            'sumx': 0, 'sumy': 0, 'sumxx': 0, 'sumxy': 0, 'sumyy': 0, 'n': 0,
            'nonref_conc1': 0, 'nonref_disc1': 0, 'nonref_conc2': 0, 'nonref_disc2': 0,
            } for k in range(len(d_indexes[1]))}

    ii = 0
    while True:

        ii += 1
        if not ii % 1000:
            print(chrom1, pos1, chrom2, pos2)

        if type(chrom1) == str and type(chrom2) == int:
            break
        if type(chrom2) == str and type(chrom1) == int:
            break

        if chrom1 < chrom2:
            print('a', kwargs1['format'], kwargs2['format'], chrom1, chrom2, pos1, pos2,)
            try:
                chrom1, pos1, A1, B1, genotypes1 = next(func1(**kwargs1))
            except StopIteration:
                break
            continue
        elif chrom2 < chrom1:
#            print('b', kwargs1['format'], kwargs2['format'], chrom1, chrom2, pos1, pos2)
            try:
                chrom2, pos2, A2, B2, genotypes2 = next(func2(**kwargs2))
            except StopIteration:
                break
            continue
        elif pos1 < pos2:
            try:
                chrom1, pos1, A1, B1, genotypes1 = next(func1(**kwargs1))
#            except: break  # tmp!!!
            except StopIteration:
                break
            continue
        elif pos2 < pos1:
            try:
                chrom2, pos2, A2, B2, genotypes2 = next(func2(**kwargs2))
            except StopIteration:
                break
            continue
        elif pos1 == pos2:
            pass
        else:
            print('bim', bim_chrom, bim_pos)
            print('extract', extract_chrom, extract_pos)
            print('gen', gen_chrom, gen_pos)
            stop

        if f_extract:
            while True:
                bool_continue = False
                bool_break = False
                if extract_chrom < chrom1:
                    print('c', format1, format2, chrom1, extract_chrom)
                    try:
                        extract_chrom, extract_pos = next(
                            parse_extract(f_extract))
                    except StopIteration:
                        break
                    continue
                elif chrom1 < extract_chrom:
                    print('d', format1, format2, chrom1, extract_chrom)
                    try:
                        chrom1, pos1, A1, B1, genotypes1 = next(
                            func1(**kwargs1))
                    except StopIteration:
                        break
                    bool_continue = True
                    break
                ## SNP to be extracted not in bed file
                elif extract_pos < pos1:
                    try:
                        extract_chrom, extract_pos = next(
                            parse_extract(f_extract))
                    except StopIteration:
                        break
                    continue
                ## SNP in bed file not to be extracted
                elif pos1 < extract_pos:
                    try:
                        chrom1, pos1, A1, B1, genotypes1 = next(
                            func1(**kwargs1))
                    except StopIteration:
                        bool_break = True
                        break
                    bool_continue = True
                    break
                ## pos1 == extract_pos
                else:
                    if bool_continue != False or bool_break != False:
                        print(
                            pos1, pos2, extract_pos, bool_continue, bool_break)
                        stop
                    break
            if bool_continue == True:
                continue
            if bool_break == True:
                break

        bool_continue = False
        ## skip insertion or deletion
        if len(B1) > 1 or len(A1) > 1:
            try:
                chrom1, pos1, A1, B1, genotypes1 = next(func1(**kwargs1))
            except StopIteration:
                break
            continue
        elif len(B2) > 1 or len(A2) > 1:
            try:
                chrom2, pos2, A2, B2, genotypes2 = next(func2(**kwargs2))
            except StopIteration:
                break
            continue
        elif (A1 == A2 or A1 == '0' or A2 == '0') and B1 == B2:
            bool_reverse = False
        elif (A1 == B2 or A1 == '0' or B2 == '0') and B1 == A2:
            bool_reverse = True
        else:
            try:
                chrom1, pos1, A1, B1, genotypes1 = next(func1(**kwargs1))
                chrom2, pos2, A2, B2, genotypes2 = next(func2(**kwargs2))
            except StopIteration:
                break
            continue

        ##
        ##
        ##
        (
            n, sum_x, sum_y, sum_xx, sum_yy, sum_xy,
#            cnt_concordance,
            ) = calc_sums_and_sum_of_squares(
                d_indexes,
                func_dosage1, func_dosage2, genotypes1, genotypes2,
                bool_reverse,
                f_discordant, d_stats,
                chrom1, pos1, A1, B1,
                chrom2, pos2, A2, B2,
                d_samples,
                )

        ## All genotypes are missing, because no QC was carried out.
        ## Or monomorphic in both sets.
        if n == 0 or (sum_x == 0 and sum_y == 0):
            try:
                chrom1, pos1, A1, B1, genotypes1 = next(func1(**kwargs1))
                chrom2, pos2, A2, B2, genotypes2 = next(func2(**kwargs2))
                if f_extract:
                    extract_chrom, extract_pos = next(parse_extract(f_extract))
                continue
            except StopIteration:
                break

        ## calculate MAF
        MAF_x = sum2maf(sum_x, n)
        MAF_y = sum2maf(sum_y, n)
        AC_x = sum_x
        AC_y = sum_y
##        MAF = math.ceil(MAF)

        ## append
        d_correl['x'][MAF_x] += sum_x
        d_correl['y'][MAF_x] += sum_y
        d_correl['xx'][MAF_x] += sum_xx
        d_correl['xy'][MAF_x] += sum_xy
        d_correl['yy'][MAF_x] += sum_yy
        d_correl['n'][MAF_x] += n

        ## calculate r2
        nom = (sum_xy-sum_x*sum_y/n)
        den_sq = (sum_xx-sum_x**2/n)*(sum_yy-sum_y**2/n)
        ## ask Deepti whether to include or exclude these SNPs
        ## if calculating overall r2 with new method
        ## instead of simple average
        if den_sq != 0:
            r2 = nom**2/abs(den_sq)
        else:
##            if (
##                nom == 0 and sum_x ==0 and sum_y == 0
##                and sum_xx == 0 and sum_yy == 0 and sum_xy == 0:
            if False:
                r2 = 1
            else:
                r2 = None
                ## One monomorphic, but the other not, so we can't calculate correlation.
                if (sum_x == 0 and sum_y > 0) or (sum_x > 0 and sum_y == 0):
                    pass
                ## Both monomorphic, so correlation can't be calculated.
                elif sum_x == 0 and sum_y == 0:
                    pass
                ## Both monomorphic.
                elif sum_x == 2*n and sum_y == 2*n:
                    pass
                ## One monomorphic.
                elif sum_x == 2*n or sum_y == 2*n:
                    pass
                ## One completely heterozygous.
                elif (sum_x == n and sum_xx == n) or (sum_y == n and sum_yy == n):
                    pass
                elif round(nom, 6) == 0:
                    print(nom)
                    print(sum_x, sum_y, sum_xy, sum_xx, sum_yy)
                    print(n, 2*n)
                    print(genotypes1)
                    print(genotypes2)
                    stoptmp
                    pass
                ## monomorphic
#                elif l_x == len(l_x)*[x] or l_y == len(l_y)*[y]:
#                    pass
                elif sum_x == 2*n or sum_y == 2*n:
                    stoptmp
                    pass
                elif (sum_x == 0 or sum_y == 0):
                    stop
                else:
                    print(genotypes1)
                    print(genotypes2)
                    print(sum_x, sum_y, sum_xx, sum_xy, sum_yy)
                    print(nom, den_sq)
                    print('x', l_x)
                    print('y', l_y)
                    print(pos1, A1, A2, B1, B2)

                    print(sum_x, sum_y, sum_xy)
                    print(nom, round(nom, 13), 13)
                    for x2 in l_x:
                        if x != x2:
                            print('x', x2)
                    for y2 in l_y:
                        if y != y2:
                            print('y', y2)
                    stop
##        elif nom == 0:
##            print(nom,den_sq)
##            print(sum_x,sum_y,n)
##            print(l_x)
##            print(l_y)
##            stop2
            pass


        if r2 is not None:
            fd_out.write('%s %s %s %s %s %s %s\n' %(chrom1, pos1, MAF_x, r2, AC_x, MAF_y, AC_y))
            d_correl['r2'][MAF_x] += r2
            d_correl['cnt'][MAF_x] += 1
            pass

        try:
            chrom1, pos1, A1, B1, genotypes1 = next(func1(**kwargs1))
            chrom2, pos2, A2, B2, genotypes2 = next(func2(**kwargs2))
            if f_extract:
                extract_chrom, extract_pos = next(parse_extract(f_extract))
        except StopIteration:
            break

        ## continue loop over SNPs
        continue

    return d_correl, d_stats


def sum2maf(sum_x, n):

    try:
        MAF = 50*sum_x/n
    except:
        print(genotypes1)
        print(genotypes2)
        MAF = 50*sum_x/n
    if MAF > 50:
        MAF = 100-MAF
    MAF = int(MAF)
        
    return MAF


def calc_sums_and_sum_of_squares(
    d_indexes,
    func_dosage1, func_dosage2, genotypes1, genotypes2, bool_reverse,
    f_discordant, d_stats,
    ## variables passed for writing to discordance file only
    chrom1, pos1, A1, B1,
    chrom2, pos2, A2, B2,
    d_samples,
    ):

    ## This is the slowest function followed by functions for parsing GTs.

    n = 0
    sum_x = 0
    sum_y = 0
    sum_xx = 0
    sum_xy = 0
    sum_yy = 0
##    l_x = [] ## tmp!!!
##    l_y = [] ## tmp!!!
#    cnt_concordance = 0
    for i, (index1, index2) in enumerate(zip(d_indexes[1], d_indexes[2])):
##        for i_byte in range(n_bytes):
##                byte = bed.read(1)
        x = func_dosage1(genotypes1, index1)
        y = func_dosage2(genotypes2, index2)
        if x == None or y == None:
            continue
        if bool_reverse:
            x = 2-x
        if f_discordant and x != y:
            print(
                chrom1, pos1, A1, B1,
                chrom2, pos2, A2, B2,
                '{0:5.3f}'.format(x), '{0:5.3f}'.format(y),
                d_samples[1][index1], d_samples[2][index2],
                file=f_discordant, flush=True, sep='\t')
##        l_x += [x] ## tmp!!!
##        l_y += [y] ## tmp!!!
        xx = x**2
        xy = x*y
        yy = y**2
        sum_x += x
        sum_y += y
        sum_xx += xx
        sum_xy += xy
        sum_yy += yy
        n += 1
        ## per sample stats
        d_stats[i]['sumx'] += x
        d_stats[i]['sumy'] += y
        d_stats[i]['sumxx'] += xx
        d_stats[i]['sumxy'] += xy
        d_stats[i]['sumyy'] += yy
        d_stats[i]['n'] += 1
        ## For the time being assume that the truth set if file set 2... Temporary!!!
        ## Actually do both ways!
        if round(y,0) != 0:
            if round(x,0) == round(y,0):
                d_stats[i]['nonref_conc1'] += 1
            else:
                d_stats[i]['nonref_disc1'] += 1
        if round(x,0) != 0:
            if round(x,0) == round(y,0):
                d_stats[i]['nonref_conc2'] += 1
            else:
                d_stats[i]['nonref_disc2'] += 1

#        if x == y:
#            cnt_concordance += 1
        ## continue loop over samples
        continue

#    return n, sum_x, sum_y, sum_xx, sum_yy, sum_xy, cnt_concordance
    return n, sum_x, sum_y, sum_xx, sum_yy, sum_xy


def parse_hap(**kwargs):

    file = kwargs['file']
#    fileinput = kwargs['fileinput']
    legend = kwargs['legend']
    line_hap = file.readline()
    line_legend = legend.readline()
    ## EOF
    if line_legend == '':
        return
    ## skip header
    if line_legend.split(' ', 1)[0] == 'id':
        line_legend = legend.readline()
    l_legend = line_legend.rstrip().split()
##    chrom, pos = l_legend[0].split(':')
##    chrom = int(chrom)
    try:
        pos = int(l_legend[1])
    except:
        print(l_legend)
        stop
    alleleA = l_legend[2]
    alleleB = l_legend[3]
    l = line_hap.rstrip().split()
    if kwargs['chrom']:
        chrom = kwargs['chrom']
    else:
        ## 0 is a magic number...
        ## instead provide a file with chromosome-file-equivalents
        chrom = int(
            os.path.basename(file.filename()).split('.')[0])  # tmp!!!
##        chrom = int(l_legend[0].split(':')[0])

    yield chrom, pos, alleleA, alleleB, l


def parse_dosage_hap(l_hap, i_sample):

    gt1 = l_hap[2*i_sample]
    gt2 = l_hap[2*i_sample+1]
    if gt1 == '0' and gt2 == '0':
        dosage = 0
    elif gt1 == '0' and gt2 == '1':
        dosage = 1
    elif gt1 == '1' and gt2 == '1':
        dosage = 2
    elif gt1 == '1' and gt2 == '0':
        dosage = 1
    else:
        print(gt1, gt2)
        print(i_sample)
        print(len(l_hap))
        stop

    return dosage


def parse_dosage_vcf(l_vcf, i_sample):

    ## assume biallelic variant...

    ## alternatively use re.split()
    ## re.split('\| |, /,str)

##    dosage = sum(
##        int(GT) for GT in l_vcf[i_sample+9].split(':')[0].split('|'))

##    dosage = sum(i*float(GP) for i,GP in enumerate(
##        l_vcf[i_sample+9].split(':')[2].split(',')))

    ## It's pretty fucking stupid doing this for each SNP and for each sample!
    for k in ('DS', 'PL', 'GL'):
        try:
            index = l_vcf[8].split(':').index(k)
            break
        except ValueError:
            continue
    else:
        assert l_vcf[8] == 'GT'
        k = 'GT'
        index = 0

    l_vcf = l_vcf[i_sample+9].split(':')
    if l_vcf[0] == './.':
        dosage = 1
    else:
        if k == 'PL':
            PL = l_vcf[index].split(',')
            if PL == ['.']:
                ## set dosage to zero if PL missing
                dosage = 0
            else:
                ## assume biallelic variant
                try:
                    assert len(PL) == 3
                except:
                    print(l_vcf, PL)
                    stop
    ##        dosage = sum(
    ##            i*prob for i, prob in enumerate(
    ##                pow(10,-int(log10likelihood)/10)
    ##                for log10likelihood in PL))
                ## calculate probabilities once
                l_probs = [
                    pow(10, -int(log10likelihood)/10) for log10likelihood in PL]
                ## evaluate sum once
                sum_prob = sum(l_probs)
                dosage = sum(i*prob/sum_prob for i, prob in enumerate(l_probs))
        elif k == 'DS':
            dosage = float(l_vcf[index])
        elif k == 'GL':
            GL = l_vcf[index].split(',')
            ## assume biallelic variant
            try:
                assert len(GL) == 3
            except:
                print(l_vcf)
                stop
            ## calculate probabilities once
            l_probs = [
                pow(10, float(log10likelihood)) for log10likelihood in GL]
            ## evaluate sum once
            sum_prob = sum(l_probs)
            dosage = sum(i*prob/sum_prob for i, prob in enumerate(l_probs))
        ## Only GTs provided.
        else:
            dosage = sum((int(_) for _ in re.split(r'\|', l_vcf[index])))

    return dosage


def parse_dosage_gen(l_gen, i_sample):

    dosage = 0
    dosage += float(l_gen[6+3*i_sample])
    dosage += 2*float(l_gen[7+3*i_sample])

    return dosage


def parse_dosage_bgl(l_bgl, i_sample):

    dosage = 0
    dosage += float(l_bgl[4+3*i_sample])
    dosage += 2*float(l_bgl[5+3*i_sample])

    return dosage


def parse_dosage_bed(bytesSNP, i_sample):

    d_dosage = {'00': 0., '01': 1., '11': 2.}

    i_byte = int(i_sample/4)
    byte = bytesSNP[i_byte]
    s8 = str(bin(byte)[2:]).zfill(8)
    if len(s8) != 8:
        stop1
    i = 3-i_sample%4
    s2 = s8[2*i+1]+s8[2*i]
    if s2 == '10':
        dosage = None  # ask Deepti whether missing should be het or N/A!!!
    else:
        dosage = d_dosage[s2]

    return dosage


def index_samples(d_args):

    ## this function is a COMPLETE mess!!! CLEAN IT UP!!!

    d_samples = {}
    d_cnt = {}
    for i in 1, 2:
        fileformat = d_args['format{}'.format(i)]
        files = d_args['files%i' %(i)]
        if fileformat == 'bed':
            l_samples = fam2samples('%s.fam' %(files[0]), d_args)
        elif fileformat == 'gen':
            if os.path.getsize(files[-1]) == 0:
                print('size 0', files[0])
                sys.exit()
            if files[-1][-3:] == '.gz':
                with gzip.open(files[-1], 'rt') as f:
                    n_fam = int((len(f.readline().rstrip().split())-5)/3)
            else:
                with open(files[-1]) as f:
                    n_fam = int((len(f.readline().rstrip().split())-5)/3)
            if d_args['fam%i' %(i)]:
                l_samples = fam2samples(d_args['fam%i' %(i)], d_args)
                if len(l_samples) != n_fam:
                    print(l_samples)
                    print(len(l_samples))
                    print(n_fam)
                    stop
            elif d_args['sample%i' %(i)]:
                with open(d_args['sample%i' %(i)]) as f:
                    l_samples = [
                        line.split()[0] for line in f.readlines()[2:]]
            else:
                l_samples = None
        elif fileformat == 'bgl':
            if files[0][-3:] == '.gz':
                with gzip.open(files[0], 'rt') as f:
                    l_samples = f.readline().rstrip().split()[3:-1:3]
            else:
                with open(files[0]) as f:
                    l_samples = f.readline().rstrip().split()[3:-1:3]
        elif fileformat == 'vcf':
            l_samples = vcf2samples(files[0])
        elif fileformat == 'hap':
            if d_args['sample%i' %(i)]:
                with open(d_args['sample{}'.format(i)]) as f:
                    l_samples = [
                        line.split()[1] for line in f.readlines()[2:]]
            else:
                l_samples = None
        else:
            print(fileformat, files)
            stop_unknown_format
        d_samples[i] = l_samples
        if l_samples:
            d_cnt[i] = len(l_samples)
        else:
            d_cnt[i] = n_fam

    d_indexes = {1: [], 2: []}
    if d_args['sampledic']:
        with open(d_args['sampledic']) as lines:
            d_sampledic = {
                line.split()[0]: line.rstrip().split()[1] for line in lines}
    else:
        if d_samples[1] and d_samples[2]:
            d_sampledic = {sample: sample for sample in set(
                set(d_samples[1]) & set(d_samples[2]))}
        elif d_cnt[1] == d_cnt[2]:
            if not d_samples[1]:
                d_samples[1] = d_samples[2]
            else:
                print(d_samples[1])
                print(d_samples[2])
                print(d_cnt)
                stoptmp
        else:
            print(d_cnt)
            print(len(d_samples[1]))
            print(len(d_samples[2]))
            stop_unknown_sample_IDs
        d_sampledic = {sample: sample for sample in d_samples[1]}

    bool_error = True
    samples1 = set(d_samples[1])
    samples2 = set(d_samples[2])
    intersect1 = set()
    intersect2 = set()
    for index1 in range(len(d_samples[1])):
        sample1 = d_samples[1][index1]
        try:
            index2 = d_samples[2].index(d_sampledic[sample1])
            intersect1.add(sample1)
            intersect2.add(d_sampledic[sample1])
        ## sample not at intersection
        except ValueError:
            print('#not found in sample dic', sample1)
            continue
        ## not in sample dic
        except KeyError:
            print('#not intersect', sample1)
            continue
        d_indexes[1] += [index1]
        d_indexes[2] += [index2]
        bool_error = False

    if bool_error:
        print(d_args['sampledic'])
        print(d_args['files1'][0])
        print(d_args['files2'][0])
        print(d_samples[1])
        print(d_samples[2])
        print(d_indexes[1])
        print(d_indexes[2])
        print('0 samples intersect')
        raise('0 samples intersect')
        sys.exit()

    print(
        'indexes', len(d_indexes[1]), len(d_indexes[2]),
        'samples', len(d_samples[1]), len(d_samples[2]),
        '\ncorresponding samples',
        '\n', d_samples[1][d_indexes[1][0]], d_samples[2][d_indexes[2][0]],
        '\n', d_samples[1][d_indexes[1][-1]], d_samples[2][d_indexes[2][-1]],
        )

    print('not found', len(samples1-intersect1), len(samples2-intersect2))
    print(list(samples1-intersect1)[:10])
    print(list(samples2-intersect2)[:10])

    print('\ncorresponding indeces (first ten in each file)')
    print(d_indexes[1][:10], d_indexes[2][:10])

    return d_samples, d_indexes


def vcf2samples(file):

    if file[-3:] == '.gz':
        f = gzip.open(file, 'rt')
    else:
        f = open(file)
    for line in f:
        if line[0] != '#':
            break
        header = line
    f.close()
    l_samples = header.rstrip().split('\t')[9:]

    return l_samples


def fam2samples(fam, d_args):

    set_samples = set()
    if d_args['sampledic']:
        with open(d_args['sampledic']) as f:
            for line in f:
                set_samples.add(line.split()[0])
                set_samples.add(line.split()[1])
    l_samples = []
    keyword = re.compile(r'(\d\d\d\d\d\d_[A-H]\d\d_)(.+\d\d\d\d\d\d\d)')
    with open(fam, 'r') as lines:
        for line in lines:
            ## Use IID rather than FID.
            sampleID = line.split()[1]
            match = result = keyword.search(sampleID)
            if match and sampleID not in set_samples:
                l_samples += [match.group(2)]
            else:
                l_samples += [sampleID]

    return l_samples


def parse_vcf(**kwargs):

#    vcf = kwargs['file']
    for l in kwargs['file']:
#    line = vcf.readline()
#        if line == '':
#            return
#    while line[0] == '#':
#        line = vcf.readline()
#        l = line.rstrip().split('\t')
#        if line[0] == '#':
#            continue
#        l = line.rstrip().split('\t')
        if l[0].startswith('chr'):
            l[0] = l[0][3:]
        try:
            chrom = int(l[0])
        except:
            chrom = l[0]
        pos = int(l[1])
        alleleA = l[3]
        alleleB = l[4]

        yield chrom, pos, alleleA, alleleB, l


def parse_bed(**kwargs):

    ## http://pngu.mgh.harvard.edu/~purcell/plink/binary.shtml

    bim = kwargs['bim']
    bed = kwargs['file']
    n_bytes = kwargs['n_bytes']

    line_bim = bim.readline()
    if line_bim == '':
        return
    l_bim = line_bim.rstrip().split()
    chrom = int(l_bim[0])
    if chrom == 23:
        chrom = 'X'
##        StopIteration
    pos = int(l_bim[3])
    alleleA = l_bim[4]
    alleleB = l_bim[5]

    ## read bed bytes corresponding to bim line
    bytesSNP = bed.read(n_bytes)

    yield chrom, pos, alleleA, alleleB, bytesSNP


def parse_extract(f_extract):

    line_extract = f_extract.readline()
    if line_extract == '':
        return
    l = line_extract.split()
    chrom = int(l[0])
    pos = int(l[1])

    yield chrom, pos


def parse_bgl(**kwargs):

    file = kwargs['file']
    fileinput = kwargs['fileinput']
    markers = kwargs['markers']

    line = file.readline()
    if line == '':
        return  # StopIteration
    if fileinput.isfirstline():
        print(fileinput.filename())
        line = file.readline()
    l = line.rstrip().split()

    if markers:
        line_markers = markers.readline()
        l_markers = line_markers.split()
        while any([
            l_markers[0] != l[0],
            l_markers[2] != l[1],
            l_markers[3] != l[2]]):
            line_markers = markers.readline()
            l_markers = line_markers.split()
        if any([
            l_markers[0] != l[0],
            l_markers[2] != l[1],
            l_markers[3] != l[2]]):
            print(line[:80])
            print(line_markers)
            print(fileinput.filename())
            print(markers.filename())
            stop
        chrom = 1  # tmp!!!
        pos = int(l_markers[1])
    else:
        l0 = l[0].split(':')
        chrom = int(l0[0])
        pos = int(l0[1])

    alleleA = l[1]
    alleleB = l[2]

    yield chrom, pos, alleleA, alleleB, l


def parse_gen(**kwargs):

    file = kwargs['file']
    fileinput = kwargs['fileinput']

    line = file.readline()
    if line == '':
        return
    l = line.rstrip().split()
    pos = int(l[2])
#    chrom = int(fileinput.filename().split('/')[1])  # tmp!!!
    ## -2 is a magic number...
    ## instead provide a file with chromosome-file-equivalents
    alleleA = l[3]
    alleleB = l[4]
    if kwargs['chrom']:
        chrom = kwargs['chrom']
    else:
        try:
            chrom = int(l[1].split(':')[0])
        except ValueError:
            try:
                chrom = int(fileinput.filename().split('/')[-2])
            except ValueError:
                try:
                    chrom = int(
                        fileinput.filename().split('/')[-1].split('.')[0])
                except ValueError:
                    chrom = int(
                        fileinput.filename().split('/')[-1].split('.')[0][3:])

    yield chrom, pos, alleleA, alleleB, l


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

    parser.add_argument('--files1', '--file1', required=True, nargs='+')
    parser.add_argument('--files2', '--file2', required=True, nargs='+')
    choices = ['bed', 'gen', 'bgl', 'vcf', 'hap']
    parser.add_argument('--format1', required=True, choices=choices)
    parser.add_argument('--format2', required=True, choices=choices)
    ## type=argparse.FileType('rt')
    parser.add_argument('--fam1', help='PLINK format')
    parser.add_argument('--fam2', help='PLINK format')
    parser.add_argument(
        '--affix', '--out', '--prefix', '--suffix', required=True)
    parser.add_argument('--sampledic', help='convert chip APP IDs to sequence EGAN IDs')
    parser.add_argument('--markers1', nargs='+', default=[])
    parser.add_argument('--markers2', nargs='+', default=[])
    ## hap, legend, sample
    parser.add_argument('--legend1', nargs='+', default=[])
    parser.add_argument('--legend2', nargs='+', default=[])
    parser.add_argument('--sample1', help='IMPUTE2 format')
    parser.add_argument('--sample2', help='IMPUTE2 format')
    ## file2chrom
    parser.add_argument('--file2chrom')
    ## chrom
    parser.add_argument('--chrom')
    parser.add_argument('--start', type=int)
    parser.add_argument('--end', type=int)
    ## optional
    parser.add_argument('--extract', help='extract only selected SNPs')
    parser.add_argument(
        '--discordant', help='print discordant samples to file')

    namespace_args = parser.parse_args()

    d_args = vars(namespace_args)

    for k, v in d_args.items():
        if v == 'None':
            d_args[k] = None

    for key in 'files1', 'files2', 'markers1', 'markers2':
        d_args[key] = sort_nicely(d_args[key])

    try:
        d_args['chrom'] = int(d_args['chrom'])
    except:
        pass

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
