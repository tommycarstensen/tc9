#!/bin/python3

## Tommy Carstensen, Wellcome Trust Sanger Institute, October 2013

## http://pngu.mgh.harvard.edu/~purcell/plink/binary.shtml

## todo: make it possible to read multiple bed files
## todo: give the user an option to submit a file that links chromosomes and IMPUTE2 output files --- files should then be sorted by chromosome

## gen,bed,bgl,vcf vs gen,bed,bgl,vcf

import argparse
import os
import math
import contextlib
import fileinput
import re
import sys
import gzip


def main():

    (
        files1,files2,format1,format2,fam1,fam2,
        fp_extract,affix,fp_sampledic) = argparser()

    d_samples,d_indexes = index_samples(
        files1,format1,files2,format2,fam1,fam2,fp_sampledic)

    with contextlib.ExitStack() as stack:

        d = open_files(stack,files1,files2,format1,format2)

        if fp_extract:
            f_extract = stack.enter_context(open(fp_extract))
        else:
            f_extract = None

        loop(affix,f_extract,d_samples,d_indexes,**d)

    return


def open_files(stack,files1,files2,format1,format2):

    ## tidy up this function and delete duplicate code!!!

    ## http://docs.python.org/3/tutorial/controlflow.html
    d = {'format1':format1,'format2':format2}
    d_files = {1:files1,2:files2}

    for i in (1,2):
        files = d_files[i]

        if d['format%i' %(i)] == 'bed':
            if len(files) != 1:
                print(files)
                stoptmp
            f_bed = stack.enter_context(open('%s.bed' %(files[0]),'rb'))
            f_bim = stack.enter_context(open('%s.bim' %(files[0])))
            ## skip magic number and mode (and check that they are correct)
            magic_number = [108,27]
            mode = [1]
            for i in range(len(magic_number)+len(mode)):
                byte = f_bed.read(1)
                if ord(byte) != list(magic_number+mode)[i]:
                    stop
            d['file%i' %(i)] = f_bed
            d['bim%i' %(i)] = f_bim

        elif d['format%i' %(i)] == 'gen':
            d['fileinput%i' %(i)] = fileinput.FileInput(files=files)
            d['file%i' %(i)] = stack.enter_context(d['fileinput%i' %(i)])

        elif d['format%i' %(i)] == 'bgl':
            d['fileinput%i' %(i)] = fileinput.FileInput(files=files)
            d['file%i' %(i)] = stack.enter_context(d['fileinput%i' %(i)])
##            ## skip header (not really necessary as this is done in parse_bgl)
##            for line in d['file%i' %(i)]:
##                break

        elif d['format%i' %(i)] == 'vcf':
            d['fileinput%i' %(i)] = fileinput.FileInput(
                files=files,openhook=hook_compressed_text)
            d['file%i' %(i)] = stack.enter_context(d['fileinput%i' %(i)])
##            ## skip header
##            for line in d['file%i' %(i)]:
##                if line[:6] == '#CHROM':
##                    break
            
        else:
            print(d['format%i' %(i)])
            stop

    return d


def hook_compressed_text(filename, mode, encoding='utf8'):

    ##http://stackoverflow.com/questions/21529163/python-gzipped-fileinput-returns-binary-string-instead-of-text-string/21529243

    ext = os.path.splitext(filename)[1]
    if ext == '.gz':
        f = gzip.open(filename, mode + 't', encoding=encoding)
    else:
        f = open(filename, mode, encoding=encoding)

    return f
    

def loop(
    affix,f_extract,d_samples,d_indexes,
    file1,file2,format1,format2,
    ## get chromosome number from gen/bgl file name
    fileinput1=None,fileinput2=None,
    ## get chromosome and position from bim file
    bim1=None,bim2=None,
    ):

    d_func = {'gen':parse_gen,'bgl':parse_bgl,'bed':parse_bed,'vcf':parse_vcf,}
    func1 = d_func[format1]
    func2 = d_func[format2]
    d_func_dosage = {
        'gen':parse_dosage_gen, 'bgl':parse_dosage_bgl,
        'bed':parse_dosage_bed, 'vcf':parse_dosage_vcf,}
    func_dosage1 = d_func_dosage[format1]
    func_dosage2 = d_func_dosage[format2]

    n_samples_intersection = len(d_indexes[1])

    ## how many bed bytes to be read for each bim line?
    n_bytes1 = math.ceil(len(d_samples[1])/4)
    n_bytes2 = math.ceil(len(d_samples[2])/4)

    kwargs1 = {'file':file1,'fileinput':fileinput1,'bim':bim1,'n_bytes':n_bytes1}
    kwargs2 = {'file':file2,'fileinput':fileinput2,'bim':bim2,'n_bytes':n_bytes2}

##    d = {'00':'1 0 0','01':'0 1 0','11':'0 0 1','10':'0.3333 0.3333 0.3333'}
    d_correl = {'xx':{},'xy':{},'yy':{},'x':{},'y':{},'n':{},'r2':{},'cnt':{}}
    for MAF in range(0,51):
        for k in d_correl.keys():
            d_correl[k][MAF] = 0

    if f_extract:
        extract_chrom,extract_pos = next(parse_extract(f_extract))

    ## parse first line
    chrom1,pos1,alleleA1,alleleB1,genotypes1 = next(func1(**kwargs1))
    chrom2,pos2,alleleA2,alleleB2,genotypes2 = next(func2(**kwargs2))

    fd_out = open('%s.tmp' %(affix),'w')

    while True:

        if chrom1 < chrom2:
            print('a',format1,format2,chrom1,chrom2,pos1,pos2,)
            try:
                chrom1,pos1,alleleA1,alleleB1,genotypes1 = next(func1(**kwargs1))
            except StopIteration:
                break
            continue
        elif chrom2 < chrom1:
##            print('b',format1,format2,chrom1,chrom2,pos1,pos2)
            try:
                chrom2,pos2,alleleA2,alleleB2,genotypes2 = next(func2(**kwargs2))
            except StopIteration:
                break
            continue
        elif pos1 < pos2:
            try:
                chrom1,pos1,alleleA1,alleleB1,genotypes1 = next(func1(**kwargs1))
            except StopIteration:
                break
            continue
        elif pos2 < pos1:
##            print(4,gen_pos,bim_pos) ## tmp!!!
            try:
                chrom2,pos2,alleleA2,alleleB2,genotypes2 = next(func2(**kwargs2))
            except StopIteration:
                break
            continue
        elif pos1 == pos2:
            pass
        else:
            print('bim',bim_chrom,bim_pos)
            print('extract',extract_chrom,extract_pos)
            print('gen',gen_chrom,gen_pos)
            stop

        if f_extract:
            while True:
                bool_continue = False
                bool_break = False
                if extract_chrom < chrom1:
                    print('c',format1,format2,chrom1,extract_chrom)
                    try:
                        extract_chrom,extract_pos = next(parse_extract(f_extract))
                    except StopIteration:
                        break
                    continue
                elif chrom1 < extract_chrom:
                    print('d',format1,format2,chrom1,extract_chrom)
                    try:
                        chrom1,pos1,alleleA1,alleleB1,genotypes1 = next(func1(**kwargs1))
                    except StopIteration:
                        break
                    bool_continue = True
                    break
                ## SNP to be extracted not in bed file
                elif extract_pos < pos1:
                    try:
                        extract_chrom,extract_pos = next(parse_extract(f_extract))
                    except StopIteration:
                        break
                    continue
                ## SNP in bed file not to be extracted
                elif pos1 < extract_pos:
                    try:
                        chrom1,pos1,alleleA1,alleleB1,genotypes1 = next(func1(**kwargs1))
                    except StopIteration:
                        bool_break = True
                        break
                    bool_continue = True
                    break
                ## pos1 == extract_pos
                else:
                    if bool_continue != False or bool_break != False:
                        print(pos1,pos2,extract_pos,bool_continue,bool_break)
                        stop
                    break
            if bool_continue == True:
                continue
            if bool_break == True:
                break
                
        bool_continue = False
        ## skip insertion or deletion
        if len(alleleB1) > 1 or len(alleleA1) > 1:
            try:
                chrom1,pos1,alleleA1,alleleB1,genotypes1 = next(func1(**kwargs1))
            except StopIteration:
                break
            continue
        elif len(alleleB2) > 1 or len(alleleA2) > 1:
            try:
                chrom2,pos2,alleleA2,alleleB2,genotypes2 = next(func2(**kwargs2))
            except StopIteration:
                break
            continue
        elif (alleleA1 == alleleA2 or alleleA1 == '0' or alleleA2 == '0') and alleleB1 == alleleB2:
            bool_reverse = False
        elif (alleleA1 == alleleB2 or alleleA1 == '0' or alleleB2 == '0') and alleleB1 == alleleA2:
            bool_reverse = True
        else:
            try:
                chrom1,pos1,alleleA1,alleleB1,genotypes1 = next(func1(**kwargs1))
                chrom2,pos2,alleleA2,alleleB2,genotypes2 = next(func2(**kwargs2))
            except StopIteration:
                break
            continue

        n = 0
        sum_x = 0
        sum_y = 0
        sum_xx = 0
        sum_xy = 0
        sum_yy = 0
        l_x = [] ## tmp!!!
        l_y = [] ## tmp!!!
        cnt_concordance = 0
        for index in range(n_samples_intersection):
            index1 = d_indexes[1][index]
            index2 = d_indexes[2][index]
##        for i_byte in range(n_bytes):
##                byte = bed.read(1)
            x = func_dosage1(genotypes1,index1)
            y = func_dosage2(genotypes2,index2)
            if x == None or y == None:
                continue
            if bool_reverse:
                x = 2-x
            l_x += [x] ## tmp!!!
            l_y += [y] ## tmp!!!
            sum_x += x
            sum_y += y
            sum_xx += x*x
            sum_xy += x*y
            sum_yy += y*y
            n += 1
            if x == y: cnt_concordance += 1
            ## continue loop over samples
            continue

        ## calculate MAF
        try:
            MAF = 50*sum_x/n
        except:
            print(genotypes1)
            print(genotypes2)
            MAF = 50*sum_x/n
        if MAF > 50:
            MAF = 100-MAF
        MAF = int(MAF)
##        MAF = math.ceil(MAF)

        ## append
        d_correl['x'][MAF] += sum_x
        d_correl['y'][MAF] += sum_y
        d_correl['xx'][MAF] += sum_xx
        d_correl['xy'][MAF] += sum_xy
        d_correl['yy'][MAF] += sum_yy
        d_correl['n'][MAF] += n

        ## calculate r2
        nom = (sum_xy-sum_x*sum_y/n)
        den_sq = (sum_xx-sum_x**2/n)*(sum_yy-sum_y**2/n)
        if den_sq == 0: ## ask Deepti whether to include or exclude these SNPs if calculating overal r2 with new method instead of simple average
            r2 = None
            if round(nom,13) == 0:
                pass
            ## monomorphic
            elif l_x == len(l_x)*[x] or l_y == len(l_y)*[y]:
                pass
            elif sum_x == 2*n or sum_y == 2*n:
                pass
            elif (sum_x == 0 or sum_y == 0):
                stop
            else:
                print('x',l_x)
                print('y',l_y)
                print(pos1,alleleA1,alleleA2,alleleB1,alleleB2)

                print(sum_x,sum_y,sum_xy)
                print(nom,round(nom,13),13)
                for x2 in l_x:
                    if x != x2:
                        print('x',x2)
                for y2 in l_y:
                    if y != y2:
                        print('y',y2)
                stop
##        elif nom == 0:
##            print(nom,den_sq)
##            print(sum_x,sum_y,n)
##            print(l_x)
##            print(l_y)
##            stop2
        else:
            r2 = nom**2/abs(den_sq)
            fd_out.write('%s %s %s %s\n' %(chrom1,pos1,MAF,r2))
##            print(r2,MAF)
            d_correl['r2'][MAF] += r2
            d_correl['cnt'][MAF] += 1

##            print(bim_chrom,bim_pos,MAF,r2)
##            print(chrom1,chrom2,pos1,pos2,MAF,r2)
##            import numpy
##            xxx = numpy.corrcoef(l_x,l_y)
##            print(xxx[0][1]**2)

##            print(r2,format1,format2)
##            print(l_x)
##            print(l_y)
            pass

##        print(r2,format1,format2)

##        fd_out.write('%s %s %s %s\n' %(chrom1,pos1,MAF,cnt_concordance/100))

##        if pos1 == 846338:
##            print(l_x)
##            print(l_y)
##            print(r2)
##            print(n,len(l_x),len(l_y))
##            stop

        try:
            chrom1,pos1,alleleA1,alleleB1,genotypes1 = next(func1(**kwargs1))
            chrom2,pos2,alleleA2,alleleB2,genotypes2 = next(func2(**kwargs2))
            if f_extract:
                extract_chrom,extract_pos = next(parse_extract(f_extract))
        except StopIteration:
            break

        if pos1 % 10000 == 0: print(chrom1,pos1) ## tmp!!!

        ## continue loop over SNPs
        continue

    fd_out.close()

    with open('%s.r2' %(affix),'w') as fd_out:
        for MAF in range(0,51):
##            nom = (d_correl['xy'][MAF]-d_correl['x'][MAF]*d_correl['y'][MAF]/d_correl['n'][MAF])
##            f1 = d_correl['xx'][MAF]-d_correl['x'][MAF]**2/d_correl['n'][MAF]
##            f2 = d_correl['yy'][MAF]-d_correl['y'][MAF]**2/d_correl['n'][MAF]
##            den_sq = abs(f1*f2)
##            r2 = nom**2/den_sq
            print(MAF)
            if d_correl['cnt'][MAF] == 0:
                continue
            r2 = d_correl['r2'][MAF]/d_correl['cnt'][MAF]
            print('\b',r2)
            fd_out.write('%s %s\n' %(MAF,r2))

    return


def parse_dosage_vcf(l_vcf,i_sample):

    ## alternatively use re.split()
    ## re.split('\| |, /,str)
##    dosage = sum(int(GT) for GT in l_vcf[i_sample+9].split(':')[0].split('|'))
    dosage = sum(i*float(GP) for i,GP in enumerate(
        l_vcf[i_sample+9].split(':')[2].split(',')))

    return dosage


def parse_dosage_gen(l_gen,i_sample):

    dosage = 0
    dosage += float(l_gen[6+3*i_sample])
    dosage += 2*float(l_gen[7+3*i_sample])

    return dosage


def parse_dosage_bgl(l_bgl,i_sample):

    dosage = 0
    dosage += float(l_bgl[4+3*i_sample])
    dosage += 2*float(l_bgl[5+3*i_sample])

    return dosage


def parse_dosage_bed(bytesSNP,i_sample):

    d_dosage = {'00':0.,'01':1.,'11':2.}

    i_byte = int(i_sample/4)
    byte = bytesSNP[i_byte]
    s8 = str(bin(byte)[2:]).zfill(8)
    if len(s8) != 8:
        stop1
    i = 3-i_sample%4
    s2 = s8[2*i+1]+s8[2*i]
    if s2 == '10':
        dosage = None ## ask Deepti whether 10 should be 1 or N/A!!!
    else:
        dosage = d_dosage[s2]

    return dosage


def index_samples(files1,format1,files2,format2,fam1,fam2,fp_sampledic):

    d_samples = {}
    for i,files,fileformat,fam in [
        (1,files1,format1,fam1),(2,files2,format2,fam2)]:
        if fileformat == 'bed':
            l_samples = fam2samples('%s.fam' %(files[0]))
        elif fileformat == 'gen':
            with open(files[0]) as f:
                n_fam = int((len(f.readline().split())-5)/3)
            l_samples = fam2samples(fam)
        elif fileformat == 'bgl':
            with open(files[0]) as f:
                l_samples = f.readline().split()[3:-1:3]
        elif fileformat == 'vcf':
            l_samples = vcf2samples(files[0])
        else:
            print(fileformat,files)
            stop
        d_samples[i] = l_samples

    d_indexes = {1:[],2:[]}
    if fp_sampledic:
        with open(fp_sampledic) as lines:
            d_sampledic = {line.split()[0]:line.rstrip().split()[1] for line in lines}
    else:
        d_sampledic = {sample:sample for sample in d_samples[1]}
    bool_error = True
    for index1 in range(len(d_samples[1])):
        sample1 = d_samples[1][index1]
        try:
            index2 = d_samples[2].index(d_sampledic[sample1])
        except ValueError:
            continue
        d_indexes[1] += [index1]
        d_indexes[2] += [index2]
        bool_error = False

    if bool_error:
        print('0 samples intersect')
        raise('0 samples intersect')
        sys.exit()

    return d_samples, d_indexes


def vcf2samples(file):

    if file[-3:] == '.gz':
        f = gzip.open(file,'rt')
    else:
        f = open(files[0])
    for line in f:
        if line[0] != '#':
            break
        header = line
    f.close()
    l_samples = header.rstrip().split('\t')[9:]

    return l_samples


def fam2samples(fam):

    l_samples = []
    keyword = re.compile(r'(\d\d\d\d\d\d_[A-H]\d\d_)(.+\d\d\d\d\d\d\d)')
    with open(fam,'r') as lines:
        for line in lines:
            sampleID = line.split()[0]
            match = result = keyword.search(sampleID)
            if not match:
                continue
            l_samples += [match.group(2)]

    return l_samples


def parse_vcf(**kwargs):

    vcf = kwargs['file']
    line = vcf.readline()
    if line == '':
        return
    while line[0] == '#':
        line = vcf.readline()
    l = line.rstrip().split('\t')
    chrom = int(l[0])
    pos = int(l[1])
    alleleA = l[3]
    alleleB = l[4]

    yield chrom,pos,alleleA,alleleB,l


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
##        StopIteration
    pos = int(l_bim[3])
    alleleA = l_bim[4]
    alleleB = l_bim[5]

    ## read bed bytes corresponding to bim line
    bytesSNP = bed.read(n_bytes)

    yield chrom,pos,alleleA,alleleB,bytesSNP


def parse_extract(f_extract):

    line_extract = f_extract.readline()
    if line_extract == '':
        return
    l = line_extract.split()
    chrom = int(l[0])
    pos = int(l[1])

    yield chrom,pos


def parse_bgl(**kwargs):

    file = kwargs['file']
    fileinput = kwargs['fileinput']

    line = file.readline()
    if line == '':
        return ## StopIteration
    if fileinput.isfirstline():
        print(fileinput.filename())
        line = file.readline()
    l = line.rstrip().split()
    l0 = l[0].split(':')
    chrom = int(l0[0])
    pos = int(l0[1])
    alleleA = l[1]
    alleleB = l[2]

    yield chrom,pos,alleleA,alleleB,l


def parse_gen(**kwargs):

    file = kwargs['file']
    fileinput = kwargs['fileinput']

    line = file.readline()
    if line == '':
        return
    l = line.rstrip().split()
    pos = int(l[2])
    chrom = int(fileinput.filename().split('/')[-2]) ## -2 is a magic number... instead provide a file with chromosome-file-equivalents
    alleleA = l[3]
    alleleB = l[4]

    yield chrom,pos,alleleA,alleleB,l


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


def remove_zero_size_files(files):

    for file in list(files):
        if not os.path.isfile(file):
            continue
        if not os.path.getsize(file):
            files.remove(file)

    return files


def argparser():

    parser = argparse.ArgumentParser()

#    parser.add_argument('--gen',dest='gen',type=str,required = True,nargs='+')
#    parser.add_argument('--bfile','--bed',dest='bfile',required = True,metavar='FILE',type=str)
    parser.add_argument('--file1',required=True,nargs='+')
    parser.add_argument('--file2',required=True,nargs='+')
    parser.add_argument(
        '--format1',required=True,choices=['bed','gen','bgl','vcf'])
    parser.add_argument(
        '--format2',required=True,choices=['bed','gen','bgl','vcf'])
    ## type=argparse.FileType('rt')
    parser.add_argument('--extract')
    parser.add_argument('--fam1')
    parser.add_argument('--fam2')
    parser.add_argument('--affix','--out','--prefix','--suffix',)
    parser.add_argument('--sampledic',)

    namespace_args = parser.parse_args()

    d = vars(namespace_args)

    for k,v in d.items():
        if v == 'None':
            d[k] = None

    files1 = d['file1']
    files1 = sort_nicely(files1)

    files2 = d['file2']
    files2 = sort_nicely(files2)

    format1 = d['format1']
    format2 = d['format2']

    fp_extract = d['extract']

    fam1 = d['fam1']
    fam2 = d['fam2']

    affix = d['affix']

    fp_sampledic = d['sampledic']

    for fp in files1+files2+[fam1,fam2,fp_extract]:
        if not fp:
            continue
        if not os.path.isfile(fp) and not os.path.isfile('%s.bed' %(fp)):
            print('does not exist', fp)
            sys.exit()

    files1 = remove_zero_size_files(files1)
    files2 = remove_zero_size_files(files2)
    if len(files1) == 0 or len(files2) == 0:
        print('all files size zero')
        sys.exit()

    return files1,files2,format1,format2,fam1,fam2,fp_extract,affix,fp_sampledic


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
