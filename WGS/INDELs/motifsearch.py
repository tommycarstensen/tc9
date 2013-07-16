#!/bin/python3

## Tommy Carstensen (tc9)
## Wellcome Trust Sanger Institute, June-July 2013

import os
import itertools
import glob
import fileinput
import sys
import scipy
from scipy import stats
import argparse
import math
import re
import contextlib
import gzip

class motifsearch():


    def __init__(self,):
        
        self.d_comp = {'A':'T','C':'G','G':'C','T':'A','N':'N',}

        return


    def main(self):

        path_fa, path_fai = self.def_paths()

        self.parse_args()

        self.files_vcf = self.sort_nicely(self.files_vcf)

##        self.motif_reduction()
##        sys.exit()

        self.d_fai = self.read_fai(path_fai)

        d_cnt = self.combinations()

        d_cnt, d_cnt_all = self.loop_vcf(path_fa,d_cnt)

        self.stats(d_cnt,d_cnt_all)

        return


    def lookup_motifs(self,motif,motifsub,d_motifs,):

        l_motifs = d_motifs[motifsub]
        l_motifs.remove(motif)

        return l_motifs


    def lookup_motifs_revcomp(self,motif,motifsub,d_motifs,):

        try:
            l_motifs = [self.get_motif_revcomp(s) for s in d_motifs[motifsub]]
        except KeyError:
            l_motifs = []
            pass
        try:
            l_motifs.remove(motif)
        except ValueError:
            pass

        return l_motifs


    def motif_reduction(self,):

        fd = open('out_motifsearch/significant_1_100_12.txt')
        lines = fd.readlines()
        fd.close()

        l_tuples = []
        d_motifs2 = {}
        ## loop over lines sorted by pvalue in reverse (sort -k2gr,2)
        for line in lines:
            l = line.rstrip().split()
            oddsratio = float(l[0])
            pvalue = float(l[1])
            motif = l[4]
            l_tuples += [(pvalue,1/oddsratio,motif)]
            ## append to list
            try:
                d_motifs2[motif[1:]] += [motif]
                d_motifs2[motif[:-1]] += [motif]
            ## initiate list
            except KeyError:
                d_motifs2[motif[1:]] = [motif]
                d_motifs2[motif[:-1]] = [motif]
        del motif
        l_tuples.sort()
        l_motifs1 = [l[2] for l in l_tuples]
        n1 = len(l_motifs1)
        print('len',n1)
        l_motifs_nonoverlap = []
        maxlen = [self.length_motif,'NNN']
        while l_motifs1:
            motif1 = l_motifs1.pop()
            motif1_revcomp = self.get_motif_revcomp(motif1)
            try:
                l_motifs1.remove(motif1_revcomp)
            except ValueError:
                pass

##            if not motif1_revcomp[1:] in d_motifs2.keys(): continue ## tmp!!!
##            if len(self.lookup_motifs_revcomp(motif1_revcomp[1:],d_motifs2,)) < 3:
##                continue ## tmp!!!

            while True:

                l_motifs2_append = self.lookup_motifs(motif1,motif1[1:],d_motifs2,)
                l_motifs2_prepend = self.lookup_motifs(motif1,motif1[:-1],d_motifs2,)
                l_motifs2_prepend_revcomp = self.lookup_motifs_revcomp(
                    motif1,motif1_revcomp[1:],d_motifs2,)
                l_motifs2_append_revcomp = self.lookup_motifs_revcomp(
                    motif1,motif1_revcomp[:-1],d_motifs2,)

                l_prepends = []
                l_appends = []
                l_remove = []
                if not any([
                    l_motifs2_append,l_motifs2_prepend,
                    l_motifs2_prepend_revcomp,l_motifs2_append_revcomp,]):
                    break
                
                print(motif1)
                print('append',l_motifs2_append,l_motifs2_append_revcomp)
                print('prepend',l_motifs2_prepend,l_motifs2_prepend_revcomp)
                stop

            continue

            while True:

                l_prepends = []
                l_appends = []
                l_remove = []
                for motif2 in l_motifs2:
                    ##
                    ## 1 of 4
                    ##
                    if motif2[1:] in motif1:
                        if len(motif1) > self.length_motif:
                            index = motif0.index(motif[1:])
                            if len(motif1) == self.length_motif+1 and index == 0:
                                l_prepends += [motif[0]]
                                l_remove += [motif]
                            elif index != 1:
##                            if len(motif0) == self.length_motif+xxx and index == yyy:
                                print('index',index)
                                print(0,motif0[index:])
                                print(1,motif)
                                print(0,motif0)
                                stop
                                pass
                            else:
                                print(index)
                                print(motif0[index:])
                                print(motif)
                                print(motif0)
                                stop1
                                l_prepends += [motif[0]]
                                l_remove += [motif]
                        else:
##                            print(1,motif0)
                            l_prepends += [motif[0]]
                            l_remove += [motif]
                    ##
                    ## 2 of 4
                    ##
                    elif motif[:-1] in motif0:
                        l_remove += [motif]
                        if len(motif0) > self.length_motif:
                            index = motif0.index(motif[:-1])
                            if len(motif0) == self.length_motif+1 and index == 0:
                                pass
                            elif len(motif0) == self.length_motif+2 and index == 0:
                                pass
                            elif len(motif0)-index == self.length_motif-1:
                                l_appends += [motif[-1]]
                                pass
                            else:
                                print('index',index)
                                print(' '*index+motif0[index:])
                                print(' '*index+motif)
                                print(motif0)
                                print(len(motif0)-index)
                                stop
                                if index == 3:
                                    print('index',index)
                                    print(motif0[index:])
                                    print(motif)
                                    print(motif0)
                                    print(index,len(motif0))
                                    stop2
                                l_appends += [motif[-1]]
                                l_remove += [motif]
    ##                    motif0 = motif0+motif[-1]
##                        print(2,motif0,motif)
                        else:
                            l_appends += [motif[-1]]
                            l_remove += [motif]
                    elif motif == motif0_revcomp:
                        l_remove += [motif]
                    ##
                    ## 3 of 4
                    ##
                    elif motif[1:] in motif0_revcomp:
                        index = motif0_revcomp.index(motif[1:])
                        if index == 1: pass
                        else:
                            print(motif0,self.get_motif_revcomp(motif))
                            stoptmp
                            l_remove += [motif]
                            if len(motif0) > self.length_motif:
                                print(motif)
                                print(motif0)
                                stop3
                            l_appends += [self.get_motif_revcomp(motif)[-1]]
                    ##
                    ## 4 of 4
                    ##
                    elif motif[:-1] in motif0_revcomp:
                        index = motif0_revcomp.index(motif[:-1])
                        if index == 0:
                            pass
                        elif len(motif0) == self.length_motif and index == 1:
                            print(motif,self.get_motif_revcomp(motif))
                            stoptmp
                            l_remove += [motif]
                            l_prepends += [self.get_motif_revcomp(motif)[0]]
                        elif len(motif0) > self.length_motif:
                            print(motif,'motif')
                            print(motif0,'motif0')
                            print(motif0_revcomp,'motif0_revcomp')
                            print(index)
                            print(self.get_motif_revcomp(motif))
                            stop4
    ##                    motif0 = self.get_motif_revcomp(motif)[0]+motif0
##                        print(4,motif0)
                        l_prepends += [self.get_motif_revcomp(motif)[0]]
                    else:
                        ## continue loop over motifs possibly matching motif0
                        continue

                ##
                ## no matches at either end
                ##
                if len(l_appends) == 0 and len(l_prepends) == 0:
                    l_motifs_nonoverlap += [motif0]
                    if len(motif0) > maxlen[0]:
                        maxlen = [len(motif0),motif0]
                    break
                ##
                ## multiple matches at both ends
                ##
                elif len(l_appends) > 1 and len(l_prepends) > 1:
                    break
                ##
                ## single match at one or both ends
                ##
                else:
##                    print('bef',motif0)
                    if len(l_appends) == 1:
                        motif0 = motif0+''.join(l_appends)
                    if len(l_prepends) == 1:
                        motif0 = ''.join(l_prepends)+motif0
##                    print('aft',motif0)
##                    print()
                    for motif in l_remove:
                        if not (motif[1:-1] in motif0 or motif[1:-1] in motif0_revcomp):
                            print(motif)
                            print(motif0)
                            print(motif0_revcomp)
                            stoptmp
                        self.remove(l_motifs,motif)
                        self.remove(l_motifs,self.get_motif_revcomp(motif))
                    continue

            ## continue pop
            n = len(l_motifs)
            if n % 10 == 0:
                print(n,'maxlen',maxlen)
            continue

        n2 = len(l_motifs_nonoverlap)
        for i in range(20):
            print(l_motifs_nonoverlap[i])
        print(n1,n2)
        print(maxlen)
        sys.exit()

        return


    def remove(self,l,x):

        try:
            l.remove(x)
        except ValueError:
            pass

        return


    def get_motif_revcomp(self,motif):

        motif_revcomp = ''
        for i in range(1,self.length_motif+1,):
            motif_revcomp += self.d_comp[motif[-i]]

        return motif_revcomp


    def stats(self,d_cnt,d_cnt_all,):

        ## Bonferroni correction
        pvalue_max = 5*(10**-2)/(4**self.length_motif)
        print('pvalue_max after Bonferroni correction',pvalue_max)

        print('do Fisher exact test on each motif')

        l = []
        fd_notobserved = open(
            'out_motifsearch/not_observed_%s_%s.txt' %(self.size_window,self.length_motif),'w',1)
        fd_infrequent = open(
            'out_motifsearch/infrequent_%s_%s.txt' %(self.size_window,self.length_motif),'w',1)
        fd1 = open(
            'out_motifsearch/significant_1_%s_%s.txt' %(self.size_window,self.length_motif),'w',1)
        fd2 = open(
            'out_motifsearch/significant_2_%s_%s.txt' %(self.size_window,self.length_motif),'w',1)
        i_motif = 0
        for motif in d_cnt['1'].keys():

            if i_motif%10000 == 0:
                print(i_motif,motif) ## tmp!!!
                sys.stdout.flush()
            ## append count
            i_motif += 1

            a = cnt_motif1 = d_cnt['1'][motif][0]
            c = cnt_motif2 = d_cnt['2'][motif][0]
            b = d_cnt_all['1']-cnt_motif1
            d = d_cnt_all['2']-cnt_motif2
            ## count the reverse complementary motif as well
            if self.bool_revcomp == True:
                motif_revcomp = self.get_motif_revcomp(motif)
                cnt_motif_revcomp1 = d_cnt['1'][motif_revcomp][0]
                cnt_motif_revcomp2 = d_cnt['2'][motif_revcomp][0]
                a += cnt_motif_revcomp1
                c += cnt_motif_revcomp2
                b -= cnt_motif_revcomp1
                d -= cnt_motif_revcomp2

            if cnt_motif1 == 0 and cnt_motif2 == 0:
                fd_notobserved.write('%s\n' %(motif))
                continue
            elif cnt_motif1 < 1000 and cnt_motif2 < 1000:
                fd_infrequent.write('%s\n' %(motif))
                continue

            if cnt_motif1 == 0:
                meandist1 = 'N/A'
            else:
                meandist1 = '%.1f' %(d_cnt['1'][motif][1]/d_cnt['1'][motif][2])
            if cnt_motif2 == 0:
                meandist2 = 'N/A'
            else:
                meandist2 = '%.1f' %(d_cnt['2'][motif][1]/d_cnt['2'][motif][2])

            obs = l_contigency = l_contigency_table = [[a,b],[c,d],]

            if cnt_motif1 > d_cnt_all['1']:
                print(motif,d_cnt['1'][motif],d_cnt_all)
                stop1tmp
            if cnt_motif2 > d_cnt_all['2']:
                print(motif,d_cnt['2'][motif],d_cnt_all)
                stop2tmp

##            oddsratio,pvalue = stats.fisher_exact(l_contigency_table)

            if c == 0:
                fd = fd1
                pvalue = 0
                oddsratio,pvalue = stats.fisher_exact(l_contigency_table)
            elif a == 0:
                fd = fd2
                oddsratio,pvalue = stats.fisher_exact(l_contigency_table)
            else:
                oddsratio = a*d/(b*c)

                ## standard error ## http://en.wikipedia.org/wiki/Odds_ratio#Statistical_inference
                SE = math.sqrt(1/a+1/b+1/c+1/d)
                Zscore = math.log(oddsratio)/SE
                pvalue = scipy.stats.zprob(Zscore)
##                p2 = scipy.special.ndtr(Zscore)

                if oddsratio > 1:
                    pvalue = 2*(1-pvalue)
                    fd = fd1
                else:
                    pvalue = 2*pvalue
                    fd = fd2

            if pvalue <= pvalue_max:
##                print(oddsratio,pvalue,cnt_motif1,cnt_motif2,motif,)
                fd.write('%s %s %i %i %s %s %s\n' %(
                    oddsratio,pvalue,cnt_motif1,cnt_motif2,motif,meandist1,meandist2))

##            print('time')
##            import time
##
##            t1 = time.time()
##            for x in range(1000):
##                l_contigency_table = [[a,b],[c,d],]
##                oddsratio,pvalue = stats.fisher_exact(l_contigency_table)
##            t2 = time.time()
##            print(t2-t1)
##
##            t1 = time.time()
##            for x in range(1000):
##                SE = math.sqrt(1/a+1/b+1/c+1/d)
##                Zscore = math.log(oddsratio)/SE
##                p1 = scipy.stats.zprob(Zscore)
##            t2 = time.time()
##            print(t2-t1)
##
##            t1 = time.time()
##            for x in range(1000):
##                chi2, p, dof, ex = scipy.stats.chi2_contingency(l_contigency_table)
##            t2 = time.time()
##            print(t2-t1)
##
##            print(pvalue,p1,pchi)
##
##            stop

##            print(oddsratio,pvalue,'SE',SE,'Z',Zscore,2*p1,2*p2,2*p1-pvalue)
##            odssratio = (
##                l_contigency_table[1][1]*l_contigency_table[0][0])/
##                (l_contigency_table[1][0]*l_contigency_table[0][1])
##                )
##            print(oddsratio,pvalue,cnt_motif1,cnt_motif2,motif,)
##            print(d_cnt_all['1']-cnt_motif1,d_cnt_all['2']-cnt_motif2,)
##            stop

        fd1.close()
        fd2.close()
        fd_notobserved.close()
        fd_infrequent.close()

        return


    def parse_rm(self,file_rm):

        line = file_rm.readline().decode("utf-8")
        l_rm = line.rstrip().split()
        chrom_rm = int(l_rm[4][3:])
        pos1_rm = int(l_rm[5])
        pos2_rm = int(l_rm[6])

        yield chrom_rm, pos1_rm, pos2_rm 


    def loop_vcf(self,path_fa,d_cnt,):

        print('loop over vcf(s)')

        d_cnt_all = {'1':0,'2':0,}

        self.cnt_newlines_minimum = self.size_window//60
        self.window_line_mod = self.size_window%60
        self.cnt_motifs_per_window = self.size_window-self.length_motif+1

        input_filenames = []
        output_filenames = ['a.txt','b.txt','c.txt',]
        with contextlib.ExitStack() as stack:
            ## FASTA sequence
            self.file_fa = stack.enter_context(open(path_fa))
            ## VCF file(s) with INDEL positions
            file_vcf = stack.enter_context(fileinput.input(files=self.files_vcf))
            ## FASTA output files to write upstream sequences to
            file_odd = stack.enter_context(open('odd%s.fa' %(self.size_window), 'w'))
            file_even = stack.enter_context(open('even%s.fa' %(self.size_window), 'w'))
            ## RepeatMasker
##            self.file_zip_rm = stack.enter_context(zipfile.ZipFile('hg19.fa.out.gz'))
            file_rm = stack.enter_context(gzip.open('hg19.fa.out.gz'))
##            inputs = [stack.enter_context(fileinput.input(files=filename)) for filename in input_filenames]
##            d_outputs = {filename:stack.enter_context(open(filename, "w")) for filename in output_filenames}
##        with open(path_fa) as self.file_fa, fileinput.input(self.files_vcf) as file_vcf:
##          with open('odd%s.fa' %(self.size_window),'w') as file_odd, open('even%s.fa' %(self.size_window),'w') as file_even:
            for i in range(3):
                file_rm.readline().decode("utf-8")
            chrom_rm, pos1_rm, pos2_rm = next(self.parse_rm(file_rm))
##            chrom,pos,l_vcf = next(self.generate_line_vcf_INDEL(file_vcf))
##            while True:
##            stop
            for chrom,pos,l_vcf in self.generate_line_vcf_INDEL(file_vcf):
                self.parse_line(
                    chrom,pos,l_vcf,d_cnt,d_cnt_all,file_odd,file_even,)

        return d_cnt, d_cnt_all


    def parse_line(
        self,chrom,pos,l_vcf,d_cnt,d_cnt_all,file_odd,file_even,):

        byte_init = self.d_fai[chrom]['start']
##                quotient = floor = (pos-1)//60
##                remainder = modulus = (pos-1)%60
        quotient,remainder = divmod(pos-1,60)
        cnt_newlines = self.cnt_newlines_minimum
        if remainder <= self.window_line_mod:
            cnt_newlines += 1
##        self.cnt_newlines_minimum = upstream//60
##        self.window_line_mod = upstream%60
##        self.cnt_motifs_per_window = upstream-self.length_motif+1
##        cnt_newlines = self.cnt_newlines_minimum
        self.file_fa.seek(
            byte_init+pos-1+quotient-self.size_window-cnt_newlines)
        seq = self.file_fa.read(self.size_window+cnt_newlines).replace('\n','')

##        ## check that parsed sequence has the size of the window
##        if len(seq) != self.size_window:
##            print(len(seq),size_window)
##            stop

##            ## check for ref seq and vcf discrepancies
##            file_fa.seek(byte_init+pos-1+quotient)
##            s = file_fa.read(80).replace('\n','')
##            if s[:len(l_vcf[3])] != l_vcf[3]:
##                print(s)
##                print(s[:len(l_vcf[3])])
##                print(l_vcf[3])
##                stop

##        ## check that ref allele is mono-allelic
##        if ',' in l_vcf[3]:
##            print(l_vcf[3],l_vcf[4])
##            stop

        if self.bool_weight_by_sample == True:
            d_cnt_samples = {}
            for i in range(9,len(l_vcf)):
                gt = l_vcf[i].split(':')[0]
                if gt == './.':
                    continue
                for s in gt.split('/'):
                    try:
                        d_cnt_samples[int(s)] += 1
                    except KeyError:
                        d_cnt_samples[int(s)] = 1

        bool1 = False
        bool2 = False
        for alleleA in l_vcf[3].split(','):
            lenA = len(alleleA)
            l_alleleB = l_vcf[4].split(',')
##                for alleleB in l_alleleB:
            for i_alleleB in range(len(l_alleleB)):
                alleleB = l_alleleB[i_alleleB]
                if self.bool_weight_by_sample == True:
                    try:
                        cnt_samples = d_cnt_samples[i_alleleB+1]
                    except KeyError:
                        cnt_samples = 0
                else:
                    cnt_samples = 1
                lenB = len(alleleB)
                lenINDEL = abs(lenB-lenA)
                if lenINDEL < 2: continue
##                if lenINDEL > 0:
##                    continue
##                if lenINDEL%12 == 0:
##                    continue
##                elif lenINDEL%3 == 0:
                elif lenINDEL%2 == 0: ## even
                    d_cnt_all['1'] += self.append_total_count(cnt_samples)
                    self.append_motif(d_cnt['1'],seq,cnt_samples,)
                    bool1 = True
                    for x in range(cnt_samples):
                        file_even.write('>%s|%s|%s|%s\n%s\n' %(
                            chrom,pos,l_vcf[3],l_vcf[4],seq,))
##                elif lenINDEL%4 == 0:
                else: ## odd
                    d_cnt_all['2'] += self.append_total_count(cnt_samples)
                    self.append_motif(d_cnt['2'],seq,cnt_samples,)
                    bool2 = True
                    for x in range(cnt_samples):
                        file_odd.write('>%s|%s|%s|%s\n%s\n' %(
                            chrom,pos,l_vcf[3],l_vcf[4],seq,))
##            if bool3n == True:
##                file_FASTA3.write('>%s|%s|%s|%s\n%s\n' %(chrom,pos,l_vcf[3],l_vcf[4],seq,))
##            if bool4n == True:
##                file_FASTA4.write('>%s|%s|%s|%s\n%s\n' %(chrom,pos,l_vcf[3],l_vcf[4],seq,))

##        if bool1 == True or bool2 == True:
##            if bool1 == True and bool2 == True:
##                v = 3
##            elif bool1 == True:
##                v = 1
##            else:
##                v = 2
##            try:
##                motif = 'GACACATGCACA'
##                motif = 'GCATATATACGT'
##                index = seq.index(motif)
##                with open('%s.graph%i' %(motif,self.size_window),'a') as f:
##                    f.write('chr%s %i %i %i\n' %(chrom,pos,v,cnt_samples))
##            except ValueError:
##                pass

        return


    def append_total_count(self,cnt_samples,):

        cnt = self.cnt_motifs_per_window
        if self.bool_weight_by_sample == True:
            cnt *= cnt_samples

        return cnt


    def combinations(self,):

        print('create combinations for hash table')
        l_nt = ['A','C','G','T',]
        g_motifs = itertools.product(l_nt,repeat=self.length_motif)
        d_cnt = {'1':{},'2':{},}
        for t_motif in g_motifs:
            s_motif = ''.join(t_motif)
            d_cnt['1'][s_motif] = [0,0,0]
            d_cnt['2'][s_motif] = [0,0,0]

        return d_cnt


    def append_motif(self,d_count_motif,seq,cnt_samples,):

        dist = len(seq)
        for motif in self.generate_motif(seq):
            try:
                cnt = 1
                if self.bool_weight_by_sample == True:
                    cnt *= cnt_samples
                d_count_motif[motif][0] += cnt
                d_count_motif[motif][1] += dist
                d_count_motif[motif][2] += 1
            except KeyError:
                if not 'N' in motif:
                    print(motif,'N not in motif')
                    sys.exit()
            dist -= 1

        return


    def generate_motif(self,seq):

        for i in range(len(seq)-self.length_motif+1):
            motif = seq[i:i+self.length_motif]

            yield motif


    def generate_line_vcf_INDEL(self,f):

        set_nt = set(['A','C','G','T',])

        for line in f:
            if line[0] == '#': continue
            l = line.rstrip().split('\t')
            ## skip SNPs
            if (
                l[3] in ['A','C','G','T',]
                and
                l[4] in [
                    ## UnifiedGenotyper
                    'A','C','G','T',
                    'A,C','A,G','A,T','C,G','C,T','G,T',
                    'A,C,G','A,C,T','A,G,T','C,G,T',
    ##                ## CombineVariants
    ##                'C,A','G,A','T,A','G,C','T,C','T,G',
    ##                'A,T,C','T,C,A','C,T,G','A,T,G','T,A,G','C,T,A','T,G,C','G,C,T',
    ##                ## UG GENOTYPE_GIVEN_ALLELES
    ##                'C,G,A','C,A,T',
    ##                'G,A,C','G,T,A','G,T,C','G,C,A',
    ##                'A,G,C',
    ##                'T,A,C','T,G,A',
                    ]
                ):
                continue
            chrom = l[0]
            pos = int(l[1])

            yield chrom,pos,l


    def read_fai(self,path_fai):

        print('parsing fai')

        d = {}
        with open(path_fai) as file_fai:
            for line_fai in file_fai:
                l = line_fai.rstrip().split()
                chrom = l[0]
                byte_length = int(l[1])
                byte_start = int(l[2])
                bytes_per_line_excl_line_break = int(l[3])
                bytes_per_line_incl_line_break = int(l[4])
                d[chrom] = {'length':byte_length,'start':byte_start}
    ##            print(line_fai.rstrip())
                if chrom == 'Y':
                    break

        return d


    def def_paths(self,):

        path1 = '/lustre/scratch111/resources'
        path2 = 'ref/Homo_sapiens/1000Genomes_hs37d5'
        path_fa = os.path.join(path1,path2,'hs37d5.fa')
        path_fai = os.path.join(path1,path2,'hs37d5.fa.fai')

        return path_fa, path_fai


    def parse_args(self,):

        parser = argparse.ArgumentParser()

        parser.add_argument(
            '--size_window','--window',
            dest='size_window',
            help='size (bp) of window upstream of INDEL',
            type=int,metavar='int',default=100,
            required = True,
            )

        parser.add_argument(
            '--length_motif','--motif',
            dest='length_motif',
            help='size (bp) of motif',
            type=int,metavar='int',default=10,
            required = True,
            )

        parser.add_argument(
            '--weight','--weight_by_sample','--bool_weight_by_sample',
            dest='bool_weight_by_sample',
            help="for Fisher's exact test weight each motif by the number of samples in which it is present",
            action='store_true',
            required = False,
            )

        parser.add_argument(
            '--revcomp',
            dest='bool_revcomp',
            help="do reverse complementary strand",
            action='store_true',
            required = False,
            )

        parser.add_argument(
            '--vcfs','--vcf',
            dest='files_vcf',
            help="vcf(s) from which to parse INDEL positions and lengths",
            nargs='+',
            required = True,
            )

        ## http://docs.python.org/2/library/argparse.html#argparse.ArgumentParser.parse_args
        ## parse arguments to argparse NameSpace
        self.namespace_args = namespace_args = parser.parse_args()

        ## http://docs.python.org/3/library/functions.html#vars
        for k,v in vars(namespace_args).items():
            setattr(self,k,v)

        if self.size_window < self.length_motif:
            print('window size cannot be smaller than motif length')
            sys.exit()

        return


    def alphanum_key(self,s):
        ## http://dave.st.germa.in/blog/2007/12/11/exception-handling-slow/
        NUM_RE = re.compile('([0-9]+)')
        return [ int(c) if c.isdigit() else c for c in NUM_RE.split(s) ]


    def sort_nicely(self,l):
        ## http://nedbatchelder.com/blog/200712/human_sorting.html
        """ Sort the given list in the way that humans expect.
        """
        l.sort(key=self.alphanum_key)
        return l


    def junk(self,):

        if sys.argv[-1] == 'WS':
            ## Weak/Strong
            d_nan = {'C':'S','G':'S','A':'W','T':'W',}
        elif sys.argv[-1] == 'MK':
            ## aMino/Keto
            d_nan = {'A':'M','C':'M','G':'K','T':'K',}
        elif sys.argv[-1] == 'RY':
            ## puRine/pYrimidine
            d_nan = {'A':'R','C':'Y','G':'R','T':'Y',}
        else:
            print(sys.argv[-1])
            sys.exit()

        return


if __name__ == '__main__':
    instance = motifsearch()
    instance.main()
