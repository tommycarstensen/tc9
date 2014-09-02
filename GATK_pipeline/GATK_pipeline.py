#!/bin/python3

## Tommy Carstensen, Wellcome Trust Sanger Institute, 2012-2014

import math
import os
import sys
import time
import re
import pwd
import argparse
import inspect
import glob
import fileinput
import itertools
import gzip
import subprocess
import contextlib


## README

## http://www.broadinstitute.org/gatk/about#typical-workflows
## http://www.broadinstitute.org/gatk/guide/topic?name=best-practices

class main():

    def main(self):

        ## parse chromsome lengths from reference sequence
        d_chrom_lens = self.parse_chrom_lens()

        self.HaplotypeCaller(d_chrom_lens,)
        self.CombineGVCFs()
        self.GenotypeGVCFs()

        ## 1000G_phase1.snps.high_confidence.b37.vcf.gz only contains chromosomes 1-22 and X
        self.VariantRecalibrator(d_chrom_lens,)

        ## phased imputation reference panels not available for chrY
        self.chroms.remove('Y')
        ## Exception in thread "main" java.lang.ArrayIndexOutOfBoundsException: -2
        self.chroms.remove('X')
        self.BEAGLE(d_chrom_lens)

        self.unite_BEAGLE(d_chrom_lens,)

        return


    def execmd(self,cmd):

        print(cmd)
        subprocess.call(cmd,shell=True)

        return


    def check_out_line(self,fp):

        bool_exit = True
        for i in [14,19,]:
            line = os.popen("tail -n%i %s | head -n1" %(i,fp)).read().strip()
            if line in [
                "Successfully completed.",
##                "Exited with exit code 127.",
##                "Exited with exit code 1.",
                ]: ## tmp!!!
                bool_exit = False
                break
        if bool_exit == True:
            print('not finished:', fp)
            sys.exit()

        return


    def BEAGLE_fileIO_checks(self,chrom):

        ##
        ## check that gprobs is identical to markers
        ##
        cmd = "awk '{print $1}' in_BEAGLE/%s/%s.*.markers" %(chrom,chrom)
        cmd += ' | sort -u > markers%s' %(chrom)
        self.execmd(cmd)

        if os.path.isfile('gprobs%s' %(chrom)):
            os.remove('gprobs%s' %(chrom))
        for file in glob.glob(
            'out_BEAGLE/%s/%s.*.like.gprobs.gz' %(chrom,chrom)):
            cmd = "zcat %s" %(file)
            cmd += " | awk 'FNR>1{print $1}'"
            cmd += ' >> gprobs%s' %(chrom)
            self.execmd(cmd)

        self.execmd('sort -u gprobs%s -o gprobs%s' %(chrom,chrom))

        cmd = 'comm -3 gprobs%s markers%s' %(chrom,chrom)
        i = int(os.popen('%s | wc -l' %(cmd)).read())
        if i > 0:
            print(os.popen('%s | head' %(cmd)).readlines())
            print(cmd)
            sys.exit()

        os.remove('markers%s' %(chrom))
        os.remove('gprobs%s' %(chrom))

        return


    def BEAGLE_unite(self,chrom):

        ##
        ##
        ##
        fd = open('BEAGLE_divide_indexes.txt','r')
        lines = fd.readlines()
        fd.close()
        d_index2pos = {}
        for line in lines:
            l = line.strip().split()
            if l[0] != chrom:
                continue
            index = int(l[1])
            pos1 = int(l[2])
            pos2 = int(l[3])
            d_index2pos[index] = [pos1,pos2,]
        fp_out = 'out_BEAGLE/%s.gprobs.gz' %(chrom,)
        for index in range(1,max(d_index2pos.keys())+1,):
            fp_in = 'out_BEAGLE/%s/%s.%i.like.gprobs.gz' %(chrom,chrom,index,)
            if index == 1:
                if not os.path.isfile(fp_out):
                    cmd = 'zcat %s | head -n1 | gzip > %s' %(fp_in,fp_out,)
                    self.execmd(cmd)
                else:
                    break
            if index % 10 == 0:
                print('BEAGLE_unite', chrom, index)
            if not os.path.isfile(fp_in):
                print('missing', fp_in)
                stop
            pos1 = d_index2pos[index][0]
            pos2 = d_index2pos[index][1]
            cmd = 'zcat %s' %(fp_in)
            cmd += " | awk 'NR>1{pos=int(substr($1,%i));" %(len(chrom)+2)
            cmd += " if(pos>%i&&pos<=%i) print $0}'" %(pos1,pos2)
            cmd += ' | gzip >> %s' %(fp_out,)
            self.execmd(cmd)

        return


    def unite_BEAGLE(self, d_chrom_lens,):

        ##
        ## 1) check input existence
        ##
        l_fp_in = self.BEAGLE_parse_output_files()
        if self.check_in('BEAGLE', l_fp_in, 'touch/BEAGLE.touch'):
            sys.exit()

        for chrom in self.chroms:
            print(chrom)
            fp = 'out_BEAGLE/%s.gprobs.gz' %(chrom)
            if os.path.isfile(fp):
                continue
            self.BEAGLE_fileIO_checks(chrom)
            self.BEAGLE_unite(chrom)

        return


    def BEAGLE_parse_output_files(self,):

        ## open file
        fd = open('BEAGLE_divide_indexes.txt','r')
        ## read all lines into memory
        lines = fd.readlines()
        ## close file
        fd.close()

        l_fp_in = []
        for line in lines:
            l = line.strip().split()
            chrom = l[0]
            if not chrom in self.chroms: continue
            index = int(l[1])
            fp_in = 'out_BEAGLE/%s/' %(chrom,)
            fp_in += '%s.%i.like.gprobs.gz' %(chrom,index,)
            l_fp_in += [fp_in]

        return l_fp_in


    def write_shell(self,fp,lines,):

        self.mkdir(fp)

        if type(lines) != list:
            print(type(lines))
            stop

        s = '\n'.join(lines)+'\n\n'
        fd = open(fp,'w')
        fd.write(s)
        fd.close()
        os.system('chmod +x %s' %(fp))

        return


    def parse_marker(self,line_m,):

        l_markers = line_m.split()
        pos_ref = int(l_markers[1])
        A_ref = l_markers[2]
        B_ref = l_markers[3]

        return pos_ref, A_ref, B_ref


    def BEAGLE_divide(self,):

        d_indexes = {}
        for chrom in self.chroms:
            d_indexes[chrom] = {}

        ## parse the input tranches file describing where to cut the data
        fp_tranches = 'out_VariantRecalibrator/VariantRecalibrator.SNP.tranches'
        minVQSLOD = self.parse_minVQSLOD(fp_tranches,self.ts_filter_level,)

        ## open the recal file
        fp_recal = 'out_VariantRecalibrator/VariantRecalibrator.SNP.recal'

        for vcf in glob.glob('out_UnifiedGenotyper/*.vcf.gz'):
            assert os.path.getmtime(vcf) < os.path.getmtime(fp_recal)

        with open(fp_recal,'r') as fd_recal:

            ## loop over the raw input variants to be recalibrated
            if os.path.islink('out_UnifiedGenotyper'):
                path = os.readlink('out_UnifiedGenotyper')
            else:
                path = 'out_UnifiedGenotyper'

            for chrom in self.chroms:
                if int(chrom) != 16: continue ## tmp!!!
                print('chrom %s' %(chrom), flush=True)
                s = '%s/%s.*.vcf.gz' %(path,chrom,)
                l_files = glob.glob(s)
                if len(l_files) == 0:
                    print('no files found')
                    print(s)
                    sys.exit()
                l_files_sorted = self.sort_nicely(l_files)
                d_indexes[chrom] = self.loop_UG_out(
                    chrom,l_files_sorted,fd_recal,minVQSLOD)

        return d_indexes


    def parse_minVQSLOD(self,fp_tranches,ts_filter_level,):

        fd = open(fp_tranches)
        lines = fd.readlines()
        fd.close()
        for line in lines:
            if line[0] == '#':
                continue
            l = line.split(',')
            if l[0] == 'targetTruthSensitivity':
                index = l.index('minVQSLod')
                continue
            targetTruthSensitivity = float(l[0])
            if targetTruthSensitivity == ts_filter_level:
                minVQSLOD = float(l[index])
                break

        return minVQSLOD


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


    def generate_line_vcf_PASS_split(self,fd_vcf,fd_recal,minVQSLOD):

        for line_vcf in fd_vcf:

            ## skip header
            if line_vcf[0] == '#':
                continue

            ## skip INDELs
            l_vcf = line_vcf.rstrip().split('\t')
            bool_indel,bool_diallelic = self.parse_variant_type(l_vcf)
            if bool_indel == True:
                continue

            ## skip LowQual
            if l_vcf[6] == 'LowQual':
                continue

##            ## skip QUAL=.
##            if l_vcf[5] == '.': continue

            if not 'line_recal' in locals(): ## tmp!!!
                ## read equivalent VR line
                while True:
                    line_recal = fd_recal.readline()
                    if line_recal[0] != '#':
                        break

            ## skip non-diallelic SNPs
            if bool_diallelic == False:
                continue

            ## Parse VQSLOD value from recal file.
            VQSLOD, l_recal = self.parse_VQSLOD(line_recal)

            ## tmp!!! temporary!!! TEMPORARY!!! temp!!!!!!!!!
            try:
                x = int(l_recal[0])
            except:
                return
            if int(l_vcf[0]) < int(l_recal[0]):
                print('a', l_recal[:2],l_vcf[:2])
                continue
            while int(l_recal[0]) < int(l_vcf[0]):
                line_recal = fd_recal.readline()
                VQSLOD, l_recal = self.parse_VQSLOD(line_recal)
            if int(l_vcf[1]) < int(l_recal[1]):
                print('c', l_recal[:2],l_vcf[:2])
                continue
            while int(l_recal[1]) < int(l_vcf[1]):
                print('d', l_recal[:2],l_vcf[:2])
                line_recal = fd_recal.readline()
                VQSLOD, l_recal = self.parse_VQSLOD(line_recal)
##            if int(l_vcf[1]) < int(l_recal[1]): ## 16.5 tmp!!!
##                print('e', l_recal[:2],l_vcf[:2])
##                continue

            ## Assert that position in recal and vcf file are identical.
            try:
                assert int(l_recal[1]) == int(l_vcf[1])
            except:
                print(line_recal)
                print(line_vcf)
                stop
            ## Skip if VQSLOD value below threshold.
            if VQSLOD < minVQSLOD:
                continue

            yield l_vcf


    def parse_VQSLOD(self,line_recal,):

        l_recal = line_recal.split('\t')
##        l_recal_INFO = l_recal[7].split(';')
##        VQSLOD = float(l_recal_INFO[1].split('=')[1])
        VQSLOD = float(re.search(r'VQSLOD=(\-?\d+.\d+)',l_recal[7]).group(1))

        return VQSLOD, l_recal


    def vcf2beagle_header(self,file):

        with gzip.open(file,'rt') as fd_vcf:
            for line_vcf in fd_vcf:
                ## skip header
                if line_vcf[0] != '#':
                    break
                line_vcf_header = line_vcf
        l_samples = line_vcf_header.strip().split('\t')[9:]
        l = ['%s %s %s' %(sample,sample,sample) for sample in l_samples]
        s = ' '.join(l)
        header = 'marker alleleA alleleB '+s+'\n'

        return header


    def chrom2file(self,file1,chrom):

        with open(file1) as f:
            for line in f:
                l = line.rstrip().split()
                if l[0] == str(chrom):
                    file2 = l[1]

        return file2


    def loop_UG_out(self,chrom,l_files,fd_recal,minVQSLOD):

        ## BEAGLE manual reads:
        ## "When imputing ungenotyped markers using a reference panel,
        ## the markers file should contain only the markers
        ## genotyped in the reference panel"
        ## Procedure would then be:
        ## 1) Refinement without reference panel.
        ## 2) Imputation of reference panel markers only.
        ## 3) updategprobs.jar to merge

        for file in l_files:
            if os.path.getsize(file) > 0:
                break
        header = self.vcf2beagle_header(file)
        d_index2pos = {}

        ##
        ## size and edge
        ##
        size = self.i_BEAGLE_size*1000000
        edge = self.i_BEAGLE_edge*1000

        ## genotype likelihood file
        fp_out_prefix = 'in_BEAGLE/%s/%s' %(chrom,chrom,)
        fp_markers = self.chrom2file(self.BEAGLE_markers, chrom)
        fp_phased = self.chrom2file(self.BEAGLE_phased, chrom)

        ##
        ## open files
        ##
        fd_phased = open(fp_phased,'r')
        fd_markers = open(fp_markers,'r')
        ##
        ## parse phased header and skip phased header
        ##
        header_phased = fd_phased.readline()

        ##
        ## initiate lines out
        ##
        lines_out_phased1 = [header_phased]
##        lines_out_phased2 = [header_phased]
        lines_out_markers1 = []
        lines_out_markers2 = []

        if os.path.isfile('%s.phased' %(fp_out_prefix,)):
            os.remove('%s.phased' %(fp_out_prefix,))
        fd = open('%s.phased' %(fp_out_prefix,),'w')
        fd.close()

        ##
        ## set smallest fragment size
        ##
        min_bps = 400000

        lines_out1 = [header]
        lines_out2 = [header]
        position = 0
        pos_init1 = position-position%size
        pos_term1 = self.calc_pos_term(pos_init1,position,size,min_bps)

        index = 1
        bool_append_markphas = False
        bool_append_to_prev = False
        bool_append_to_next = False
        pos_prev = 0
        pos_prev_EOF = None
        bool_BOF2 = False
        bool_EOF1 = False

        ##
        ## parse 1st markers line
        ##
        line_m = fd_markers.readline()
        pos_ref, A_ref, B_ref = self.parse_marker(line_m)
        line_p = fd_phased.readline()

        cnt_variants = 0
        with fileinput.input(
            files=l_files,openhook=self.hook_compressed_text) as fd_vcf:
            while True:
                try:
                    l_vcf = next(self.generate_line_vcf_PASS_split(
                        fd_vcf,fd_recal,minVQSLOD,))
                    cnt_variants += 1
                except StopIteration:
                    break

                position = int(l_vcf[1])
                alleleA = l_vcf[3]
                alleleB = l_vcf[4]

                ##
                ## avoid multiple comparisons of large integers
                ## by means of booleans instead of nesting
                ## do this immediately after getting current position
                ##
                if bool_BOF2 == False and position > pos_term1-edge:
                    bool_BOF2 = True
                    pos_init2 = pos_term1
                    ## centromere
                    if position > pos_term1+edge:
                        bool_EOF1 = True
                        if pos_prev-pos_init1 < min_bps:
                            bool_append_to_prev = True
                elif bool_BOF2 == True and bool_EOF1 == False and position > pos_term1+edge:
                    bool_EOF1 = True
                else:
                    pass


                ##
                ## loop to determine
                ## whether markers and phased should be appended or not
                ##
                while True:
                    ##
                    ## 1) same position
                    ##
                    if pos_ref == position:
                        ##
                        ## 1a) markers identical
                        ##
                        if (
                            A_ref == alleleA
                            and
                            B_ref == alleleB
                            ):
                            bool_append_markphas = True
                            pass
                        ##
                        ## 1b) markers different
                        ##
                        else:
                            bool_append_markphas = False
                            pos_ref_prev = pos_ref
                            while True:
                                ## read markers and phased
                                line_p = fd_phased.readline()
                                line_m = fd_markers.readline()
                                ## last marker is in the genotype probability file and not the markers file
                                ## e.g. 7:159128574 in fula4x
                                if line_m == '':
                                    bool_append_markphas = False
                                    break
                                (pos_ref, A_ref, B_ref
                                 ) = self.parse_marker(line_m)
                                if pos_ref > pos_ref_prev:
                                    break
                                else:
                                    continue
                                continue
                            pass
                        break
                    ##
                    ## 2) continue loop over genotype likelihoods ("panel 2")
                    ##
                    elif position < pos_ref:
                        bool_append_markphas = False
                        break
                    ##
                    ## 3) continue loop over markers ("panel 0")
                    ##
                    ## elif position > pos_ref:
                    else:
                        ## INDEL
                        if pos_ref == pos_prev:
                            lines_out_markers1 += [line_m]
                            lines_out_phased1 += [line_p]
                            if bool_BOF2 == True and pos_ref > pos_term1-edge:
                                lines_out_markers2 += [line_m]
        ##                            lines_out_phased2 += [line_p]
                        ## SNP unique to panel 0
        ##                    elif pos_prev < pos_ref:
                        else:
                            lines_out_markers1 += [line_m]
                            lines_out_phased1 += [line_p]
                            if bool_BOF2 == True:
                                lines_out_markers2 += [line_m]
        ##                            lines_out_phased2 += [line_p]
                        ## append marker not present in genotype likehoods file
                        ## read markers and phased
                        line_p = fd_phased.readline()
                        line_m = fd_markers.readline()
                        if line_m == '':
                            bool_append_markphas = False
                            break
                        pos_ref, A_ref, B_ref = self.parse_marker(
                            line_m)
                        continue
                    continue

                line_beagle = self.vcf2beagle(l_vcf,)

                ##
                ##
                ##
                ## BOF2
                if bool_BOF2 == True:
                    lines_out2 += [line_beagle]
                    ## match
                    if bool_append_markphas == True:
                        lines_out_markers2 += [line_m]
    ##                    lines_out_phased2 += [line_p]
                    ## mismatch
                    else: ## ms23/dg11 2013mar12
                        lines_out_markers2 += ['%s:%s %s %s %s\n' %(
                            chrom,position,position,alleleA,alleleB,)]
                if bool_EOF1 == False:
                    lines_out1 += [line_beagle]
                    ## match
                    if bool_append_markphas == True:
                        lines_out_markers1 += [line_m]
                        lines_out_phased1 += [line_p]
                    ## mismatch
                    else: ## ms23/dg11 2013mar12
                        lines_out_markers1 += ['%s:%s %s %s %s\n' %(
                            chrom,position,position,alleleA,alleleB,)]
                ## EOF1
                else:
                    ## append to previous file if less than 1000 variants
                    if bool_append_to_prev == False and len(lines_out1) < 1000:
                        ## 2nd append to prev
                        if index != 1 and pos_init1-size == d_index2pos[index-1][0]:
                            bool_append_to_prev = True
                        ## append to next
                        else:
                            bool_append_to_next = True
                    ## append to next file if first fragment
                    if bool_append_to_prev == True and index == 1:
                        bool_append_to_prev = False
                        bool_append_to_next = True
                    ## append to next file if gap (500kbp) between
                    ## first position of current lines and
                    ## last position of previous file
                    if bool_append_to_prev == True:
                        cmd = 'tail -n1 %s.%i.like | cut -d " " -f1' %(
                            fp_out_prefix,index-1)
                        pos_prev_EOF = int(os.popen(cmd).read().split(':')[1])
                        pos_curr_BOF = int(lines_out1[1].split()[0].split(':')[1])
                        if pos_curr_BOF-pos_prev_EOF > 500000:
                            bool_append_to_prev = False
                            bool_append_to_next = True
                    if bool_append_to_prev == True:
                        print('append prev')
                        mode = 'a'
                        (
                            lines_out1, lines_out_markers1,
                            ) = self.remove_duplicate_lines(
                                index, fp_out_prefix,
                                lines_out1, lines_out_markers1)
                        ##
                        index -= 1
                        if index != 0:
                            pos_init1 = d_index2pos[index][0]
                    elif bool_append_to_next == True:
                        print('append next')
                    else:
                        mode = 'w'
                    print('%2s, pos_curr %9i, pos_init %9i, pos_term %9i, i %3i, n %5i' %(
                        chrom,position,pos_init1,pos_term1,index,len(lines_out1)))
                    sys.stdout.flush()
                    if bool_append_to_next == False:
                        ## write/append current lines to file
                        fd_out = open('%s.%i.like' %(fp_out_prefix,index,),mode)
                        fd_out.writelines(lines_out1)
                        fd_out.close()
                        fd = open('%s.%i.markers' %(fp_out_prefix,index,),mode)
                        fd.writelines(lines_out_markers1)
                        fd.close()
    ##                    fd = open('%s.%i.phased' %(fp_out_prefix,index,),mode)
                        fd = open('%s.phased' %(fp_out_prefix,),'a')
                        fd.writelines(lines_out_phased1)
                        fd.close()
                        ## replace current lines with next lines
                        lines_out1 = lines_out2
    ##                    lines_out_phased1 = lines_out_phased2
                        lines_out_phased1 = []
                        if bool_append_markphas == True:
                            lines_out_phased1 += [line_p]
                        lines_out_markers1 = lines_out_markers2
                    ## append to current lines
                    ## elif bool_append_to_next == True
                    else:
                        lines_out1 += [line_beagle]
                        if bool_append_markphas == True:
                            lines_out_markers1 += [line_m]
                            lines_out_phased1 += [line_p]
                        else: ## ms23/dg11 2013mar12
                            lines_out_markers1 += ['%s:%s %s %s %s\n' %(
                                chrom,position,position,alleleA,alleleB,)]
                    ## reset next lines
                    lines_out2 = [header]
    ##                lines_out_phased2 = [header_phased]
                    lines_out_markers2 = []
                    if bool_append_to_next == False:
                        ## append index and positions to dictionary
                        d_index2pos[index] = [pos_init1,pos_term1,]
                        index += 1
                    ## set pos init/term
                    if bool_append_to_next == False:
                        pos_init1 = position-position%size
                        pos_term1 = self.calc_pos_term(pos_init1,position,size,min_bps)
                    else:
                        pos_term1 += size

                    ## reset booleans
                    bool_BOF2 = False
                    bool_EOF1 = False
                    bool_append_to_prev = False
                    bool_append_to_next = False

                    ## end of if bool_EOF == True
                    pass

                pos_prev = position

                ## read markers and phased
                if bool_append_markphas == True:
                    line_m = fd_markers.readline()
                    line_p = fd_phased.readline()
                    if line_m:
                        pos_ref, A_ref, B_ref = self.parse_marker(line_m)
                    else:
                        bool_append_markphas = False

                ## continue loop over genotype likelihoods
                continue

        ## less than minimum number of base pairs
        ## or minimum number of variants
        if position-pos_init1 < min_bps or len(lines_out1) < 1000:
            (
                lines_out1, lines_out_markers1
                ) = self.remove_duplicate_lines(
                    index, fp_out_prefix,
                    lines_out1, lines_out_markers1)
            mode = 'a'
            index -= 1
            pos_init1 = d_index2pos[index][0]
        else:
            mode = 'w'
        print('%2s, pos_curr %9i, pos_init %9i, pos_term %9i, i %3i, n %5i' %(
            chrom,position,pos_init1,pos_term1,index,len(lines_out1)))

        ## write remaining few lines to output files
        ## without checking for duplicate positions
        lines_out_markers1 += fd_markers.readlines()
        lines_out_phased1 += fd_phased.readlines()

        ##
        ## write/append lines
        ##
        fd_out = open('%s.%i.like' %(fp_out_prefix,index,),mode)
        fd_out.writelines(lines_out1)
        fd_out.close()

        fd = open('%s.%i.markers' %(fp_out_prefix,index,),mode)
        fd.writelines(lines_out_markers1)
        fd.close()

##        fd = open('%s.%i.phased' %(fp_out_prefix,index,),mode)
        fd = open('%s.phased' %(fp_out_prefix,),'a')
        fd.writelines(lines_out_phased1)
        fd.close()

        ## append position to dictionary
        d_index2pos[index] = [pos_init1,pos_term1,]

        ##
        ## close all files after looping over last line
        ##
        fd_markers.close()
        fd_phased.close()

        with open('summaryVR.txt','a') as file_summary:
            file_summary.write('%s %i\n' %(chrom,cnt_variants))
##        print('%s %i\n' %(chrom,cnt_variants),file='summaryVR.txt')

        ##
        ## file i/o checks
        ##
##        self.BEAGLE_divide_fileIO_checks(chrom,fp_phased,)
    
        return d_index2pos


    def hook_compressed_text(self, filename, mode):

        ext = os.path.splitext(filename)[1]
        if ext == '.gz':
            f = gzip.open(filename, mode + 't')
        else:
            f = open(filename, mode)

        return f


    def vcf2beagle(self,l_vcf):

        ## this function is very slow; especially split, sum and pow

        ## append chrom:pos alleleA alleleB
        line_beagle = '%s:%s %s %s' %(l_vcf[0],l_vcf[1],l_vcf[3],l_vcf[4],)

        index = l_vcf[8].split(':').index('PL')
        for s_vcf in l_vcf[9:]:
            ## variant not called
            if s_vcf[:3] == './.':
                line_beagle += ' 0.3333 0.3333 0.3333'
                continue
            l_probs = []
            l_log10likelihoods = s_vcf.split(':')[index].split(',')
            for log10likelihood in l_log10likelihoods:
                log10likelihood = int(log10likelihood)
                if log10likelihood == 0:
                    prob = 1
                elif log10likelihood > 50:
                    prob = 0
                else:
                    prob = pow(10,-log10likelihood/10)
                l_probs += [prob]
            ## append normalized probabilities
            for prob in l_probs:
                line_beagle += ' %6.4f' %(prob/sum(l_probs))
            if ',' in l_vcf[4]:
                print(line_beagle)
                print(s_vcf)
                print(l_log10likelihoods)
                print(l_probs)
                print(line_beagle)
                stop_tmp_decide_on_triallelic
        line_beagle += '\n'

        return line_beagle


    def parse_variant_type(self,l_vcf):

        ## skip deletions
        if not l_vcf[3] in ['A','C','G','T',]:
            return True,None
        ## dialleic SNP
        if l_vcf[4] in ['A','C','G','T',]:
            return False,True
        ## triallelic and other non-diallelic SNPs
##        elif l_vcf[4] in [
##            'A,C','A,G','A,T','C,G','C,T','G,T',
##            'A,C,G','A,C,T','A,G,T','C,G,T',
##            ]:
##            return False,False,False
##        ## UG DISCOVERY
##        elif l_vcf[4] in iter(','.join(tup) for tup in itertools.combinations('ACGT',2)):
        ## UG GENOTYPE_GIVEN_ALLELES
        elif l_vcf[4] in iter(','.join(tup) for tup in itertools.permutations('ACGT',2)):
            return False,False
##        ## UG DISCOVERY
##        elif l_vcf[4] in iter(','.join(tup) for tup in itertools.combinations('ACGT',3)):
        ## UG GENOTYPE_GIVEN_ALLELES
        elif l_vcf[4] in iter(','.join(tup) for tup in itertools.permutations('ACGT',3)):
            return False,False
        ## skip insertions
        else:
            return True,None

        return


    def BEAGLE_divide_fileIO_checks(self,chrom,fp_phased,):

        '''Compare VR out and BEAGLE in'''

        cmd = 'cat %s' %(fp_phased)
        cmd += "| awk 'NR>1{print $2}' | sort -u > panel0in%s" %(chrom)
        self.execmd(cmd)

        cmd = "awk 'FNR>1{print $2}' in_BEAGLE/%s/%s.phased" %(chrom,chrom)
        cmd += "| sort -u > panel0out%s" %(chrom)
        self.execmd(cmd)

        fp_tranches = 'out_VariantRecalibrator/VariantRecalibrator.SNP.tranches'
        minVQSLOD = self.parse_minVQSLOD(fp_tranches,self.ts_filter_level,)
        fp_recal = 'out_VariantRecalibrator/VariantRecalibrator.SNP.recal'

##        with contextlib.ExitStack() as stack:
##            fd_recal = stack.enter_context(open(fp_recal))
##            fd_out = stack.enter_context(open('panel2in%s' %(chrom),'w'))
##            fd_UG = stack.enter_context(fileinput.fileinput(files_sorted))
##            l_files = self.sort_nicely(
##                glob.glob('out_UnifiedGenotyper/%s.*.vcf.gz' %(chrom)))
##            assert len(l_files) > 0
##                l_files_sorted = self.sort_nicely(l_files)
##            fd_vcf = fileinput.input(
##                files=l_files,openhook=self.hook_compressed_text)
##            for l_vcf in self.generate_line_vcf_PASS_split(
##                fd_vcf,fd_recal,minVQSLOD,))

        with open(fp_recal,'r') as fd_recal, open('panel2in%s' %(chrom),'w') as fd_out:
            for line in fd_recal:
                if line[0] == '#': continue
                VQSLOD, l_recal = self.parse_VQSLOD(line)
                if VQSLOD < minVQSLOD:
                    continue
                ## slow conversion to integer just to make sure it is an integer
                pos = int(l_recal[1])
                fd_out.write('%s:%i\n' %(chrom,pos))
        self.execmd('sort panel2in%s -o panel2in%s' %(chrom,chrom))

        cmd = "awk 'FNR>1{print $1}' in_BEAGLE/%s/%s.*.like" %(chrom,chrom)
        cmd += ' | sort -u > panel2out%s' %(chrom)
        self.execmd(cmd)

        cmd = "awk '{print $1}' in_BEAGLE/%s/%s.*.markers" %(chrom,chrom)
        cmd += ' | sort -u > markers%s' %(chrom)
        self.execmd(cmd)

        ## check that like (panel2out) is identical to PBI (panel2in)
        cmd = 'comm -3 panel2out%s panel2in%s' %(chrom,chrom)
        i = int(os.popen('%s | wc -l' %(cmd)).read())
        if i > 0:
            print(os.popen('%s | head' %(cmd)).readlines())
            print(cmd)
            sys.exit()

        ## check that panel0out is a subset of panel0in
        cmd = 'comm -23 panel0out%s panel0in%s' %(chrom,chrom)
        i = int(os.popen('%s | wc -l' %(cmd)).read())
        if i > 0:
            print(os.popen('%s | head' %(cmd)).readlines())
            print(cmd)
            sys.exit()

        ## check that like (panel2out) is a true subset of markers
        cmd = 'comm -32 panel2out%s markers%s' %(chrom,chrom)
        i = int(os.popen('%s | wc -l' %(cmd)).read())
        if i > 0:
            print(os.popen('%s | head' %(cmd)).readlines())
            print(cmd)
            sys.exit()

        ## check that phased (panel0out) is a true subset of markers
        cmd = 'comm -32 panel0out%s markers%s' %(chrom,chrom)
        i = int(os.popen('%s | wc -l' %(cmd)).read())
        if i > 0:
            print(os.popen('%s | head' %(cmd)).readlines())
            print(cmd)
            sys.exit()

        ## check that markers is identical to panel2in+panel0in
        cmd = 'cat panel2out%s panel0out%s | sort -u > panelsout' %(chrom,chrom)
        self.execmd(cmd)
        cmd = 'comm -3 panelsout markers%s' %(chrom)
        i = int(os.popen('%s | wc -l' %(cmd)).read())
        if i > 0:
            print(os.popen(cmd).readlines()[:10])
            print(cmd)
            sys.exit()

        for affix in ['panel2out','panel2in','panel0in','panel0out','markers',]:
            os.remove('%s%s' %(affix,chrom,))
        os.remove('panelsout')

        return


    def remove_duplicate_lines(
        self, index, fp_out_prefix,
        lines_out1, lines_out_markers1):

        ##
        ## remove headers
        ##
        lines_out1 = lines_out1[1:]
##        lines_out_phased1 = lines_out_phased1[1:]

        ##
        ## get terminal positions in previous files
        ##
        cmd = 'tail -n1 %s.%i.like | cut -d " " -f1' %(
            fp_out_prefix,index-1)
        pos_prev_like = int(os.popen(cmd).read().split(':')[1])

        cmd = 'tail -n1 %s.%i.markers' %(
            fp_out_prefix,index-1)
        l = os.popen(cmd).read().strip().split()
        pos_prev_markers = int(l[1])
        alleleA_prev = l[2]
        alleleB_prev = l[3]

##        cmd = 'tail -n1 %s.%i.phased | cut -d " " -f2' %(
##            fp_out_prefix,index-1)
##        pos_prev_phased = int(os.popen(cmd).read().split(':')[1])

        ##
        ## index position in current lines
        ## and remove duplicate lines
        ##
        bool_found = False
        for i in range(len(lines_out1)):
            if pos_prev_like == int(lines_out1[i].split()[0].split(':')[1]):
                bool_found = True
                break
        if bool_found == True:
            lines_out1 = lines_out1[i+1:]

        bool_found = False
        for i in range(len(lines_out_markers1)):
            l = lines_out_markers1[i].split()
            alleleA = l[2]
            alleleB = l[3]
            if pos_prev_markers == int(l[1]):
                bool_found = True
                break
        if bool_found == True:
            lines_out_markers1 = lines_out_markers1[i+1:]

##        bool_found = False
##        for i in xrange(len(lines_out_phased1)):
##            if pos_prev_phased == int(lines_out_phased1[i].split()[1].split(':')[1]):
##                bool_found = True
##                break
##        if bool_found == True:
##            lines_out_phased1 = lines_out_phased1[i+1:]

        return lines_out1, lines_out_markers1


    def calc_pos_term(self,pos_init1,position,size,min_bps):

        ## acrocentric chromosome telomeres (13,14,15,21,22,Y)
        ## and various centromeres (e.g. chromosome 6)
        if size-min_bps < position-pos_init1:
            pos_term1 = pos_init1+2*size
        ## non-acrocentric chromosome
        ## and various centromeres
        else:
            pos_term1 = pos_init1+1*size

        return pos_term1


    def BEAGLE(self,d_chrom_lens,):

        ## http://faculty.washington.edu/browning/beagle/beagle_3.3.2_31Oct11.pdf

        ## http://www.broadinstitute.org/gatk/guide/article?id=43
        ## Interface with BEAGLE imputation software - GSA
        ## IMPORTANT: Due to BEAGLE memory restrictions,
        ## it's strongly recommended that BEAGLE be run on a separate chromosome-by-chromosome basis.
        ## In the current use case, BEAGLE uses RAM in a manner approximately proportional to the number of input markers.

        ## "For haplotype phase inference and imputation of missing data
        ## with default BEAGLE options,
        ## memory usage increases with the number of markers."
        ## but not bloody linear...
        memMB = 3900 ## 100 samples
        memMB = 4900 ## 320 samples
        ## fgrep CPU */stdout/BEAGLE/*.out | sort -k5nr,5 | head -n1
        ## zulu_20121208/stdout/BEAGLE/BEAGLE.19.1.out:    CPU time   :  80270.55 sec.
        queue = 'long'

        ##
        ## 1) check input existence
        ##
        mode = 'SNP'
        fp_in_recal = 'out_VariantRecalibrator/VariantRecalibrator.%s.recal' %(mode)
        fp_in_tranches = 'out_VariantRecalibrator/VariantRecalibrator.%s.tranches' %(mode)
        if self.check_in(
            'VariantRecalibrator',[fp_in_recal,fp_in_tranches,],
            'touch/VariantRecalibrator.SNP.touch'):
            sys.exit()
        stoptmp

        ##
        ## 2) touch
        ##
        bool_return = self.touch('BEAGLE')
        if bool_return == True: return

        ##
        ## write shell script
        ##
        self.BEAGLE_write_shell_script(memMB)

        ##
        ## split into fragments prior to imputation
        ##
        d_indexes = {}

        if not (
            os.path.isfile('BEAGLE_divide_indexes.txt')
            and
            os.path.isdir('in_BEAGLE')
            ):
            d_indexes = self.BEAGLE_divide()
            s = ''
            for chrom,d_index2pos in d_indexes.items():
                for index,[pos1,pos2,] in d_index2pos.items():
                    s += '%s %i %i %i\n' %(chrom,index,pos1,pos2,)
            fd = open('BEAGLE_divide_indexes.txt','w')
            fd.write(s)
            fd.close()
        else:
            print('BEAGLE_divide_indexes.txt exists')
            d_indexes = {}
            for chrom in self.chroms:
                d_indexes[chrom] = {}
            fd = open('BEAGLE_divide_indexes.txt','r')
            lines = fd.readlines()
            fd.close()
            for line in lines:
                l = line.strip().split()
                chrom = l[0]
                if not chrom in self.chroms:
                    continue
                index = int(l[1])
                pos1 = int(l[2])
                pos2 = int(l[3])
                d_indexes[chrom][index] = [pos1,pos2,]
                continue
            pass

        ##
        ## execute shell script
        ##
        for chrom in self.chroms:

            print('bsub BEAGLE %s' %(chrom))

##            J = '%s%s[%i-%i]' %('BEAGLE',chrom,1,max(d_indexes[chrom].keys()),)
            for index in d_indexes[chrom].keys():

                ## finished?
                fn_out = 'out_BEAGLE/%s/%s.%i.like.gprobs.gz' %(
                    chrom,chrom,index)
                if os.path.isfile(fn_out):
                    continue

                ## started?
                fn = 'LSF/BEAGLE/%s.%i.out' %(chrom,index)
                if os.path.isfile(fn):
                    if os.path.getsize(fn):
                        with open(fn) as f:
                            if 'iteration' in f.readlines()[-1]:
                                continue

                print(chrom,index)

                J = '%s.%s[%s-%s]' %('BEAGLE',chrom,index,index,)
                LSF_affix = '%s/%s.%%I' %('BEAGLE',chrom)
                cmd = self.bsub_cmd(
                    'BEAGLE', J, memMB=memMB, LSF_affix=LSF_affix,
                    chrom=chrom, queue=queue,)
                os.system(cmd)

        return


    def init_java(self, jar, memMB, java='java', bool_checkpoint=False):

        s = '%s -Djava.io.tmpdir=%s' %(java,'tmp')
        ## set maximum heap size
        s += ' -Xmx%im' %(memMB)
        if bool_checkpoint:
            s += ' -XX:-UsePerfData -Xrs '
        s += ' \\\n -jar %s' %(jar)

        return s


    def BEAGLE_write_shell_script(self,memMB,):

        fp_out = 'out_BEAGLE/$CHROMOSOME/$CHROMOSOME.${LSB_JOBINDEX}.like'

        ## initiate shell script
        lines = ['#!/bin/bash\n']

        ## parse chromosome from command line
        lines += ['CHROMOSOME=$1\n']

        lines += ['if [ -s %s.gprobs.gz ]; then exit; fi\n' %(fp_out)] ## redundant

##        ## init cmd
##        lines += ['cmd="']

        ##
        ## initiate BEAGLE
        ##
        s_java = self.init_java(self.fp_software_beagle,memMB)
        lines += ['%s \\' %(s_java)]

        lines += self.body_BEAGLE(fp_out,)

        ## term cmd
        lines += self.term_cmd('BEAGLE',['%s.gprobs.gz' %(fp_out)],)

        ## write shell script
        self.write_shell('shell/BEAGLE.sh',lines,)

        return


    def body_BEAGLE(self,fp_out,):

        fp_like = 'in_BEAGLE/$CHROMOSOME/$CHROMOSOME.${LSB_JOBINDEX}.like'
##        fp_phased = 'in_BEAGLE/$CHROMOSOME/$CHROMOSOME.${LSB_JOBINDEX}.phased'
        fp_phased = 'in_BEAGLE/$CHROMOSOME/$CHROMOSOME.phased'
        fp_markers = 'in_BEAGLE/$CHROMOSOME/$CHROMOSOME.${LSB_JOBINDEX}.markers'

        lines = []

##like=<unphased likelihood data file> where <unphased likelihood data file> is
##the name of a genotype likelihoods file for unphased, unrelated data
##(see Section 2.2). You may use multiple like arguments if data from different
##cohorts are in different files.
        lines += [' like=%s \\' %(fp_like)]
####arguments for phasing and imputing data ...
#### Arguments for specifying files
## phased=<phased unrelated file> where <phased unrelated file> is the name of a
## Beagle file containing phased unrelated data (see Section 2.1).
## You may use multiple phased arguments if data from different cohorts are in
## different files.
##        lines += [' phased=in_BEAGLE/ALL.chr$CHROMOSOME.phase1_release_v3.20101123.filt.renamed.bgl \\']
        lines += [' phased=%s \\' %(fp_phased)]
####  unphased=<unphased data file>                     (optional)
####  phased=<phased data file>                         (optional)
####  pairs=<unphased pair data file>                   (optional)
####  trios=<unphased trio data file>                   (optional)
####  like=<unphased likelihood data file>              (optional)
##markers=<markers file> where <markers file> is the name of the markers file containing
## marker identifiers, positions, and alleles described in Section 2.4.
## The markers argument is optional if you specify only one Beagle file,
## and is required if you specify more than one Beagle file.
        lines += [' markers=%s \\' %(fp_markers)]
####missing=<missing code> where <missing code> is the character or sequence of characters used to represent a missing allele (e.g. missing=-1 or missing=?).
#### The missing argument is required.
##        s += ' missing=? '

##nsamples=<number of samples> where <number of samples> is positive integer
##giving the number of haplotype pairs to sample for each individual during each
##iteration of the phasing algorithm. The nsamples argument is optional. The
##default value is nsamples=4. If you are phasing an extremely large sample
## (say > 4000 individuals), you may want to use a smaller nsamples parameter
## (e.g. 1 or 2) to reduce computation time.  If you are phasing a small sample
## (say < 200 individuals), you may want to use a larger nsamples parameter
## (say 10 or 20) to increase accuracy.
        lines += [' nsamples=%i \\' %(int(self.i_BEAGLE_nsamples))]

        lines += [' niterations=10 \\'] ## default 10

        lines += [' omitprefix=true \\'] ## default false

        lines += [' verbose=false \\'] ## default false

##        lines += [' lowmem=false \\'] ## default false
        lines += [' lowmem=true \\'] ## default false

        ## non-optional output prefix
        lines += [' out=%s \\' %(fp_out)]

        return lines


    def parse_chrom_lens(self):

        d_chrom_lens = {}

        ## 1000G chromosome ranges
        fn = '%s.fai' %(self.reference_sequence)
        fd = open(fn)
        lines = fd.readlines()
        fd.close()
        for line in lines:
            l = line.strip().split()
            chrom = l[0]
            chrom_len = int(l[1])
##            if not chromosome in ['X','Y',]:
##                chrom = int(chromosome)
            d_chrom_lens[chrom] = chrom_len
            ## break, when we reach the Y chromosome
            if chrom == 'Y':
                break

        return d_chrom_lens


    def check_in(self, analysis_type, l_fp_in, fp_touch,):

        d_l_fp_out = {}
        
        with open(fp_touch) as fd:
            s = fd.read()
        l_fp_out = s.split('\n')
        d_l_fp_out['touch'] = l_fp_out

        ## todo: use os.walk here instead...
        for dirname in ['',]:
            d_l_fp_out[dirname] = []
            l = os.listdir(os.path.join(dirname,'out_%s' %(analysis_type)))
            for s in l:
                path1 = os.path.join('out_%s' %(analysis_type),s)
                path2 = os.path.join(dirname,path1)
                ## append files in chromosomal subdirectories
                if os.path.isdir(path2):
                    l = os.listdir(path2)
                    for fn in l:
                        d_l_fp_out[dirname] += [os.path.join(path1,fn)]
                ## append files in main dir
                elif os.path.isfile(path2):
                    d_l_fp_out[dirname] += [path1]
                else:
                    print(path2)
                    print(os.path.realpath(path2))
                    print(os.path.isfile(os.path.realpath(path2)))
                    stop_not_expected

        bool_exit = False
        for dirname,l_fp_out in d_l_fp_out.items():
            if len(set(l_fp_in)-set(l_fp_out)) > 0:
                print('%s and possibly %s other files not generated.' %(
                    list(set(l_fp_in)-set(l_fp_out))[0],
                    len(set(l_fp_in)-set(l_fp_out))-1,))
                print('dirname', dirname)
                print('%s has not run to completion. Exiting.' %(analysis_type))
                bool_exit = True
#                print(inspect.stack()[1])
##                sys.exit()

        return bool_exit


    def CombineGVCFs(self):

        '''Merge gVCFs prior to GenotypeGVCFs'''

        analysis_type = T = 'CombineGVCFs'
        memMB = 15900
        queue = 'long'

        ## 1) check input existence / check that previous jobs finished
        l_vcfs_in = [
            'out_HaplotypeCaller/%s.vcf.gz' %(
                os.path.splitext(os.path.basename(bam))[0])
            for bam in sorted(self.bams)]
        if self.check_in(
            'HaplotypeCaller', ['%s.tbi' %(vcf) for vcf in l_vcfs_in],
            'touch/HaplotypeCaller.touch'):
            sys.exit()

        ## write shell script
        self.shell_CombineGVCFs(T, memMB)

        for chrom in self.chroms:

            ## 2) check output existence / check that job did not start
            if self.touch('%s.%s' %(analysis_type, chrom)):
                return

            self.mkdir('lists')
            l_combined = []
            for fn_list in glob.glob('lists/CombineGVCFs.%s.*.list' %(chrom)):
                with open(fn_list) as f:
                    l_combined += f.read().rstrip().split('\n')

            l_vcfs_in = list(sorted(list(set(l_vcfs_in)-set(l_combined))))

            for i, vcf in enumerate(
                l_vcfs_in,
                self.gVCF_limit*len(glob.glob(
                    'lists/CombineGVCFs.%s.*.list' %(chrom)))):
                if i%self.gVCF_limit == 0:
                    fn_out = 'lists/CombineGVCFs.{chrom}.{i}.list'.format(
                        chrom=chrom, i=i//self.gVCF_limit)
                    assert not os.path.isfile(fn_out)
                    fd_out = open(fn_out, 'w')
                fd_out.write('{}\n'.format(vcf))
            fd_out.close()

            self.mkdir('LSF/%s' %(T))
            for i in range(len(glob.glob(
                'lists/CombineGVCFs.%s.*.list' %(chrom)))):
                ## skip if job initiated
                if os.path.isfile('out_CombineGVCFs/%s.%i.vcf.gz' %(chrom,i)):
                    continue
                cmd = self.bsub_cmd(
                    T, 'CgVCFs.{}.{}'.format(chrom, i),
                    memMB=memMB, queue=queue,
                    LSF_affix = '%s/%s.%s' %(T, chrom, i)
                    )
                cmd += ' {} {}'.format(chrom, i)
                self.execmd(cmd)

        return


    def shell_CombineGVCFs(self, T, memMB):

        lines = ['#!/bin/bash\n']
        lines += ['chrom=$1']
        lines += ['i=$2']
        lines += ['out=out_CombineGVCFs/$chrom.$i.vcf.gz']
        lines += ['## exit if job started']
        lines += ['if [ -s $out ]; then exit; fi\n']
        ## make output directory
        lines += ['mkdir -p $(dirname $out)\n']

        lines += self.init_GATK_cmd(T, memMB,)
        lines += [' -L $chrom \\']
        lines += [' -V lists/CombineGVCFs.$i.list \\']
        lines += [' -o $out \\']

        ## terminate shell script
        lines += self.term_cmd(T, ['$out.tbi'],)

        ## write shell script
        self.write_shell('shell/%s.sh' %(T), lines,)

        return


    def GenotypeGVCFs(self):

        '''Convert gVCFs to VCFs'''

        analysis_type = T = 'GenotypeGVCFs'
        memMB = 7900
        queue = 'basement'

        ## write shell script
        self.shell_GenotypeGVCFs(T, memMB)
        self.mkdir('LSF/%s' %(T))

        for chrom in self.chroms:

            ## 1) check input existence / check that previous jobs finished
            l_vcfs_in = [
                'out_CombineGVCFs/%s.%i.vcf.gz' %(chrom,i)
                for i in range(len(self.bams)//split+1)]
            if self.check_in(
                'CombineGVCFs', ['%s.tbi' %(vcf) for vcf in l_vcfs_in],
                'touch/CombineGVCFs.touch'):
                ## continue loop over chromosomes
                continue

            ## 2) check output existence / check that job did not start
            if self.touch('%s.%s' %(analysis_type, chrom)):
                continue
        
            self.mkdir('lists')
            with open(
                'lists/GenotypeGVCFs.{chrom}.list'.format(
                    chrom=chrom), 'w') as f:
                for vcf in l_vcfs_in:
                    f.write('{}\n'.format(vcf))

                cmd = self.bsub_cmd(
                    T, 'GgVCFs.{}'.format(chrom),
                    memMB=memMB, queue=queue,
                    LSF_affix = '%s/%s' %(T, chrom)
                    )
                cmd += ' {}'.format(chrom)
                self.execmd(cmd)

        return


    def shell_GenotypeGVCFs(self, T, memMB):

        lines = ['#!/bin/bash\n']
        lines += ['chrom=$1']
        lines += ['out=out_GenotypeGVCFs/$chrom.vcf.gz']
        lines += ['## exit if job started']
        lines += ['if [ -s $out ]; then exit; fi\n']
        lines += ['## exit if job finished']
        lines += ['if [ -s $out.tbi ]; then exit; fi\n']
        ## make output directory
        lines += ['mkdir -p $(dirname $out)\n']

        lines += self.init_GATK_cmd(T, memMB,)
        lines += [' -L $chrom \\']
        lines += [' -V lists/{}.$chrom.list \\'.format(T)]
        lines += [' -o out_{}/$chrom.vcf.gz \\'.format(T)]

        ## terminate shell script
        lines += self.term_cmd(T, ['$out.tbi'])

        ## write shell script
        self.write_shell('shell/%s.sh' %(T), lines,)

        return


    def VariantRecalibrator(self,d_chrom_lens):

        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_variantrecalibration_VariantRecalibrator.html

        T = analysis_type = 'VariantRecalibrator'
        memMB = 20900
        queue = 'yesterday'
        num_threads = 4

        ##
        ## 1) check input existence (vcf)
        ##
        l_vcfs_in = [
            'out_HaplotypeCaller/%s.vcf.gz' %(
                os.path.splitext(os.path.basename(bam))[0])
            for bam in self.bams]
        if self.check_in(
            'HaplotypeCaller',l_vcfs_in, 'touch/HaplotypeCaller.touch'):
            sys.exit(0)
        stop1

        d_resources = {'SNP':self.fp_resources_SNP,'INDEL':self.fp_resources_INDEL,}

        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_variantrecalibration_VariantRecalibrator.html#--mode
        for mode in ['SNP','INDEL',]:

            memMB = {'SNP':20900,'INDEL':8900}[mode]

            ## 2) touch / check output
            bool_continue = self.touch('%s.%s' %(analysis_type,mode))
            if bool_continue == True:
                continue

            ## Define file paths.
            fp_tranches = 'out_VariantRecalibrator/VariantRecalibrator.%s.tranches' %(mode)
            fp_recal = 'out_VariantRecalibrator/VariantRecalibrator.%s.recal' %(mode)

            ## Initiate GATK walker.
            lines = self.init_GATK_cmd(analysis_type,memMB,)

            lines += [' --num_threads %i \\' %(num_threads)]

            ## required, in
            lines += [' --input %s \\' %(vcf) for vcf in l_vcfs_in]
            ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_variantrecalibration_VariantRecalibrator.html#--use_annotation
            if mode == 'SNP':
                lines += [
                    ' --use_annotation QD \\',
                    ' --use_annotation HaplotypeScore \\',
                    ' --use_annotation MQRankSum \\',
                    ' --use_annotation ReadPosRankSum \\',
                    ' --use_annotation MQ \\',
                    ' --use_annotation FS \\',
                    ' --use_annotation DP \\',
                    ]
            elif mode == 'INDEL':
                lines += [' -an DP -an FS -an ReadPosRankSum -an MQRankSum \\',]

            ##
            ## required, out
            ##
            lines += [' --recal_file %s \\' %(fp_recal)]
            lines += [' --tranches_file %s \\' %(fp_tranches)]

            ##
            ## Optional Parameters.
            ##
            ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_variantrecalibration_VariantRecalibrator.html#--mode
            lines += [' --mode %s \\' %(mode)]

            ## http://gatkforums.broadinstitute.org/discussion/1259/what-vqsr-training-sets-arguments-should-i-use-for-my-specific-project
            if mode == 'INDEL':
                ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_variantrecalibration_VariantRecalibrator.html#--maxGaussians
                lines += [' --maxGaussians 4 \\'] ## default 8
                ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_variantrecalibration_VariantRecalibrator.html#--minNumBadVariants
                lines += [' --minNumBadVariants 1000 \\'] ## default 1000

            ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_variantrecalibration_VariantRecalibrator.html#--resource
            fd = open(d_resources[mode],'r')
            lines_resources = fd.readlines()
            fd.close()
            lines += [' %s \\' %(line.strip()) for line in lines_resources]

            ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_variantrecalibration_VariantRecalibrator.html#--TStranche
            l_TStranches = []
    ##        l_TStranches += [99.70+i/20. for i in range(6,0,-1,)]
            l_TStranches += [99+i/10. for i in range(10,0,-1,)]
            l_TStranches += [90+i/2. for i in range(18,-1,-1,)]
##            l_TStranches += [70+i for i in range(19,-1,-1,)]
            s_TStranches = ''
            for TStranche in l_TStranches:
                s_TStranches += '--TStranche %.1f ' %(TStranche)
            lines += [' %s \\' %(s_TStranches)]

            ##
            ## GATKwalker, optional, out
            ##
            lines += [' --rscript_file out_%s/%s.%s.plots.R \\' %(T,T,mode,)]

            ## Terminate command and rerun pipeline.
            lines += self.term_cmd(
                '%s.%s' %(analysis_type,mode),[fp_tranches,fp_recal,],)

            self.write_shell('shell/%s.%s.sh' %(analysis_type,mode,),lines,)

            J = 'VR.%s' %(mode)
            cmd = self.bsub_cmd(
                '%s.%s' %(analysis_type,mode), J, memMB=memMB, queue=queue,
                LSF_affix='%s/%s' %(analysis_type,mode,),
                num_threads = num_threads)
            self.execmd(cmd)

        return


    def touch(self,analysis_type):

        bool_return = False
        fn_touch = 'touch/%s.touch' %(analysis_type)
        if os.path.isfile(fn_touch):
            if self.verbose == True:
                print('in progress or completed:', analysis_type)
            bool_return = True
        else:
            self.execmd('touch %s' %(fn_touch))

        return bool_return


    def write_brestart(self,):

        with open('brestart.sh','w') as f:
            f.write('sleep 30\n')
            f.write("IFS=$'\\n'\n")
            f.write('jobID=$1\n')
            f.write('project=$2\n')
            f.write('memMB=$3\necho memMB $memMB\n')
            f.write('pwd=$(pwd)\n')
            f.write('bhist=$(bhist -l $jobID)\n')
            f.write('i=$(echo $bhist | fgrep TERM_RUNLIMIT | wc -l)\n') ## exit code 140, run limit
            f.write('j=$(echo $bhist | fgrep "Exited with exit code 16" | wc -l)\n') ## exit code 16, pid taken
            f.write('k=$(echo $bhist | fgrep "Checkpoint failed" | wc -l)\n')
            f.write('l=$(echo $bhist | fgrep "Done successfully" | wc -l)\n')
            f.write('if [ $l -eq 0 ]; then\n')
            f.write('if [ $i -eq 0 -a $j -eq 0 -a $k -eq 0 ]; then echo $bhist >> bhist_unexpectederror_or_success.tmp; exit; fi\n')
            f.write('fi\n')
            f.write('s=$(brestart -G $project -M$memMB $pwd/checkpoint/$jobID)\n')
            f.write('echo s $s\n')
            f.write('if [ $k -ne 0 ]; then echo $s >> checkpointfailed_brestartout.tmp; fi\n')
            f.write('''jobID=$(echo $s | awk -F "[<>]" '{print $2}')\n''')
            f.write('echo jobID $jobID\n')
            f.write('echo memMB $memMB\n')
#            f.write("bsub -R 'select[mem>900] rusage[mem=900]' -M900 \\\n")
            f.write("bsub -R 'select[mem>'$memMB'] rusage[mem='$memMB']' -M$memMB \\\n")
            f.write(' -o tmp_brestart2.out -e tmp_brestart2.err \\\n')
            f.write(' -G $project -q normal -w "ended($jobID)" \\\n')
            f.write(' bash brestart.sh $jobID $project $memMB\n')

        return


    def HaplotypeCaller(self,d_chrom_lens,):

        T = analysis_type = 'HaplotypeCaller'
        queue = 'long'
        memMB = 3900
#        queue = 'basement'

        ## 1) touch/lock
        if self.touch(analysis_type):
            return

        ## 2) write shell script
        self.shell_HC(analysis_type, memMB)

        ## Create folders.
        self.mkdir('LSF/%s' %(T))
        self.mkdir('out_HaplotypeCaller/')
        self.mkdir('touch/HaplotypeCaller/')
        self.mkdir('tmp/')

        ## 3) execute shell script
        for bam in self.bams:
            basename = os.path.splitext(os.path.basename(bam))[0]

            ## finished?
            if os.path.isfile('out_HaplotypeCaller/%s.vcf.gz.tbi' %(basename)):
                continue

            ## in progress?
            pathLSF = 'LSF/HaplotypeCaller/%s' %(basename)
            ## job finished
            if os.path.isfile('out_%s/%s.vcf.gz.tbi' %(T, basename)):
                continue
            if os.path.isfile('%s.err' %(pathLSF)):
                ## file changed within the past 5 minutes?
                if time.time()-os.path.getmtime('%s.err' %(pathLSF)) < 300:
                    continue
                os.remove('%s.err' %(pathLSF))
            if os.path.isfile('%s.out' %(pathLSF)):
                os.remove('%s.out' %(pathLSF))

            J = '%s %s' %('HC', basename)
            LSF_affix = '%s/%s' %(T, basename)
            cmd = self.bsub_cmd(
                analysis_type, J, LSF_affix=LSF_affix,
                memMB=memMB, queue=queue, num_threads=self.nct,
                bam=bam)
            self.execmd(cmd)

        return


    def mkdir(self, path):

        if not os.path.splitext(path)[1]:
            dirname = path
        else:
            dirname = os.path.dirname(path)
        if not os.path.isdir(dirname):
            os.makedirs(dirname)

        return


    def shell_HC(self, analysis_type, memMB):

        lines = ['#!/bin/bash\n']
        lines += ['BAM=$1']
        lines += ['BAMBASENAME=$(basename $BAM | rev | cut -d "." -f2- | rev)']
        lines += ['out=out_HaplotypeCaller/$BAMBASENAME.vcf.gz']
##        ## exit if job started
##        lines += ['if [ -s $out ]; then exit; fi\n']
        ## exit if job finished
        lines += ['if [ -s $out.tbi ]; then exit; fi\n']
        ## make output directory
        lines += ['mkdir -p $(dirname $out)\n']

        ## initiate GATK command
        lines += self.init_GATK_cmd(
            analysis_type, memMB, bool_checkpoint=self.bool_checkpoint)

        ## append GATK command options
        lines += self.body_HaplotypeCaller()

        ## terminate shell script
##        lines += self.term_cmd(
##            analysis_type, ['$out.tbi'], extra='tabix -p vcf $out')
        lines += self.term_cmd(analysis_type, ['$out.tbi'])

        ## write shell script
        self.write_shell('shell/%s.sh' %(analysis_type),lines,)

        return


    def body_HaplotypeCaller(self,):

        '''Write GATK HaplotypeCaller specific command line arguments.'''

        ## https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_haplotypecaller_HaplotypeCaller.html

        lines = []

        ##
        ## Inherited arguments
        ##

        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_CommandLineGATK.html#--input_file
        lines += [' --input_file $BAM \\']

        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_CommandLineGATK.html#--intervals
        if self.intervals:
            lines += ['--intervals %s \\' %(self.intervals)]
            ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_CommandLineGATK.html#--interval_set_rule
            lines += ['--interval_set_rule INTERSECTION \\']
            pass

        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_CommandLineGATK.html#--variant_index_parameter
        ## http://gatkforums.broadinstitute.org/discussion/3893/calling-variants-on-cohorts-of-samples-using-the-haplotypecaller-in-gvcf-mode
        lines += [' --variant_index_parameter 128000 \\']

        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_CommandLineGATK.html#--variant_index_type
        lines += [' --variant_index_type LINEAR \\']

        ##
        ## Optional Inputs
        ##

        ## https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_haplotypecaller_HaplotypeCaller.html#--alleles
        if self.alleles:
            lines += [' --alleles %s \\' %(self.alleles)]

        ## dbSNP file. rsIDs from this file are used to populate the ID column of the output.
        ## https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_haplotypecaller_HaplotypeCaller.html#--dbsnp
        if self.dbsnp:
            lines += [' --dbsnp %s \\' %(self.dbsnp)]

        ##
        ## Optional Outputs
        ##

        ## https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_haplotypecaller_HaplotypeCaller.html#--out
        lines += [' --out $out \\']

        ##
        ## Optional Parameters
        ##

        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_haplotypecaller_HaplotypeCaller.html#--genotyping_mode
        lines += [' -gt_mode %s \\' %(self.genotyping_mode)] ## default value DISCOVERY

        ## https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_haplotypecaller_HaplotypeCaller.html#--standard_min_confidence_threshold_for_calling
        ## https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_haplotypecaller_HaplotypeCaller.html#--standard_min_confidence_threshold_for_emitting
        if self.coverage > 10:
            lines += [' -stand_call_conf 30 \\']
            lines += [' -stand_emit_conf 30 \\']
        elif len(self.bams) > 100:
            lines += [' -stand_call_conf 10 \\']
            lines += [' -stand_emit_conf 10 \\']
        else:
            print(self.project)
            stop
            lines += [' -stand_call_conf 4 \\']
            lines += [' -stand_emit_conf 4 \\']

        ##
        ## Advanced Parameters
        ##
        
        ## http://www.broadinstitute.org/gatk/gatkdocs/#VariantAnnotatorannotations
        ## https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_haplotypecaller_HaplotypeCaller.html#--annotation
        s_annotation = ''
        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_annotator_Coverage.html
        ## http://gatkforums.broadinstitute.org/discussion/2318/undocumented-change-in-2-4-a-depthofcoverage
        s_annotation += ' --annotation Coverage'
        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_annotator_FisherStrand.html
        s_annotation += ' -A FisherStrand'
        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_annotator_HaplotypeScore.html
        s_annotation += ' -A HaplotypeScore'
        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_annotator_MappingQualityRankSumTest.html
        s_annotation += ' -A MappingQualityRankSumTest'
        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_annotator_QualByDepth.html
        s_annotation += ' -A QualByDepth'
        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_annotator_RMSMappingQuality.html
        s_annotation += ' -A RMSMappingQuality'
        s_annotation += ' -A ReadPosRankSumTest'
        lines += [' %s \\' %(s_annotation)]

        lines += [' --emitRefConfidence GVCF \\']

        lines += ['\n']

        return lines


    def init_GATK_cmd(self,analysis_type,memMB,bool_checkpoint=False):

        s = ''
        ## Java version alert: starting with release 2.6, GATK now requires Java 1.7. See Version Highlights for 2.6 for details.
        ## http://www.broadinstitute.org/gatk/guide/article?id=2846
        s_java = self.init_java(
            self.fp_GATK, memMB, java=self.java, bool_checkpoint=bool_checkpoint)
        s += ' %s \\' %(s_java)
        lines = [s]

        ## CommandLineGATK, required, in
        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_CommandLineGATK.html#--analysis_type
        lines += [' --analysis_type %s \\' %(analysis_type)]
        ## CommandLineGATK, optional, in
        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_CommandLineGATK.html#--reference_sequence
        lines += [' --reference_sequence %s \\' %(self.reference_sequence)]
        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_CommandLineGATK.html#-nct
        if self.nct != 1:
            lines += [' --num_cpu_threads_per_data_thread %s \\' %(self.nct)]


        return lines


    def bsub_cmd(
        self,
        analysis_type,
        J,
        queue='normal',memMB=4000,
        LSF_affix=None,
        chrom=None,
        chromlen=None,
        index=None,
        bool_checkpoint=False,
        num_threads = None,
        bam=None,
        ):

        if not LSF_affix:
            LSF_affix = '%s/%s' %(analysis_type,analysis_type,)

        cmd = 'bsub -J"%s" -q %s' %(J,queue,)
        cmd += ' -G %s' %(self.project)
        cmd += " -M%i -R'select[mem>%i] rusage[mem=%i]'" %(
            memMB,memMB,memMB,)
        cmd += ' -o %s/LSF/%s.out' %(os.getcwd(), LSF_affix)
        cmd += ' -e %s/LSF/%s.err' %(os.getcwd(), LSF_affix)
        if num_threads:
            cmd += ' -n%i -R"span[hosts=1]"' %(num_threads)
        if bool_checkpoint:
            cmd += ' -k "%s method=blcr 710"' %(
                os.path.join(os.getcwd(),'checkpoint'))
            cmd += ' -r'
        if bool_checkpoint:
            cmd += ' cr_run'
        cmd += ' bash %s/shell/%s.sh' %(os.getcwd(),analysis_type,)
        if chrom:
            cmd += ' %s' %(chrom)
        if chromlen:
            cmd += ' %s' %(chromlen)
        if index:
            cmd += ' %s' %(index)
        if bam:
            cmd += ' %s' %(bam)

        return cmd


    def term_cmd(self, analysis_type, l_fp_out, extra=None):

        if type(l_fp_out) != list:
            print(l_fp_out)
            stop


        lines = ['\n']
        lines += ['if [ $? -eq 0 ]; then']
        for fp_out in l_fp_out:
            fp_out = fp_out
            ## previous command exited cleanly
            lines += ['if [ ! -s %s ]; then exit; fi' %(fp_out)]
            lines += ['echo %s >> touch/%s.touch' %(fp_out, analysis_type,)]

        if extra:
            lines += ['%s\n' %(extra)]

        lines += ['bash ./rerun.sh']
        lines += ['fi']

        ## write continuation shell script
        ## do not continue as part of previous command
        ## as this will influence CPU statistics
        ## and risk job of hitting CPU walltime
        s = "bsub -R 'select[mem>1500] rusage[mem=1500]' -M1500 \\\n"
        s += ' -o LSF/rerun.out \\\n'
        s += ' -e LSF/rerun.err \\\n'
        s += ' -G %s \\\n' %(self.project)
        s += ' bash ./rerun_python.sh'
        fd = open('rerun.sh','w')
        fd.write(s)
        fd.close()
        self.execmd('chmod +x rerun.sh')

        s = ' python'
        s += ' %s' %(sys.argv[0])
        for k,v in vars(self.namespace_args).items():
            s += ' --%s %s' %(k,str(v))
        fd = open('rerun_python.sh','w')
        fd.write(s)
        fd.close()
        self.execmd('chmod +x rerun_python.sh')

        return lines


    def is_file(self, str_):
        if not os.path.isfile(str_) and not os.path.islink(str_):
            msg = '%s is neither a readable file nor a symbolic link' % str_
            raise argparse.ArgumentTypeError(msg)
        return str_


    def is_file_or_dir(self, str_):
        print(str_)
        if not any([
            os.path.isfile(str_),os.path.islink(str_),os.path.isdir(str_)]):
            msg = '%s is neither a readable file nor a directory' % str_
            raise argparse.ArgumentTypeError(msg)
        return str_


    def add_arguments(self,parser):

        ## required arguments

        parser.add_argument(
            '--fp_bams','--bam','--bams','--input',
            help='Path to BAM and/or directory containing BAMs',
            nargs='+', required=True, type=self.is_file_or_dir)

        parser.add_argument('--coverage', required=True, type=float)

        parser.add_argument(
            '--fp_GATK', '--GATK', '--gatk', '--jar', required = True,
            help='File path to GATK',)

        parser.add_argument('--project', required=True)

        parser.add_argument('--arguments','--args')

        parser.add_argument('--java', required=True)

        ##
        ## CommandLineGATK arguments
        ##

        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_CommandLineGATK.html#--intervals
        parser.add_argument(
            '--intervals','-L',
            help='Additionally, one may specify a rod file to traverse over the positions for which there is a record in the file (e.g. -L file.vcf).',
            )

        ## http://www.broadinstitute.org/gatk/guide/article?id=1975
        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_CommandLineGATK.html#-nct
        parser.add_argument(
            '--nct', '--num_cpu_threads_per_data_thread',
            type=int, default=1)

        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_CommandLineGATK.html#--reference_sequence
        parser.add_argument(
            '--reference_sequence', '-R', required=True, type=self.is_file)

        ##
        ## HaplotypeCaller specific arguments
        ##

        ## Optional Inputs
        
        ## https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_haplotypecaller_HaplotypeCaller.html#--alleles
        parser.add_argument(
            '--alleles', default='',
            help='The set of alleles at which to genotype when --genotyping_mode is GENOTYPE_GIVEN_ALLELES.',)

        ## https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_haplotypecaller_HaplotypeCaller.html#--dbsnp
        parser.add_argument('--dbsnp', '-D')

        ## Optional Parameters

        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_genotyper_UnifiedGenotyper.html#--genotyping_mode
        parser.add_argument(
            '--genotyping_mode','-gt_mode',
            help='Specifies how to determine the alternate alleles to use for genotyping.',
            default='DISCOVERY')

        ##
        ## CombineGVCFs/GenotypeGVCFs related arguments
        ##

        ## http://gatkforums.broadinstitute.org/discussion/4074/file-number-limit-for-genotypegvcfs
        parser.add_argument('--gVCF_limit', default=200, type=int)

        ##
        ## VariantRecalibrator resources
        ##

        parser.add_argument(
            '--fp_resources_SNP', '--resources_SNP', '--VR_snp',
            help='File path to a file with -resource lines to append to GATK VR',)

        parser.add_argument(
            '--fp_resources_INDEL', '--resources_INDEL', '--VR_indel',
            help='File path to a file with -resource lines to append to GATK VR',)

##        parser.add_argument(
##            '--hapmap',
##            help='File path to hapmap vcf to be used by VariantCalibrator (e.g. hapmap_3.3.b37.sites.vcf)',
##            )
##

##        parser.add_argument(
##            '--omni',
##            help='File path to omni vcf to be used by VariantCalibrator (e.g. 1000G_omni2.5.b37.sites.vcf)',
##            )

        ##
        ## ApplyRecalibration
        ##

        parser.add_argument(
            '--ts_filter_level','--ts', type=float, required=True,)

        ##
        ## BEAGLE
        ##
        parser.add_argument(
            '--fp_software_beagle', '--beagle','--BEAGLE','--BEAGLEjar',
            help='File path to BEAGLE .jar file (e.g. beagle_3.3.2.jar)',
            required=True,
            )

        parser.add_argument(
            '--i_BEAGLE_size',
            help='Size (Mbp) of divided parts.',
            type=int,default=2, ## CPU bound (2Mbp=22hrs,3090MB) with lowmem option...
            )

        parser.add_argument(
            '--i_BEAGLE_edge',
            help='Window size (kbp) at either side of the divided part to avoid edge effects.',
            type=int,default=150,
            )

        s_help = 'Phased file to be divided (e.g.'
        s_help += ' ALL.chr$CHROM.phase1_release_v3.20101123.filt.renamed.bgl'
        parser.add_argument(
            '--BEAGLE_phased', '--phased',
            help=s_help,
            required = True,
            )

        s_help = 'Markers file to be divided (e.g.'
        s_help += ' ALL.chr$CHROM.phase1_release_v3.20101123.filt.renamed.markers'
        parser.add_argument(
            '--BEAGLE_markers', '--markers',
            help=s_help, required=True,
            )

        parser.add_argument('--i_BEAGLE_nsamples', default=20)

        ##
        ## optional arguments
        ##
        parser.add_argument(
            '--checkpoint', dest='bool_checkpoint', action='store_true', default=False)

        parser.add_argument(
            '--chroms', type=str, nargs='+',
            default=[str(i+1) for i in range(22)]+['X','Y',])

        return parser


    def parse_arguments(self,):

        parser = argparse.ArgumentParser()

        parser = self.add_arguments(parser)

        ## parse arguments to argparse NameSpace
        self.namespace_args = namespace_args = parser.parse_args()

        ## setatrr
        for k,v in vars(namespace_args).items():
            setattr(self,k,v)

        if self.fp_GATK is None and self.fp_options is None:
            parser.error('--GATK or --arguments')

        s_arguments = ''
        for k,v in vars(namespace_args).items():
            s_arguments += '%s %s\n' %(k,v)

        if self.arguments == None or self.arguments == 'None':
            self.arguments = '%s.arguments' %(self.project)
            fd = open(self.arguments,'w')
            fd.write(s_arguments)
            fd.close()
        else:
            fd = open(self.arguments,'r')
            lines = fd.readlines()
            fd.close()
            for line in lines:
                l = line.strip().split()
                k = l[0]
                v = l[1]
                setattr(self,k,v)

        if self.ts_filter_level not in [None,'None',]:
            if self.ts_filter_level < 1:
                self.ts_filter_level *= 100.

        self.bams = []
        for fp in self.fp_bams:
            if os.path.isdir(fp):
                self.bams += glob.glob(os.path.join(fp,'*.bam'))
            elif os.path.isfile(fp):
                self.bams += [fp]
            else:
                stop_take_care_of_symlinks

        return


    def __init__(self,):

        ##
        ## parse command line arguments
        ##
        self.parse_arguments()

        self.verbose = True

        return


if __name__ == '__main__':
    self = main()
    self.main()

## TODO

## todo20140902: tc9: run HC per chromosome
## todo20140902: tc9: make CombineGVCFs start each time 200 samples and a chromosome finishes? different/random 200 samples or the same 200 samples for each chromosome..?!
## todo20140902: tc9: make it possible to add on extra samples (e.g. NA12878) and run GenotypeGVCFs after running HC on the new sample bam.

## DISK USAGE

## todo20140216: tc9: write all output (BEAGLE) to gzipped files


## CPU

## todo20130320: tc9: make the function BEAGLE_divide run in "parallel"; i.e. run a process for each chromosome
## todo20130424: tc9: run BEAGLE_divide in parallel for each chromosome


## MEMORY

## todo20130204: tc9: make memory sample size dependent... only tested on 3 datasets with 100 samples each...

## todo20130320: tc9: test memory requirements when more than 100 samples
