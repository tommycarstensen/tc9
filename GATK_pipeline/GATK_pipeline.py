#!/bin/python3

## Tommy Carstensen, Wellcome Trust Sanger Institute, 2012-2013

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

## README

## http://www.broadinstitute.org/gatk/about#typical-workflows
## http://www.broadinstitute.org/gatk/guide/topic?name=best-practices

class main():

    def main(self):

        ## define list of chromosomes
        l_chroms = [str(i) for i in range(1,22+1,)]+['X','Y',]

        self.init(l_chroms,)

        ## parse chromsome lengths from reference sequence
        d_chrom_lens = self.parse_chrom_lens()

        ##
        ## write shell scripts
        ##
        self.UnifiedGenotyper(l_chroms,d_chrom_lens,)

##        l_chroms = ['22'] ## tmp!!!
        self.VariantRecalibrator(l_chroms,d_chrom_lens,)

##        l_chroms = ['22','X','Y',] ## tmp!!!
        ## phased imputation reference panels not available for chrY
        l_chroms.remove('Y')
        ## Exception in thread "main" java.lang.ArrayIndexOutOfBoundsException: -2
        l_chroms.remove('X')
##        l_chroms = ['22'] ## tmp!!!
        self.BEAGLE(l_chroms,d_chrom_lens)

        self.unite_BEAGLE(l_chroms,d_chrom_lens,)

        return


    def execmd(self,cmd):

        print(cmd)
        os.system(cmd)

        return


    def init(self,l_chroms,):

        l_dn = self.mkdirs(l_chroms,)

        self.check_stderr(l_chroms,l_dn,)

        return


    def check_stderr(self,l_chroms,l_dn,):

        print('checking that stderr is empty')
        bool_exit = False
        for dn in l_dn:
            if dn[:3] != 'out': continue
            l_fn = os.listdir('LSF/%s' %(dn[4:],))
            for fn in l_fn:
                if fn[-4:] != '.err': continue ## use glob instead
                fp = 'LSF/%s/%s' %(dn[4:],fn)
                if os.path.getsize(fp) > 0:
                    print('error:', fp)
                    bool_exit = True
        if bool_exit == True: sys.exit(0)

##        ## it takes too long to search all the output files
##        ## instead do a tail on each and every file and check that
##        ## os.popen("tail -n19 stdout/*/*.out | head -n1").read().strip()
##        ## ==
##        ## "Successfully completed."
##        print 'checking in stdout that no jobs were terminated prematurely'
##        s = os.popen('fgrep TERM stdout/*/*.out').read().strip()
##        if s != '':
##            print s
##            bool_exit = True
##        if bool_exit == True: sys.exit(0)
##
##        ##
##        ## faster to just check for file output... this step is redundant...
##        ##
##        print 'checking in stdout that all jobs were successfully completed'
##        l = os.listdir('stdout')
##        for s in l:
##            if os.path.isdir('stdout/%s' %(s)):
##                for fn in os.listdir('stdout/%s' %(s)):
##                    fp = os.path.join('stdout',s,fn)
##                    self.check_out_line(fp)
##            ## elif os.path.isfile('stdout/%s' %(s)):
##            else:
##                continue
##                fp = os.path.join('stdout',s)
##                self.check_out_line(fp)

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

##        cmd = "awk '{print $1}' out_BEAGLE/%s/%s.*.like.r2" %(chrom,chrom)
        cmd = "awk 'FNR>1{print $1}' out_BEAGLE/%s/%s.*.like.gprobs" %(chrom,chrom)
        cmd += ' | sort -u > gprobs%s' %(chrom)
        self.execmd(cmd)

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
            if l[0] != chrom: continue
            index = int(l[1])
            pos1 = int(l[2])
            pos2 = int(l[3])
            d_index2pos[index] = [pos1,pos2,]
        fp_out = 'out_BEAGLE/%s.gprobs' %(chrom,)
        for index in range(1,max(d_index2pos.keys())+1,):
            fp_in = 'out_BEAGLE/%s/%s.%i.like.gprobs' %(chrom,chrom,index,)
            if not os.path.isfile(fp_out):
                cmd = 'head -n1 %s > %s' %(fp_in,fp_out,)
                self.execmd(cmd)
            if index % 10 == 0:
                print('BEAGLE_unite', chrom, index)
            if not os.path.isfile(fp_in):
                print('missing', fp_in)
                stop
            pos1 = d_index2pos[index][0]
            pos2 = d_index2pos[index][1]
            cmd = "awk 'NR>1{pos=int(substr($1,%i));" %(len(chrom)+2)
            cmd += "if(pos>%i&&pos<=%i) print $0}'" %(
                pos1,pos2,)
            cmd += ' %s >> %s' %(fp_in,fp_out,)
            self.execmd(cmd)

        return


    def unite_BEAGLE(self,l_chroms,d_chrom_lens,):

        ##
        ## 1) check input existence
        ##
        l_fp_in = self.BEAGLE_parse_output_files(l_chroms)
        bool_exit = self.check_in('BEAGLE',l_fp_in,)

        for chrom in l_chroms:
            print(chrom)
            fp = 'out_BEAGLE/%s.gprobs' %(chrom)
            if os.path.isfile(fp): continue
            self.BEAGLE_fileIO_checks(chrom)
            self.BEAGLE_unite(chrom)

        return


    def BEAGLE_parse_output_files(self,l_chroms,):

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
            if not chrom in l_chroms: continue
            index = int(l[1])
            fp_in = 'out_BEAGLE/%s/' %(chrom,)
            fp_in += '%s.%i.like.gprobs' %(chrom,index,)
            l_fp_in += [fp_in]

        return l_fp_in


    def write_shell(self,fp,lines,):

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


    def BEAGLE_divide(self,l_chroms,):

        d_indexes = {}
        for chrom in l_chroms:
            d_indexes[chrom] = {}

        ## parse the input tranches file describing where to cut the data
        fp_tranches = 'out_VariantRecalibrator/VariantRecalibrator.SNP.tranches'
        minVQSLOD = self.parse_minVQSLOD(fp_tranches,self.f_ts_filter_level,)

        ## open the recal file
        fp_recal = 'out_VariantRecalibrator/VariantRecalibrator.SNP.recal'
        with open(fp_recal,'r') as fd_recal:

            ## loop over the raw input variants to be recalibrated
            if os.path.islink('out_UnifiedGenotyper'):
                path = os.readlink('out_UnifiedGenotyper')
            else:
                path = 'out_UnifiedGenotyper'

            for chrom in l_chroms:
                l_files = glob.glob('%s/%s.*.vcf' %(path,chrom,))
                if len(l_files) == 0:
                    print('no files found')
                    print('%s/%s.*.vcf' %(path,chrom,))
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
            l_vcf = line_vcf.split()
            bool_indel,bool_diallelic = self.parse_variant_type(l_vcf)
            if bool_indel == True:
                continue

            ## skip LowQual
            if l_vcf[6] == 'LowQual': continue

##            ## skip QUAL=.
##            if l_vcf[5] == '.': continue

            ## read equivalent VR line
            while True:
                line_recal = fd_recal.readline()
                if line_recal[0] != '#':
                    break

            ## skip non-diallelic SNPs
            if bool_diallelic == False:
                continue

            ## skip if VQSLOD value below threshold
            VQSLOD, l_recal = self.parse_VQSLOD(line_recal)
            ## tmp!!! check. remove after testing.
            if l_recal[0] != l_vcf[0] or l_recal[1] != l_vcf[1]:
                print('VRrecal')
                print(l_recal)
                print('UGvcf')
                print(line_vcf)
                print(l_recal[1],l_vcf[1])
                stoptmp
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

        with open(file) as fd_vcf:
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
        fp_phased = self.fp_BEAGLE_phased.replace('$CHROMOSOME',chrom)
        fp_markers = self.fp_BEAGLE_markers.replace('$CHROMOSOME',chrom)
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

        index = 1 ## LSF does not allow 0 for LSB_JOBINDEX! "Bad job name. Job not submitted." ## e.g. bsub -J"test[0-2]" echo A
##        pos_ref = 0
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
        with fileinput.input(files=l_files) as fd_vcf:
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
                    pos_ref, A_ref, B_ref = self.parse_marker(line_m)

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
    
        return d_index2pos


    def vcf2beagle(self,l_vcf):

        ## append chrom:pos alleleA alleleB
        line_beagle = '%s:%s %s %s' %(l_vcf[0],l_vcf[1],l_vcf[3],l_vcf[4],)

        for s_vcf in l_vcf[9:]:
            ## variant not called
            if s_vcf == './.':
                line_beagle += ' 0.3333 0.3333 0.3333'
                continue
            l_probs = []
            l_log10likelihoods = s_vcf.split(':')[-1].split(',')
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
##grep "100022465\|100043981\|100047620" /lustre/scratch113/projects/agv/users/tc9/HiSeq/pipeline/uganda4x/out_UnifiedGenotyper/UnifiedGenotyper.10.*.vcf
##UnifiedGenotyper.10.11.vcf:10    100022465       rs78402300      C       A,T     384.41  .       AC=5,4;AF=0.025,0.020;AN=198;BaseQRankSum=-1.132;DB;DP=459;Dels=0.00;FS=0.000;HaplotypeScore=0.2870;InbreedingCoeff=0.2595;MLEAC=4,3;MLEAF=0.020,0.015;MQ=57.71;MQ0=0;MQRankSum=1.004;QD=13.73;ReadPosRankSum=-0.013;SB=-1.944e+02      GT:AD:DP:GQ:PL  0/0:5,0,0:5:9:0,9,120,9,120,120 0/2:2,0,5:7:55:165,171,241,0,70,55      0/0:5,0,0:5:15:0,15,186,15,186,186      0/0:5,0,0:5:15:0,15,191,15,191,191      0/0:4,0,0:4:12:0,12,147,12,147,147      0/0:5,0,0:5:15:0,15,175,15,175,175      0/0:6,0,0:6:18:0,18,238,18,238,238      0/0:3,0,0:3:9:0,9,123,9,123,123 0/0:3,
##UnifiedGenotyper.10.11.vcf:10    100043981       rs4919198       T       A,G     14156.03        .       AC=40,119;AF=0.202,0.601;AN=198;BaseQRankSum=-3.709;DB;DP=530;Dels=0.00;FS=0.427;HaplotypeScore=0.4226;InbreedingCoeff=0.4117;MLEAC=34,129;MLEAF=0.172,0.652;MQ=58.46;MQ0=0;MQRankSum=0.574;QD=28.37;ReadPosRankSum=0.935;SB=-6.345e+03 GT:AD:DP:GQ:PL  2/2:0,0,4:4:12:158,158,158,12,12,0      0/2:8,0,1:9:14:14,38,335,0,297,294      1/2:0,4,3:7:94:258,106,94,152,0,143    0/2:6,0,2:8:48:48,66,287,0,221,215       0/2:2,0,2:4:58:58,64,133,0,70,64        2/2:0,0,9:9:27:346,346,346,27,27,0      2/2:0,0,3:3:9:90,90,90,9,9,0    2/2:0,0,10:10:30:360,360,360,3
##UnifiedGenotyper.10.11.vcf:10    100047620       rs74779652      C       A,G     1442.06 .       AC=2,28;AF=0.010,0.141;AN=198;BaseQRankSum=-6.605;DB;DP=555;Dels=0.00;FS=1.046;HaplotypeScore=0.4547;InbreedingCoeff=0.1284;MLEAC=2,25;MLEAF=0.010,0.126;MQ=57.69;MQ0=0;MQRankSum=0.864;QD=9.81;ReadPosRankSum=-0.200;SB=-8.715e+02     GT:AD:DP:GQ:PL  0/0:8,0,0:8:24:0,24,281,24,281,281      0/0:12,0,0:12:36:0,36,497,36,497,497    0/0:10,0,0:10:27:0,27,342,27,342,342    0/0:6,0,0:6:18:0,18,242,18,242,242      0/0:1,0,0:1:3:0,3,42,3,42,42    0/1:4,4,0:8:99:131,0,102,140,114,255    0/0:6,0,0:6:18:0,18,216,18,216,216      0/2:6,0,2:8:42:42,57,236,0,179,173    

##grep "100022465\|100043981\|100047620" /lustre/scratch113/projects/agv/users/tc9/HiSeq/pipeline/uganda4x/out_BEAGLE/10.gprobs | cut -d " " -f-10
##10:100022465 C A 1 0 0 1 0 0 1
##10:100043981 T G 0 0.0056 0.9944 0 0.9954 0.0046 0
##10:100047620 C G 1 0 0 1 0 0 1

##grep "100022465\|100043981\|100047620" /nfs/team149/resources/beagle_1000_Genomes.phase1_release_v3/ALL.chr10.phase1_release_v3.20101123.filt.markers
##rs78402300      100022465       C       A
##rs4919198       100043981       T       G
##rs74779652      100047620       C       G
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

        cmd = 'cat %s' %(fp_phased)
        cmd += "| awk 'NR>1{print $2}' | sort -u > panel0in%s" %(chrom)
        self.execmd(cmd)

        cmd = "awk 'FNR>1{print $2}' in_BEAGLE/%s/%s.phased" %(chrom,chrom)
        cmd += "| sort -u > panel0out%s" %(chrom)
        self.execmd(cmd)

        minVQSLOD = self.parse_minVQSLOD(fp_tranches,self.f_ts_filter_level,)
        fp_recal = 'out_VariantRecalibrator/VariantRecalibrator.SNP.recal'
        with open(fp_recal,'r') as fd_recal, open('panel2in%s' %(chrom),'w') as fd_out:
            for line in fd_recal:
                if line[0] == '#': continue
                VQSLOD, l_recal = self.parse_VQSLOD(line)
                if VQSLOD < minVQSLOD:
                    continue
                ## slow conversion to integer just to make sure it is an integer
                pos = int(l_recal[1])
                fd_out.write('%s:%i\n' %(chrom,pos))

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


    def BEAGLE(self,l_chroms,d_chrom_lens,):

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
        memMB = 3900
        ## fgrep CPU */stdout/BEAGLE/*.out | sort -k5nr,5 | head -n1
        ## zulu_20121208/stdout/BEAGLE/BEAGLE.19.1.out:    CPU time   :  80270.55 sec.
        queue = 'long'

        ##
        ## 1) check input existence
        ##
        mode = 'SNP'
        fp_in_recal = 'out_VariantRecalibrator/VariantRecalibrator.%s.recal' %(mode)
        fp_in_tranches = 'out_VariantRecalibrator/VariantRecalibrator.%s.tranches' %(mode)
        bool_exit = self.check_in('VariantRecalibrator',[fp_in_recal,fp_in_tranches,],)

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
            d_indexes = self.BEAGLE_divide(l_chroms)
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
            for chrom in l_chroms:
                d_indexes[chrom] = {}
            fd = open('BEAGLE_divide_indexes.txt','r')
            lines = fd.readlines()
            fd.close()
            for line in lines:
                l = line.strip().split()
                chrom = l[0]
                if not chrom in l_chroms: continue
                index = int(l[1])
                pos1 = int(l[2])
                pos2 = int(l[3])
                d_indexes[chrom][index] = [pos1,pos2,]

        ##
        ## execute shell script
        ##
        for chrom in l_chroms:

            print('bsub BEAGLE %s' %(chrom))

##            J = '%s%s[%i-%i]' %('BEAGLE',chrom,1,max(d_indexes[chrom].keys()),)
            for index in d_indexes[chrom].keys():
                fn_out = 'out_BEAGLE/%s/%s.%i.like.gprobs' %(chrom,chrom,index)
                if os.path.isfile(fn_out):
                    continue
                J = '%s.%s[%s-%s]' %('BEAGLE',chrom,index,index,)
                std_suffix = '%s/%s.%%I' %('BEAGLE',chrom)
                cmd = self.bsub_cmd(
                    'BEAGLE',J,memMB=memMB,std_suffix=std_suffix,chrom=chrom,
                    queue=queue,)
                os.system(cmd)

        return


    def init_java(self, jar, memMB, java='java', bool_checkpoint=False):

        s = '%s -Djava.io.tmpdir=%s' %(java,'tmp')
        ## set maximum heap size
        s += ' -Xmx%im' %(memMB)
        if bool_checkpoint:
            s += ' -XX:-UsePerfData -Xrs '
        s += ' -jar %s' %(jar)

        return s


    def BEAGLE_write_shell_script(self,memMB,):

        fp_out = self.d_out['BEAGLE']

        ## initiate shell script
        lines = ['#!/bin/bash\n']

        ## parse chromosome from command line
        lines += ['CHROMOSOME=$1\n']

        lines += ['if [ -s %s.gprobs ]; then\nexit\nfi\n' %(fp_out)] ## redundant

##        ## init cmd
##        lines += ['cmd="']

        ##
        ## initiate BEAGLE
        ##
        s_java = self.init_java(self.fp_software_beagle,memMB)
        lines += ['%s \\' %(s_java)]

        lines += self.body_BEAGLE(fp_out,)

        ## term cmd
        lines += self.term_cmd('BEAGLE',['%s.gprobs' %(fp_out)],)

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
##        s += ' markers=/lustre/scratch107/projects/uganda/users/tc9/in_BEAGLE/ALL.chr%s.phase1_release_v3.20101123.filt.markers ' %(chrom)
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
        fn = '%s.fai' %(self.fp_FASTA_reference_sequence)
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


    def mkdirs(self,l_chroms,):

        l_dn = [
            'touch',
            'LSF',
            'shell',
            'out_UnifiedGenotyper',
            'out_VariantRecalibrator',
            'in_BEAGLE','out_BEAGLE',
            ]

        ## create subdirs
        for dn in l_dn:
            if not os.path.isdir(dn) and not os.path.islink(dn):
                os.mkdir(dn)
            if dn != 'touch' and (dn[:4] == 'out_' or dn[:3] == 'in_'):
                if not os.path.isdir(os.path.join('touch',dn)):
                    os.mkdir(os.path.join('touch',dn))
        for dn in l_dn:
            if dn[:3] != 'out': continue
            if not os.path.isdir('LSF/%s' %(dn[4:],)):
                os.mkdir('LSF/%s' %(dn[4:],))

        for chrom in l_chroms:
            for dn in ['in_BEAGLE','out_BEAGLE',]:
                if not os.path.isdir('%s/%s' %(dn,chrom)):
                    os.mkdir('%s/%s' %(dn,chrom))
                if not os.path.isdir('touch/%s/%s' %(dn,chrom)):
                    os.mkdir('touch/%s/%s' %(dn,chrom))

        ## java (i.e. BEAGLE) tmp dir
        if not os.path.isdir('tmp'):
            os.mkdir('tmp')

        return l_dn


    def check_in(self,analysis_type,l_fp_in,):
        
##        fd = open('%s.touch' %(analysis_type),'r')
##        s = fd.read()
##        fd.close()
##        l_fp_out = s.split('\n')

        d_l_fp_out = {}
        for dirname in ['','touch',]:
            d_l_fp_out[dirname] = []
            l = os.listdir(os.path.join(dirname,'out_%s' %(analysis_type)))
            for s in l:
                path1 = os.path.join('out_%s' %(analysis_type),s)
                path2 = os.path.join(dirname,path1)
                ## append files in chromosomal subdirectories (e.g. BEAGLE)
                if os.path.isdir(path2):
                    l = os.listdir(path2)
                    for fn in l:
                        d_l_fp_out[dirname] += [os.path.join(path1,fn)]
                ## append files in main dir (e.g. UnifiedGenotyper)
                elif os.path.isfile(path2):
                    d_l_fp_out[dirname] += [path1]

        bool_exit = False
        for dirname,l_fp_out in d_l_fp_out.items():
            if len(set(l_fp_in)-set(l_fp_out)) > 0:
                print('%s and possibly %s other files not generated.' %(
                    list(set(l_fp_in)-set(l_fp_out))[0],
                    len(set(l_fp_in)-set(l_fp_out))-1,))
                print('dirname', dirname)
##                for fp in list(set(l_fp_in)-set(l_fp_out)):
##                    os.system('touch touch/%s' %(fp))
##                    os.system('touch %s' %(fp))
##                if dirname == 'touch':
##                    for fp in list(set(l_fp_in)-set(l_fp_out)):
##                        if os.path.isfile('%s' %(fp)):
##                            os.system('touch touch/%s' %(fp))
##                stop
                print('%s has not run to completion. Exiting.' %(analysis_type))
                bool_exit = True
#                print(inspect.stack()[1])
                sys.exit()

        return bool_exit


    def determine_TS_level(self,fp_tranches,):

        cmd = 'ts_filter_level=$('
        ## loop over lines
        cmd += 'cat %s' %(fp_tranches)
        ## init awk
        cmd += " | awk '"
        ## BEGIN
        cmd += ' BEGIN{FS=","}'
        ## skip header and loop over rows
        cmd += ' NR>3{'
        ## if below novelTiTv ratio
        cmd += ' if($5<2.1) {'
        ## if first line
        cmd += ' if(NR==4) {print $1; exit}'
        ## else if not first line
        cmd += ' else {'
        cmd += ' novelTiTv2=$5; targetTruthSensitivity2=$1;'
        cmd += ' dy=(novelTiTv2-novelTiTv1);'
        cmd += ' dx=(targetTruthSensitivity2-targetTruthSensitivity1);'
        cmd += ' slope=dy/dx;'
        cmd += ' intercept=novelTiTv1-slope*targetTruthSensitivity1;'
        cmd += ' print (2.1-intercept)/slope; exit}'
        ## else if above novelTiTv ratio
        cmd += ' } else {targetTruthSensitivity1=$1;novelTiTv1=$5}'
        ## continue loop over rows
        cmd += ' }'
        ## term awk
        cmd += " '"
        cmd += ')'

        return cmd


    def get_fps_in(self,l_chroms,d_chrom_lens,bps_per_interval,):

        l_vcfs_in = []
        for chrom in l_chroms:
            intervals = int(math.ceil(
                d_chrom_lens[chrom]/float(bps_per_interval)))
            for interval in range(1,intervals+1,):
                fp_in = 'out_%s/%s.%i.vcf' %(
                    'UnifiedGenotyper',chrom,interval,)
##                if os.path.isfile(fp_in):
##                    if os.path.getsize(fp_in) == 0:
##                        continue
                l_vcfs_in += [fp_in]

        return l_vcfs_in


    def VariantRecalibrator(self,l_chroms,d_chrom_lens):

        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_variantrecalibration_VariantRecalibrator.html

        T = analysis_type = 'VariantRecalibrator'
        ##fula_20120704/stdout/VariantRecalibrator/VariantRecalibrator.out:    Max Memory :     12840 MB
        ##zulu_20121208/stdout/VariantRecalibrator/VariantRecalibrator.out:    Max Memory :     13710 MB
        ##uganda_20130113/stdout/VariantRecalibrator/VariantRecalibrator.out:    Max Memory :     13950 MB
        ##pipeline/ethiopia4x/stdout/VariantRecalibrator/VariantRecalibrator.out:    Max Memory :     15092 MB
        memMB = 15900
##        memMB = 27500 ## tmp!!! union4x
        ##fula_20120704/stdout/VariantRecalibrator/VariantRecalibrator.out:    CPU time   :   9635.38 sec.
        ##zulu_20121208/stdout/VariantRecalibrator/VariantRecalibrator.out:    CPU time   :  10015.61 sec.
        ##uganda_20130113/stdout/VariantRecalibrator/VariantRecalibrator.out:    CPU time   :  11657.10 sec.
        ##pipeline/zulu1x/stdout/VariantRecalibrator/VariantRecalibrator.out:    CPU time   :  11309.62 sec.
        queue = 'normal'
        queue = 'yesterday'

        if os.path.isfile('VariantRecalibrator.touch'): return

        ##
        ## 1) check input existence (vcf)
        ##
        l_vcfs_in = self.get_fps_in(
            l_chroms,d_chrom_lens,self.i_UG_or_HC_size,)
        bool_exit = self.check_in('UnifiedGenotyper',l_vcfs_in,)

        ##
        ## 1b) check input size (vcf.idx)
        ##
        bool_exit = False
        for vcf_in in l_vcfs_in:
            fp_vcf_idx = '%s.idx' %(vcf_in)
            if not os.path.isfile(fp_vcf_idx):
                continue
            if os.path.getsize(fp_vcf_idx) == 0:
                print('zero size:', fp_vcf_idx)
                bool_exit = True
        if bool_exit == True: sys.exit(0)

        d_resources = {'SNP':self.fp_resources_SNP,'INDEL':self.fp_resources_INDEL,}

        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_variantrecalibration_VariantRecalibrator.html#--mode
        for mode in ['SNP','INDEL',]:

            ##
            ## 2) touch
            ##
            bool_return = self.touch('%s.%s' %(analysis_type,mode))
            if bool_return == True: continue

            fp_tranches = 'out_VariantRecalibrator/VariantRecalibrator.%s.tranches' %(mode)
            fp_recal = 'out_VariantRecalibrator/VariantRecalibrator.%s.recal' %(mode)

            ##
            ## init GATK walker
            ##
            lines = self.init_GATK_cmd(analysis_type,memMB,)

            ## GATKwalker, required, in
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
            ## GATKwalker, required, out
            ##
            lines += [' --recal_file %s \\' %(fp_recal)]
            lines += [' --tranches_file %s \\' %(fp_tranches)]

            ##
            ## GATKwalker, optional, in
            ##
            
            ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_variantrecalibration_VariantRecalibrator.html#--mode
            lines += [' --mode %s \\' %(mode)]

            ## http://gatkforums.broadinstitute.org/discussion/1259/what-vqsr-training-sets-arguments-should-i-use-for-my-specific-project
            if mode == 'INDEL':
                ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_variantrecalibration_VariantRecalibrator.html#--maxGaussians
                lines += [' --maxGaussians 4 \\'] ## default 10
                ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_variantrecalibration_VariantRecalibrator.html#--percentBadVariants
                lines += [' --percentBadVariants 0.01 \\'] ## default 0.03
                ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_variantrecalibration_VariantRecalibrator.html#--minNumBadVariants
                lines += [' --minNumBadVariants 1000 \\'] ## default 2500

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
            l_TStranches += [70+i for i in range(19,-1,-1,)]
            s_TStranches = ''
            for TStranche in l_TStranches:
                s_TStranches += '--TStranche %.1f ' %(TStranche)
            lines += [' %s \\' %(s_TStranches)]

            ##
            ## GATKwalker, optional, out
            ##
            lines += [' --rscript_file out_%s/%s.%s.plots.R \\' %(T,T,mode,)]

            ##
            ## term GATK walker
            ##
            lines += self.term_cmd(
                '%s.%s' %(analysis_type,mode),[fp_tranches,fp_recal,],)

            self.write_shell('shell/%s.%s.sh' %(analysis_type,mode,),lines,)

            J = 'VR'
            cmd = self.bsub_cmd(
                '%s.%s' %(analysis_type,mode),J,memMB=memMB,queue=queue,
                std_suffix='%s/%s' %(analysis_type,mode,),)
            self.execmd(cmd)

        return


    def init_GATK_cmd(self,analysis_type,memMB,bool_checkpoint=False):

##        s = 'cmd="'
        s = ''
        ## Java version alert: starting with release 2.6, GATK now requires Java 1.7. See Version Highlights for 2.6 for details.
        ## http://www.broadinstitute.org/gatk/guide/article?id=2846
        if '2.7' in self.fp_GATK or '2.6' in self.fp_GATK or '2.8' in self.fp_GATK or 'nc6' in self.fp_GATK:
            java = '/software/team149/opt/jre1.7.0_45/bin/java'
        elif '2.5' in self.fp_GATK:
            java = 'java'
        else:
            print('unknown GATK version for %s' %(analysis_type))
            print(self.fp_GATK)
            sys.exit()
        s_java = self.init_java(
            self.fp_GATK, memMB, java=java, bool_checkpoint=bool_checkpoint)
        s += ' %s \\' %(s_java)
        lines = [s]

        ## CommandLineGATK, required, in
        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_CommandLineGATK.html#--analysis_type
        lines += [' --analysis_type %s \\' %(analysis_type)]
        ## CommandLineGATK, optional, in
        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_CommandLineGATK.html#--reference_sequence
        lines += [' --reference_sequence %s \\' %(self.fp_FASTA_reference_sequence)]

        return lines


    def touch(self,analysis_type):

        bool_return = False
        fn_touch = '%s.touch' %(analysis_type)
        if os.path.isfile(fn_touch):
                if self.verbose == True:
                    print('in progress or completed:', analysis_type)
                bool_return = True
        else:
            self.execmd('touch %s' %(fn_touch))

        return bool_return


    def UnifiedGenotyper(self,l_chroms,d_chrom_lens,):

        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_genotyper_UnifiedGenotyper.html

        T = analysis_type = 'UnifiedGenotyper'
        ## fula_20120704/stdout/UnifiedGenotyper.16.9.out:    CPU time   :   9420.62 sec.
        ## zulu_20121208/stdout/UnifiedGenotyper/UnifiedGenotyper.10.12.out:    CPU time   :  12106.00 sec.
##        queue = 'normal'
        ## ../ethiopia8x/stdout/UnifiedGenotyper/UnifiedGenotyper.1.7.out:    CPU time   :  63554.28 sec.
        queue = 'long'
        ## zulu_20121208/stdout/UnifiedGenotyper/UnifiedGenotyper.16.5.out:    Max Memory :      2030 MB
        if self.fp_intervals:
            memMB = 5900
        elif self.project == 'ug2g':
            ## LSF/UnifiedGenotyper/UnifiedGenotyper.22.3.out:    Max Memory :             12075 MB
            ## LSF/UnifiedGenotyper/UnifiedGenotyper.22.6.out:    Max Memory :             8817 MB
            memMB = 12900 ## tmp!!!
            queue = 'yesterday'
            queue = 'long'
            queue = 'basement'
            queue = 'small'
            queue = 'normal'
        else:
            memMB = 2900 ## 2600 if 250 samples at 4x (helic), 2100 if 120 samples at 4x or 8x (ethiopia) ## 1900=2000-default
            memMB = 3900

        ##
        ## 1) touch
        ##
        bool_return = self.touch(analysis_type)
        if bool_return == True: return

        fp_out = self.d_out[analysis_type]

        ## initiate shell script
        lines = ['#!/bin/bash']

        ## commands prior to GATK command
        lines += self.init_UG_or_HC(
            l_chroms,d_chrom_lens,analysis_type)

        ## initiate GATK command
        lines += self.init_GATK_cmd(
            'UnifiedGenotyper',memMB,bool_checkpoint=self.bool_checkpoint,)

        ## append GATK command options
        lines += self.body_UnifiedGenotyper(fp_out,)

        ## terminate shell script
        lines += self.term_cmd(analysis_type,[fp_out],)

        ## write shell script
        self.write_shell('shell/UnifiedGenotyper.sh',lines,)

##        ## execute shell script
##        for chrom in l_chroms:
##            intervals = int(math.ceil(
##                d_chrom_lens[chrom]/float(self.i_UG_or_HC_size)))
##            J = '%s%s[%i-%i]' %('UG',chrom,1,intervals,)
##            std_suffix = '%s/%s.%s.%%I' %(T,T,chrom)
##            cmd = self.bsub_cmd(
##                analysis_type,J,std_suffix=std_suffix,chrom=chrom,
##                memMB=memMB,queue=queue,
##                bool_checkpoint=bool_checkpoint)
##            self.execmd(cmd)

        ## execute shell script
        for chrom in l_chroms:
##            if int(chrom) != 21: continue ## tmp!!!
            for i in range(
                1,1+math.ceil(d_chrom_lens[chrom]/float(self.i_UG_or_HC_size))):
                J = '%s%s.%i' %('UG',chrom,i,)
                std_suffix = '%s/%s.%i' %(T,chrom,i)
                cmd = self.bsub_cmd(
                    analysis_type,J,std_suffix=std_suffix,chrom=chrom,
                    memMB=memMB,queue=queue,
                    bool_checkpoint=self.bool_checkpoint, index=i)
                self.execmd(cmd)

        return


    def HaplotypeCaller(self,l_chroms,d_chrom_lens,):

        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_genotyper_UnifiedGenotyper.html

        T = analysis_type = 'HaplotypeCaller'
        ## fula_20120704/stdout/UnifiedGenotyper.16.9.out:    CPU time   :   9420.62 sec.
        ## zulu_20121208/stdout/UnifiedGenotyper/UnifiedGenotyper.10.12.out:    CPU time   :  12106.00 sec.
##        queue = 'normal'
        ## ../ethiopia8x/stdout/UnifiedGenotyper/UnifiedGenotyper.1.7.out:    CPU time   :  63554.28 sec.
        queue = 'long'
        ## zulu_20121208/stdout/UnifiedGenotyper/UnifiedGenotyper.16.5.out:    Max Memory :      2030 MB
        if self.fp_intervals:
            memMB = 5900
        else:
            memMB = 2900 ## 2600 if 250 samples at 4x (helic), 2100 if 120 samples at 4x or 8x (ethiopia) ## 1900=2000-default
            memMB = 3900

        ##
        ## 1) touch
        ##
        bool_return = self.touch(analysis_type)
        if bool_return == True: return

        fp_out = self.d_out[analysis_type]

        ## initiate shell script
        lines = ['#!/bin/bash']

        ## commands prior to GATK command
        lines += self.init_UG_or_HC(
            l_chroms,d_chrom_le1ns,analysis_type)

        ## initiate GATK command
        lines += self.init_GATK_cmd(
            analysis_type,memMB,bool_checkpoint=self.bool_checkpoint)

        ## append GATK command options
        lines += self.body_HaplotypeCaller(fp_out,)

        ## terminate shell script
        lines += self.term_cmd(analysis_type,[fp_out],)

        ## write shell script
        self.write_shell('shell/%s.sh' %(analysis_type),lines,)

        ## execute shell script
        for chrom in l_chroms:
            intervals = int(math.ceil(
                d_chrom_lens[chrom]/float(self.i_UG_or_HC_size)))
            J = '%s%s[%i-%i]' %('UG',chrom,1,intervals,)
            std_suffix = '%s/%s.%%I' %(T,chrom)
            cmd = self.bsub_cmd(
                analysis_type,J,std_suffix=std_suffix,chrom=chrom,
                memMB=memMB,queue=queue,)
            self.execmd(cmd)

        return


    def bsub_cmd(
        self,
        analysis_type,
        J,
        queue='normal',memMB=4000,
        std_suffix=None,
        chrom=None,
        chromlen=None,
        index=None,
        bool_checkpoint=False,
        ):

        if not std_suffix:
            std_suffix = '%s/%s' %(analysis_type,analysis_type,)

        cmd = 'bsub -J"%s" -q %s' %(J,queue,)
        cmd += ' -G %s' %(self.project)
        cmd += " -M%i -R'select[mem>%i] rusage[mem=%i]'" %(
            memMB,memMB,memMB,)
        cmd += ' -o %s/LSF/%s.out' %(os.getcwd(),std_suffix)
        cmd += ' -e %s/LSF/%s.err' %(os.getcwd(),std_suffix)
        if bool_checkpoint:
            cmd += ' -k "%s method=blcr 710"' %(
                os.path.join(os.getcwd(),'checkpoint'))
            cmd += ' -r'
        if bool_checkpoint:
            cmd += ' cr_run'
        cmd += ' bash %s/shell/%s.sh' %(os.getcwd(),analysis_type,)
        if chrom != None:
            cmd += ' %s' %(chrom)
        if chromlen != None:
            cmd += ' %s' %(chromlen)
        if index != None:
            cmd += ' %s' %(index)

        return cmd


    def body_UnifiedGenotyper(self,fp_out,):

        lines = []

        ##
        ## required
        ##
        ## File to which variants should be written
        lines += [' --out %s \\' %(fp_out)]

        ##
        ## CommandLineGATK, optional
        ##

        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_CommandLineGATK.html#--intervals
        lines += [
            ' --intervals $CHROMOSOME:$(((${LSB_JOBINDEX}-1)*%i+1))-$posmax \\' %(
                self.i_UG_or_HC_size)]
        if self.fp_intervals:
            lines += [' --intervals %s \\' %(self.fp_intervals)]
            ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_CommandLineGATK.html#--interval_set_rule
            lines += [' --interval_set_rule INTERSECTION \\']

        ##
        ## UnifiedGenotyper, optional
        ##
        ## dbSNP file. rsIDs from this file are used to populate the ID column of the output.
        lines += [' --dbsnp %s \\' %(self.fp_vcf_dbsnp)]
        ## Selecting an appropriate quality score threshold
        ## A common question is the confidence score threshold
        ## to use for variant detection. We recommend:
        ## Deep (> 10x coverage per sample) data:
        ## we recommend a minimum confidence score threshold of Q30.
        ## Shallow (< 10x coverage per sample) data:
        ## because variants have by necessity lower quality
        ## with shallower coverage we recommend
        ## a minimum confidence score of Q4 in projects
        ## with 100 samples or fewer
        ## and Q10 otherwise.
        if self.project == 'ug2g':
            lines += [' -stand_call_conf 10 \\']
            lines += [' -stand_emit_conf 10 \\']
        else:
            print(self.project)
            stop
            lines += [' -stand_call_conf 4 \\']
            lines += [' -stand_emit_conf 4 \\']
        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_genotyper_UnifiedGenotyper.html#--genotype_likelihoods_model
        lines += [' --genotype_likelihoods_model SNP \\'] ## default value SNP
        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_genotyper_UnifiedGenotyper.html#--genotyping_mode
        lines += [' -gt_mode %s \\' %(self.genotyping_mode)] ## default value DISCOVERY
        lines += [' -out_mode %s \\' %(self.output_mode)] ## default value EMIT_VARIANTS_ONLY
        if self.fp_alleles:
            if os.path.isdir(self.fp_alleles):
                s = '%s/$CHROMOSOME.$LSB_JOBINDEX.vcf' %(self.fp_alleles)
                lines += [' --alleles %s \\' %(s)]
            elif os.path.isfile(self.fp_alleles):
                lines += [' --alleles %s \\' %(self.fp_alleles)]
            else:
                print(self.fp_alleles, 'not found')
                sys.exit()

        ## http://www.broadinstitute.org/gatk/gatkdocs/#VariantAnnotatorannotations
        s_annotation = ''

        ## http://gatkforums.broadinstitute.org/discussion/2318/undocumented-change-in-2-4-a-depthofcoverage
        if '2.4' in self.fp_GATK or '2.5' in self.fp_GATK or '2.7' in self.fp_GATK or '2.8' in self.fp_GATK or 'nc6' in self.fp_GATK:
            ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_annotator_Coverage.html
            s_annotation += ' --annotation Coverage'
        elif (
            '2.3' in self.fp_GATK
            or
            '2.2' in self.fp_GATK
            or
            '2.1' in self.fp_GATK
            ):
            ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_annotator_DepthOfCoverage.html
            s_annotation += ' --annotation DepthOfCoverage'
        else:
            print(self.fp_GATK)
            print('Unknown version of GATK. Please see why the version number must be known here:')
            print('http://gatkforums.broadinstitute.org/discussion/2318/undocumented-change-in-2-4-a-depthofcoverage')
            sys.exit()

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
        lines += [s_annotation+' \\']

        ##
        ## append input files
        ##
        for fp_bam in self.fp_bams:
            for s in os.listdir(fp_bam):
                if s[-4:] != '.bam':
                    continue
                ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_CommandLineGATK.html#--input_file
                path = os.path.join(fp_bam,s,)
                lines += [' --input_file %s \\' %(path)]

        return lines


    def body_HaplotypeCaller(self,fp_out,):

        lines = []

        ##
        ## required
        ##
        ## File to which variants should be written
        lines += [' --out %s \\' %(fp_out)]

        ##
        ## CommandLineGATK, optional
        ##

        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_CommandLineGATK.html#--intervals
        lines += [
            '--intervals $CHROMOSOME:$(((${LSB_JOBINDEX}-1)*%i+1))-$posmax \\' %(
                self.i_UG_or_HC_size)]
        if self.fp_intervals:
            lines += ['--intervals %s \\' %(self.fp_intervals)]
            ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_CommandLineGATK.html#--interval_set_rule
            lines += ['--interval_set_rule INTERSECTION \\']

        ##
        ## HaplotypeCaller, optional
        ##
        ## dbSNP file. rsIDs from this file are used to populate the ID column of the output.
        lines += [' --dbsnp %s \\' %(self.fp_vcf_dbsnp)]
        ## Selecting an appropriate quality score threshold
        ## A common question is the confidence score threshold
        ## to use for variant detection. We recommend:
        ## Deep (> 10x coverage per sample) data:
        ## we recommend a minimum confidence score threshold of Q30.
        ## Shallow (< 10x coverage per sample) data:
        ## because variants have by necessity lower quality
        ## with shallower coverage we recommend
        ## a minimum confidence score of Q4 in projects
        ## with 100 samples or fewer
        ## and Q10 otherwise.
        lines += [' -stand_call_conf 4 \\']
        lines += [' -stand_emit_conf 4 \\']
        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_haplotypecaller_HaplotypeCaller.html#--genotyping_mode
        lines += [' -gt_mode %s \\' %(self.genotyping_mode)] ## default value DISCOVERY
        lines += [' -out_mode %s \\' %(self.output_mode)] ## default value EMIT_VARIANTS_ONLY
        if self.fp_alleles:
            if os.path.isdir(self.fp_alleles):
                s = '%s/$CHROMOSOME.$LSB_JOBINDEX.vcf' %(self.fp_alleles)
                lines += [' --alleles %s \\' %(s)]
            elif os.path.isfile(self.fp_alleles):
                lines += [' --alleles %s \\' %(self.fp_alleles)]
            else:
                print(self.fp_alleles, 'not found')
                sys.exit()

        ## http://www.broadinstitute.org/gatk/gatkdocs/#VariantAnnotatorannotations
        s_annotation = ''

        ## http://gatkforums.broadinstitute.org/discussion/2318/undocumented-change-in-2-4-a-depthofcoverage
        if '2.4' in self.fp_GATK or '2.5' in self.fp_GATK or '2.7' in self.fp_GATK or '2.8' in self.fp_GATK:
            ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_annotator_Coverage.html
            s_annotation += ' --annotation Coverage'
        elif '2.3' in self.fp_GATK:
            ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_annotator_DepthOfCoverage.html
            s_annotation += ' --annotation DepthOfCoverage'
        else:
            print('Unknown version of GATK. Please see why the version number must be known here:')
            print('http://gatkforums.broadinstitute.org/discussion/2318/undocumented-change-in-2-4-a-depthofcoverage')
            sys.exit()

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
        lines += [s_annotation]

        ##
        ## append input files
        ##
        for fp_bam in self.fp_bams:
            for s in os.listdir(fp_bam):
                if s[-4:] != '.bam':
                    continue
                ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_CommandLineGATK.html#--input_file
                lines += [' --input_file %s/%s \\' %(
                    fp_bam,s,)]

        return lines


    def term_cmd(self,analysis_type,l_fp_out,):

        if type(l_fp_out) != list:
            print(l_fp_out)
            stop

        ## cont cmd
##        lines = [';']
        lines = ['\n']
        for fp_out in l_fp_out:
            fp_out = fp_out
##            lines += ['echo %s >> %s.touch;' %(fp_out,analysis_type,)]
##            lines += ['touch touch/%s;' %(fp_out,)]
            lines += ['echo %s >> %s.touch' %(fp_out,analysis_type,)]
            lines += ['touch touch/%s' %(fp_out,)]

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
            s += ' --%s %s' %(k,str(v).replace('$','\\$'))
        fd = open('rerun_python.sh','w')
        fd.write(s)
        fd.close()
        self.execmd('chmod +x rerun_python.sh')

        ## cont cmd
        lines += ['bash ./rerun.sh']
##        ## term cmd
##        lines += ['"']
##
##        lines += ['echo $cmd']
##        lines += ['eval $cmd']

        return lines


    def init_UG_or_HC(self,l_chroms,d_chrom_lens,analysis_type):

        lines = []

        lines += ['\n## parse chromosome from command line']
        lines += ['CHROMOSOME=$1']
        lines += ['LSB_JOBINDEX=$2']

        lines += self.bash_chrom2len(l_chroms,d_chrom_lens,)

        ##
        ## do not allow interval to exceed the length of the chromosome
        ## otherwise it will raise an I/O error (v. 1.4-15)
        ##
        lines += ['posmax=$(((${LSB_JOBINDEX}+0)*%i))' %(self.i_UG_or_HC_size)]
        lines += ['if [ $posmax -gt $LENCHROMOSOME ]']
        lines += ['then posmax=$LENCHROMOSOME']
        lines += ['fi\n']

        fn = 'out_%s/' %(analysis_type)
        fn += '${CHROMOSOME}.${LSB_JOBINDEX}.vcf'

        ## job finished
        lines += ['if [ -s %s.idx ]; then exit; fi\n' %(fn)]

##        ## job started
##        lines += ['if [ -s %s ]; then' %(fn)]
##        lines += ['exit']
##        lines += ['fi\n']

        return lines


    def bash_chrom2len(self,l_chroms,d_chrom_lens,):

        lines = []
        lines += ['\n## define arrays']
        lines += ['CHROMOSOMES=(%s)' %(' '.join(l_chroms))]
        lines += ['LENCHROMOSOMES=(%s)' %(' '.join(
            [str(d_chrom_lens[chrom]) for chrom in l_chroms],
            ),)]

        lines += ['\n## find chromosome index in array of chromosomes']
        lines += ['for ((index=0; index<${#CHROMOSOMES[@]}; index++))\ndo']
        lines += ['if [ "${CHROMOSOMES[$index]}" = "$CHROMOSOME" ]\nthen']
        lines += ['break\nfi\ndone\n']

        lines += ['\n## find length of current chromosome passed via the command line']
        lines += ['LENCHROMOSOME=${LENCHROMOSOMES[$index]}\n']

        return lines


    def parse_arguments(self,):

        parser = argparse.ArgumentParser()

        ##
        ## add arguments
        ##

        parser.add_argument(
            '--bam','--bams','--bamdir','--fp_bams',
            dest='fp_bams',
            help='Path to directory containing improved BAMs',
            nargs='+',
            required = True,
            )

        parser.add_argument(
            '--GATK','--fp_GATK',
            dest='fp_GATK',
            help='File path to GATK (e.g. /software/varinf/releases/GATK/GenomeAnalysisTK-1.4-15-gcd43f01/GenomeAnalysisTK.jar)',
            required = True,
            )

        ## wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.3/b37/human_g1k_v37.fasta.gz
        parser.add_argument(
            '--FASTA','--reference','--reference-sequence','--reference_sequence','--fp_FASTA_reference_sequence',
            dest='fp_FASTA_reference_sequence',
            help='File path to reference sequence in FASTA format (e.g. /lustre/scratch111/resources/vrpipe/ref/Homo_sapiens/1000Genomes/human_g1k_v37.fasta)',
            required = True,
            )

        parser.add_argument(
            '--project','-P','-G',
            dest='project',
            help='Project',
            required = True,
            )

        parser.add_argument(
            '--arguments','--args','--options','--opts','--fp_arguments',
            dest='fp_arguments',
            required = False,
            )

        parser.add_argument(
            '--name','--dataset',
            dest='name',
            help='Name of dataset',
            required = False,
            )

        ##
        ## UnifiedGenotyper / HaplotypeCaller
        ##

        parser.add_argument(
            '--i_UG_or_HC_size',
            dest='i_UG_or_HC_size',
            help='Size (bp) of divided parts.',
            type=int,default=10*10**6,
            required = False,
            )

        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_genotyper_UnifiedGenotyper.html#--genotyping_mode
        parser.add_argument(
            '--genotyping_mode','-gt_mode',
            help='Specifies how to determine the alternate alleles to use for genotyping.',
            default='DISCOVERY', required = False,
            )

        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_genotyper_UnifiedGenotyper.html#--output_mode
        parser.add_argument(
            '--output_mode','-out_mode',
            help='Specifies which type of calls we should output.',
            default='EMIT_VARIANTS_ONLY', required = False,
            )

        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_genotyper_UnifiedGenotyper.html#--alleles
        parser.add_argument(
            '--fp_alleles','--alleles','-alleles',
            help='The set of alleles at which to genotype when --genotyping_mode is GENOTYPE_GIVEN_ALLELES.',
            default='',
            required = False,
            )

        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_CommandLineGATK.html#--intervals
        parser.add_argument(
            '--fp_intervals','--intervals','-L',
            help='Additionally, one may specify a rod file to traverse over the positions for which there is a record in the file (e.g. -L file.vcf).',
            required = False,
            )

        ##
        ## VariantRecalibrator resources
        ##

        parser.add_argument(
            '--resources','--VariantRecalibrator','--fp_resources_SNP','--fp_resources',
            dest='fp_resources_SNP',
            help='File path to a file with -resource lines to append to GATK VariantRecalibrator',
            required = False,
            )

        parser.add_argument(
            '--resources_INDEL','--VariantRecalibrator_INDEL','--fp_resources_INDEL',
            dest='fp_resources_INDEL',
            help='File path to a file with -resource lines to append to GATK VariantRecalibrator',
            required = False,
            )

        ## wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.3/b37/hapmap_3.3.b37.vcf.gz
##        parser.add_argument(
##            '--hapmap',
##            dest='fp_resource_hapmap',
##            help='File path to hapmap vcf to be used by VariantCalibrator (e.g. /lustre/scratch107/projects/uganda/users/tc9/in_GATK/hapmap_3.3.b37.sites.vcf)',
##            type=str,default=None,
##            required = False,
##            )
##

        ## wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.3/b37/1000G_omni2.5.b37.vcf.gz
##        parser.add_argument(
##            '--omni',
##            dest='fp_resource_omni',
##            help='File path to omni vcf to be used by VariantCalibrator (e.g. /lustre/scratch107/projects/uganda/users/tc9/in_GATK/1000G_omni2.5.b37.sites.vcf)',
##            type=str,default=None,
##            required = False,
##            )

        ## dbSNP file. rsIDs from this file are used to populate the ID column of the output.
        ## wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.3/b37/dbsnp_137.b37.vcf.gz
        s_help = 'File path to dbsnp vcf to be used by VariantCalibrator'
        s_help += ' (e.g. /lustre/scratch107/projects/uganda/users/tc9/in_GATK/dbsnp_135.b37.vcf)'
##        s_help += '\nYou can get the vcf file with this command:'
##        s_help += '\ncurl -u gsapubftp-anonymous: ftp.broadinstitute.org/bundle/1.5/b37/dbsnp_135.b37.vcf.gz -o dbsnp_135.b37.vcf.gz; gunzip dbsnp_135.b37.vcf.gz'
        parser.add_argument(
            '--fp_vcf_dbsnp','--dbsnp',
            help=s_help, required = True,
            )

        ##
        ## ApplyRecalibration
        ##

        parser.add_argument(
            '--f_ts_filter_level','--ts','--f_ApplyRecalibration_ts_filter_level','--ts_filter_level',
            dest='f_ts_filter_level',
            type=float,default=None,required=True,
            )

        ##
        ## BEAGLE
        ##
        parser.add_argument(
            '--beagle','--BEAGLE','--BEAGLEjar','--fp_software_beagle',
            dest='fp_software_beagle',
            help='File path to BEAGLE .jar file (e.g. /nfs/team149/Software/usr/share/beagle_3.3.2.jar)',
            required=True,
            )

        parser.add_argument(
            '--i_BEAGLE_size',
            dest='i_BEAGLE_size',
            help='Size (Mbp) of divided parts.',
            type=int,default=2, ## CPU bound (2Mbp=22hrs,3090MB) with lowmem option...
            required = False,
            )

        parser.add_argument(
            '--i_BEAGLE_edge',
            dest='i_BEAGLE_edge',
            help='Window size (kbp) at either side of the divided part to avoid edge effects.',
            type=int,default=150,
            required = False,
            )

        s_help = 'Phased file to be divided (e.g.'
        s_help += ' ALL.chr$CHROM.phase1_release_v3.20101123.filt.renamed.bgl'
        parser.add_argument(
            '--fp_BEAGLE_phased',
            dest='fp_BEAGLE_phased',
            help=s_help,
            required = True,
            )

        s_help = 'Markers file to be divided (e.g.'
        s_help += ' ALL.chr$CHROM.phase1_release_v3.20101123.filt.renamed.markers'
        parser.add_argument(
            '--fp_BEAGLE_markers',
            dest='fp_BEAGLE_markers',
            help=s_help, required=True,
            )

        parser.add_argument(
            '--i_BEAGLE_nsamples',
            dest='i_BEAGLE_nsamples',
            help='',
            default=20,
            required = False,
            )

        ##
        ## optional arguments
        ##
        parser.add_argument(
            '--checkpoint', dest='bool_checkpoint', action='store_true',
            default=False)

        ##
        ##
        ##
        
        ## http://docs.python.org/2/library/argparse.html#argparse.ArgumentParser.parse_args
        ## parse arguments to argparse NameSpace
        self.namespace_args = namespace_args = parser.parse_args()

        ## http://docs.python.org/2/library/functions.html#vars
        for k,v in vars(namespace_args).items():
            setattr(self,k,v)

        if self.fp_GATK is None and self.fp_options is None:
            parser.error('--GATK or --arguments')

        bool_not_found = False
        s_arguments = ''
        for k,v in vars(namespace_args).items():
            s_arguments += '%s %s\n' %(k,v)
            if k[:3] != 'fp_':
                continue
            fp = v
            ## argument not specified
            if fp == None or fp == 'None' or fp == '': continue
            if fp == self.fp_bams:
                l_fp = self.fp_bams
            else:
                if '$CHROMOSOME' in fp:
                    l_fp = [fp.replace('$CHROMOSOME',str(chrom)) for chrom in range(1,22+1)]
                else:
                    l_fp = [fp]
            for fp in l_fp:
##                if not f(fp):
                if not(any([os.path.isfile(fp),os.path.isdir(fp)])):
                    print('file path does not exist:', fp)
                    bool_not_found = True
        if bool_not_found == True:
            sys.exit(0)

        if self.fp_arguments == None or self.fp_arguments == 'None':
            self.fp_arguments = '%s.arguments' %(self.project)
            fd = open(self.fp_arguments,'w')
            fd.write(s_arguments)
            fd.close()
        else:
            fd = open(self.fp_arguments,'r')
            lines = fd.readlines()
            fd.close()
            for line in lines:
                l = line.strip().split()
                k = l[0]
                v = l[1]
                setattr(self,k,v)

        if self.f_ts_filter_level not in [None,'None',]:
            if self.f_ts_filter_level < 1:
                self.f_ts_filter_level *= 100.

        if self.name == None:
            self.name = os.path.basename(self.fp_bams[0])

        return


    def __init__(self,):

        ##
        ## parse command line arguments
        ##
        self.parse_arguments()

        self.d_out = {
            'UnifiedGenotyper':'out_UnifiedGenotyper/$CHROMOSOME.${LSB_JOBINDEX}.vcf',
            'HaplotypeCaller':'out_HaplotypeCaller/$CHROMOSOME.${LSB_JOBINDEX}.vcf',
            'BEAGLE':'out_BEAGLE/$CHROMOSOME/$CHROMOSOME.${LSB_JOBINDEX}.like',
            }

        self.verbose = True

        return


if __name__ == '__main__':
    self = main()
    self.main()

##
## todo20120809: tc9/dg11: add a ReduceReads step before UnifiedGenotyper for speed purposes
## and for better variant calling? cf. slide 14 of https://www.dropbox.com/sh/e31kvbg5v63s51t/ajQmlTL6YH/ReduceReads.pdf
## "VQSR Filters are highly empowered by calling all samples together"
## "Reduced BAMs provides better results for large scale analysis projects ( > 100 samples) because it doesn't require batching."
## http://www.broadinstitute.org/gatk/guide/topic?name=best-practices
## "Even for single samples ReduceReads cuts the memory requirements, IO burden, and CPU costs of downstream tools significantly (10x or more) and so we recommend you preprocess analysis-ready BAM files with ReducedReads."
##


## TODO

##
## todo20120809: tc9/dg11: add a ReduceReads step before UnifiedGenotyper for speed purposes
## and for better variant calling? cf. slide 14 of https://www.dropbox.com/sh/e31kvbg5v63s51t/ajQmlTL6YH/ReduceReads.pdf
## "VQSR Filters are highly empowered by calling all samples together"
## "Reduced BAMs provides better results for large scale analysis projects ( > 100 samples) because it doesn't require batching."
##

## todo20130207: tc9: get rid of orphant functions

## todo20130210: tc9: download markers from http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes.phase1_release_v3/ instead of asking user for their location by default... Make those arguments non-required... --- alternatively download as specified here http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes.phase1_release_v3/READ_ME_beagle_phase1_v3

## todo20130514: tc9: split up UnifiedGenotyper output in mode SNP and INDEL instead of BOTH


## CPU

## todo20130320: tc9: make the function BEAGLE_divide run in "parallel"; i.e. run a process for each chromosome
## todo20130424: tc9: run BEAGLE_divide in parallel for each chromosome


## MEMORY

## todo20130204: tc9: make memory sample size dependent... only tested on 3 datasets with 100 samples each...

## todo20130320: tc9: test memory requirements when more than 100 samples
