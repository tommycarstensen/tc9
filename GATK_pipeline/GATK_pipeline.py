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

tmp_pgs_prob = True ## tmp!!!
version=3
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

        self.VariantRecalibrator(l_chroms,d_chrom_lens,)

        if version == 2:
            self.ApplyRecalibration(l_chroms,d_chrom_lens,)
            self.ProduceBeagleInput()

        ## phased imputation reference panels not available for chrY
        l_chroms.remove('Y')
        ## Exception in thread "main" java.lang.ArrayIndexOutOfBoundsException: -2
        l_chroms.remove('X')
        self.BEAGLE(l_chroms,)

##        ## todo: move BEAGLE unite from IMPUTE2 to here...
##        self.tmp_fill_holes = True
##        self.IMPUTE2(l_chroms,d_chrom_lens,)
####        self.IMPUTE2_unite(l_chroms,d_chrom_lens,)
##
##        self.tmp_fill_holes = False
##        self.IMPUTE2_without_BEAGLE(l_chroms,d_chrom_lens,) ## tmp
####        self.IMPUTE2_without_BEAGLE_unite(l_chroms,d_chrom_lens,) ## tmp

        return


    def IMPUTE2_without_BEAGLE_unite(self,l_chroms,d_chrom_lens,):

        ## tmp function!!!

        for chrom in l_chroms:
            index_max = int(math.ceil(
                    d_chrom_lens[chrom]/float(self.i_IMPUTE2_size)))
            cmd = 'cat'
            for index in range(1,index_max+1):
                fp_gen = 'out_IMPUTE2_without_BEAGLE/%s/%s.%i.gen' %(chrom,chrom,index)
                if not os.path.isfile(fp_gen): continue
                cmd += ' %s' %(fp_gen)
            fp_out = 'out_IMPUTE2_without_BEAGLE/%s.gen' %(chrom)
            if os.path.isfile(fp_out): continue
            cmd += ' > %s' %(fp_out)
            self.execmd(cmd)

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
            l_fn = os.listdir('stderr/%s' %(dn[4:],))
            for fn in l_fn:
                fp = 'stderr/%s/%s' %(dn[4:],fn)
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


    def IMPUTE2_fileIO_checks(self,chrom):

##        cmd = """awk 'FNR==1{print FILENAME":"$1}' out_IMPUTE2/*/*.gen """
##        cmd += """| awk -F ":" '{if($3%5000000>100000) print $3%5000000,$0}'"""
##        s = os.popen(cmd).read().strip()
##        if s != '':
##            print cmd
##            print 'WARNING'
##            print s

        ## in
        cmd = 'cat in_IMPUTE2/%s.gen' %(chrom)
        cmd += " | awk '{"
        cmd += ' max=0;'
        cmd += ' for (i=6; i<=NF; i++) {if($i>max) {max=$i}};'
        if tmp_pgs_prob == True:
            cmd += 'print $1'
        else:
            cmd += " if(max>=0.9) {print $1}"
        cmd += " }'"
        cmd += ' > in%s' %(chrom)
        self.execmd(cmd)

        fp_legend = self.fp_impute2_legend
        fp_legend = fp_legend.replace('$CHROMOSOME','${CHROMOSOME}')
        fp_legend = fp_legend.replace('${CHROMOSOME}',chrom)
        cmd = 'zcat %s' %(fp_legend)
        cmd += """ | awk -v CHROM=%s 'NR>1{print CHROM":"$2}'""" %(chrom)
        cmd += ' >> in%s' %(chrom)
        self.execmd(cmd)

        cmd = 'sort -u in%s -o in%s' %(chrom,chrom)
        self.execmd(cmd)

        ## out
        cmd = 'cat out_IMPUTE2/%s/%s.*.gen' %(chrom,chrom)
        cmd += """ | awk -v CHROM=%s '""" %(chrom)
        cmd += '{'
##        cmd += ' max=0;'
##        cmd += ' for (i=6; i<=NF; i++) {if($i>max) {max=$i}};'
##        cmd += ' if(max>0.9) {print CHROM":"$3}'
        cmd += 'print CHROM":"$3'
        cmd += " }'"
        cmd += ' | sort -u > out%s' %(chrom)
        self.execmd(cmd)

        cmd = 'comm -3 in%s out%s' %(chrom,chrom)
        i = int(os.popen('%s | wc -l' %(cmd)).read())
        if i > 0:
            print(os.popen('%s | head' %(cmd)).readlines())
            print(i)
            print(cmd)
            sys.exit()

        os.remove('in%s' %(chrom))
        os.remove('out%s' %(chrom))

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


    def IMPUTE2(self,l_chroms,d_chrom_lens,):

        ## http://mathgen.stats.ox.ac.uk/impute/impute_v2.html
        ## http://mathgen.stats.ox.ac.uk/impute/impute_v2_instructions.pdf
        ## http://mathgen.stats.ox.ac.uk/impute/example_one_phased_panel.html
        ## http://www.stats.ox.ac.uk/~marchini/software/gwas/file_format.html

        ##
        ## 1) check input existence
        ##
        l_fp_in = self.IMPUTE2_parse_input_files()
        bool_exit = self.check_in('BEAGLE',l_fp_in,)

        ##
        ## 2) touch
        ##
        bool_return = self.touch('IMPUTE2')
        if bool_return == True: return

        ##
        ## BEAGLE gprobs > IMPUTE2 gen
        ##
        for chrom in l_chroms:
            print('gprobs2gen', chrom)
            fp_in = 'out_BEAGLE/%s.gprobs' %(chrom)
            fp_out = 'in_IMPUTE2/%s.gen' %(chrom)
            if os.path.isfile(fp_out): continue
            self.BEAGLE_fileIO_checks(chrom)
            self.BEAGLE_unite(chrom)
            self.gprobs2gen(fp_in,fp_out)

##        print('exiting IMPUTE2 temporarily')
##        sys.exit() ## tmp!!!

        fp_in = 'in_IMPUTE2/$CHROMOSOME.gen'
        fp_out = self.d_out['IMPUTE2']

        ## write shell script to list/memory
        lines = self.IMPUTE2_write_shell_script(fp_in,fp_out)

        ## term cmd
        lines += self.term_cmd('IMPUTE2',[fp_out],)

        ## write shell script to file
        self.write_shell('shell/IMPUTE2.sh',lines,)

        ##
        ## execute shell script
        ##
        for chrom in l_chroms:
            print('bsub IMPUTE2', chrom)
            self.bsub_IMPUTE2('IMPUTE2',chrom,d_chrom_lens,)

        return


    def bsub_IMPUTE2(self,dn_suffix,chrom,d_chrom_lens,):

        ## fgrep CPU */stdout/IMPUTE2/*.out | awk '{print $5}' | sort -nr | head -n1
        ## 41013.13
        queue = 'long'

        s_legend = self.fp_impute2_legend
        s_legend = s_legend.replace('$CHROMOSOME','${CHROMOSOME}')
        s_legend = s_legend.replace('${CHROMOSOME}',chrom)
        cmd = 'cat in_IMPUTE2/%s.gen' %(chrom)
        cmd += " | awk '{cnt[int($3/%i)]++}" %(self.i_IMPUTE2_size)
        cmd += " END{for(j in cnt) print j+1,cnt[j]}'"
        cmd += ' | sort -k1n,1'
        d_index2variants = dict([s.split() for s in os.popen(
            cmd).read().strip().split('\n')])
        for index,variants in d_index2variants.items():
            ## memory requirement/usage different if not 100 samples?
            dev = fudgefactor = 800
            intersect = 255
##            slope = 4247.268 ## -pgs ???
##            slope = 8200 ## -pgs_prob ???
            slope = 8600 ## -pgs_prob and -fill_holes ???
            memMB = dev+intersect+int(slope*int(variants)/100000.)
            fp_in = 'in_%s/%s.gen' %(dn_suffix,chrom)
            fp_out = 'out_%s/%s/%s.%s.gen' %(dn_suffix,chrom,chrom,index)
            if os.path.getsize(fp_in) == 0: continue ## redundant
            if os.path.isfile(fp_out):
                if os.path.getsize(fp_out) > 0:
                    continue ## redundant
##            index_max = int(math.ceil(
##                    d_chrom_lens[chrom]/float(self.i_IMPUTE2_size)))
##            if os.path.isfile('%s_summary' %(fp_out)):
##                if not 'ERROR' in open('%s_summary' %(fp_out)).read():
##                    continue
            if os.path.isfile('stdout/%s/%s.%s.%s.out' %(dn_suffix,dn_suffix,chrom,index)):
                continue
            J = '%s.%s.%s' %(dn_suffix,chrom,index,)
            std_suffix = '%s/%s.%s.%s' %(dn_suffix,dn_suffix,chrom,index,)
            cmd = self.bsub_cmd(
                dn_suffix,J,std_suffix=std_suffix,
                chrom=chrom, chromlen=d_chrom_lens[chrom], index = index,
                memMB=memMB,queue=queue,)
            os.system(cmd)

        return


    def gprobs2gen(self,fp_in,fp_out,):
    
        cmd = 'cat %s' %(fp_in)
        cmd += " | awk '"
        cmd += 'BEGIN{FS=":"}'
        cmd += 'NR>1{print $1,$2}'
        cmd += "'"
        cmd += " | awk '"
        cmd += '{print $1":"$2,$1":"$2,$0}'
        cmd += "'"
        cmd += ' | cut -d " " -f 3 --complement'
        cmd += ' > %s' %(fp_out,)
        print(cmd)
        self.execmd(cmd)

        return


    def IMPUTE2_parse_input_files(self,):

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
            index = int(l[1])
            fp_in = 'out_BEAGLE/%s/' %(chrom,)
            fp_in += '%s.%i.like.gprobs' %(chrom,index,)
            l_fp_in += [fp_in]

        return l_fp_in

    
    def IMPUTE2_write_shell_script(self,fp_in,fp_out,):

        ## initiate shell script
        lines = ['#!/bin/bash\n']

        ## parse chromosome and chromosome length from command line
        lines += ['CHROMOSOME=$1']
        lines += ['CHROM=$1']
##        lines += self.bash_chrom2len(l_chroms,d_chrom_lens,)
        lines += ['LENCHROMOSOME=$2']
        lines += ['LSB_JOBINDEX=$3']

        lines += ['if [ -s %s ]; then\nexit\nfi\n' %(fp_out)] ## redundant

        ##
        ## define files
        ##
        ## Input file options
        ## http://mathgen.stats.ox.ac.uk/impute/ALL_1000G_phase1integrated_v3_impute.tgz
        s_haps = self.fp_impute2_hap.replace('$CHROMOSOME','${CHROMOSOME}')
        s_legend = self.fp_impute2_legend.replace('$CHROMOSOME','${CHROMOSOME}')

        lines += ['if [ $CHROMOSOME = "X" ]; then']
        lines += ['map=%s' %(self.fp_impute2_map.replace(
            '$CHROMOSOME','\\${CHROMOSOME}_PAR1'))]
        lines += ['else']
        lines += ['map=%s' %(self.fp_impute2_map.replace(
            '$CHROMOSOME','\\${CHROMOSOME}'))]
        lines += ['fi\n']

        ## http://mathgen.stats.ox.ac.uk/impute/impute_v2.html
        ## IMPUTE2 requires that you specify an analysis interval
        ## in order to prevent accidental whole-chromosome analyses.
        ## If you want to impute a region larger than 7 Mb
        ## (which is not generally recommended),
        ## you must activate the -allow_large_regions flag.
        lines += ['posmax=$(((${LSB_JOBINDEX}+0)*%i))' %(self.i_IMPUTE2_size)]
        lines += ['if [ $posmax -gt $LENCHROMOSOME ]; then']
        lines += ['posmax=$LENCHROMOSOME']
        lines += ['fi\n']

        ## init cmd
        lines += ['cmd="']

        lines += self.body_IMPUTE2(
            fp_in,s_legend,s_haps,fp_out,)

        return lines


    def body_IMPUTE2(self,s_gen,s_legend,s_haps,fp_out,):

        ##
        ## initiate impute2 command
        ##
        lines = ['%s \\' %(self.fp_software_IMPUTE2)]

        ##
        ## Required arguments
        ## http://mathgen.stats.ox.ac.uk/impute/required_arguments.html
        ##
        ## http://mathgen.stats.ox.ac.uk/impute/input_file_options.html#-g
        lines += ['-g %s \\' %(s_gen)]
        ## http://mathgen.stats.ox.ac.uk/impute/required_arguments.html#-m
        ## http://mathgen.stats.ox.ac.uk/impute/input_file_options.html#-m
        ## Ask Deepti whether a different recombination map file should be used
        lines += ['-m $map \\']

        ##
        ## Basic options
        ## http://mathgen.stats.ox.ac.uk/impute/impute_v2.html#basic_options
        ##
        ## http://mathgen.stats.ox.ac.uk/impute/basic_options.html#-int
        lines += ['-int $(((${LSB_JOBINDEX}-1)*%i+1)) $posmax \\' %(self.i_IMPUTE2_size)]
        lines += ['-buffer 250',]
        lines += ['-Ne 20000',]
        lines += ['-call_thresh 0.9',]

        ##
        ## Output file options
        ## http://mathgen.stats.ox.ac.uk/impute/output_file_options.html
        ##
        lines += ['-o %s \\' %(fp_out)]
        ## http://mathgen.stats.ox.ac.uk/impute/output_file_options.html#-pgs
        lines += ['-pgs \\']
        if tmp_pgs_prob == True:
            lines = lines[:-1] ## ERROR: The -pgs_prob flag is not compatible with the -pgs flag.
            lines += ['-pgs_prob \\']
            lines += ['-prob_g %s \\' %(s_gen)]
            if self.tmp_fill_holes == True:
                lines += ['-fill_holes \\']

        ##
        ## Input file options
        ##
        ## http://mathgen.stats.ox.ac.uk/impute/input_file_options.html#-l
        lines += ['-l %s \\' %(s_legend,)]
        ## http://mathgen.stats.ox.ac.uk/impute/input_file_options.html#-h
        lines += ['-h %s \\' %(s_haps,)]

        ##
        ## MCMC OPTIONS
        ##
        ## http://mathgen.stats.ox.ac.uk/impute/impute_v2.html#mcmc_options
        lines += ['-iter 30 \\']
        lines += ['-burnin 10 \\']
        lines += ['-k 80 \\']
        lines += ['-k_hap 500 \\']

        return lines


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


    def ProduceBeagleInput(self,):

        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_beagle_ProduceBeagleInput.html

        ## "After variants are called and possibly filtered,
        ## the GATK walker ProduceBeagleInputWalker will take the resulting VCF as input,
        ## and will produce a likelihood file in BEAGLE format."

        T = analysis_type = 'ProduceBeagleInput'
        ##fula_20120704/stdout/ProduceBeagleInput/ProduceBeagleInput.out:    Max Memory :      2154 MB
        ##uganda_20130113/stdout/ProduceBeagleInput/ProduceBeagleInput.out:    Max Memory :      2168 MB
        ##zulu_20121208/stdout/ProduceBeagleInput/ProduceBeagleInput.out:    Max Memory :      2171 MB
        memMB = 1900
        ##fula_20120704/stdout/ProduceBeagleInput/ProduceBeagleInput.out:    CPU time   :   2483.24 sec.
        ##zulu_20121208/stdout/ProduceBeagleInput/ProduceBeagleInput.out:    CPU time   :   4521.89 sec.
        ##uganda_20130113/stdout/ProduceBeagleInput/ProduceBeagleInput.out:    CPU time   :   4948.94 sec.
        queue = 'normal'

        ##
        ## 1) check file in existence
        ##
        fp_in = self.d_in['ProduceBeagleInput']
        l_fp_in = []
        for mode in ['SNP',]:
            fp_in = 'out_ApplyRecalibration/ApplyRecalibration.recalibrated.filtered.%s.vcf' %(mode)
            l_fp_in += [fp_in]
        bool_exit = self.check_in('ApplyRecalibration',l_fp_in,)

##        l_fp_in = ['out_SelectVariants/SelectVariants.vcf'] ## tmp!!!

        ##
        ## 2) touch
        ##
        bool_return = self.touch(analysis_type)
        if bool_return == True: return

        fp_out = self.d_out['ProduceBeagleInput']

        ##
        ## 3) init script
        ##
        lines = ['#!/bin/bash']

        lines += self.init_cmd(analysis_type,memMB,)

        ## GATKwalker, required, out
        lines += ['--out %s \\' %(fp_out,)]
        ## GATKwalker, required, in
        lines += ['--variant %s \\' %(l_fp_in[0])]

##        l_fp_out = ['out_%s/%s.%i.bgl' %(T,T,chrom) for chrom in l_chrom]
        lines += self.term_cmd(analysis_type,[fp_out],)

        self.write_shell('shell/ProduceBeagleInput.sh',lines,)

        ## execute shell script
        J = 'PBI'
        cmd = self.bsub_cmd(analysis_type,J,memMB=memMB,)
        self.execmd(cmd)

        return


    def parse_marker(self,line_markers,):

        l_markers = line_markers.split()
        pos_phased = int(l_markers[1])
        alleleA_phased = l_markers[2]
        alleleB_phased = l_markers[3]

        return pos_phased, alleleA_phased, alleleB_phased


    def new_BEAGLE_divide(self,l_chroms,):

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
                l_files = glob.glob('%s/UnifiedGenotyper.%s.*.vcf' %(path,chrom,))
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

            ## read equivalent VR line
            while True:
                line_recal = fd_recal.readline()
                if line_recal[0] != '#':
                    break

            ## skip non-diallelic SNPs
            if bool_diallelic == False:
                continue

            ## skip if VQSLOD value below threshold
            l_recal = line_recal.split('\t')
            l_recal_INFO = l_recal[7].split(';')
            VQSLOD = float(l_recal_INFO[1].split('=')[1])
            if VQSLOD < minVQSLOD:
                continue

            ## tmp!!!
            if l_recal[0] != l_vcf[0] or l_recal[1] != l_vcf[1]:
                print('VRrecal')
                print(l_recal)
                print('UGvcf')
                print(line_vcf)
                stoptmp

            yield l_vcf


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

        header = self.vcf2beagle_header(l_files[0])
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
##        pos_phased = 0
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
        line_phased = fd_phased.readline()
        line_markers = fd_markers.readline()
        pos_phased, alleleA_phased, alleleB_phased = self.parse_marker(
            line_markers)

        with fileinput.input(files=l_files) as fd_vcf:
            while True:
                try:
                    l_vcf = next(self.generate_line_vcf_PASS_split(
                        fd_vcf,fd_recal,minVQSLOD,))
                except StopIteration:
                    break

                line_beagle = self.vcf2beagle(l_vcf)

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
                    if pos_phased == position:
                        ##
                        ## 1a) markers identical
                        ##
                        if (
                            alleleA_phased == alleleA
                            and
                            alleleB_phased == alleleB
                            ):
                            bool_append_markphas = True
                            pass
                        ##
                        ## 1b) markers different
                        ##
                        else:
                            bool_append_markphas = False
                            pos_phased_prev = pos_phased
                            while True:
                                ## read markers and phased
                                line_phased = fd_phased.readline()
                                line_markers = fd_markers.readline()
                                ## last marker is in the genotype probability file and not the markers file
                                ## e.g. 7:159128574 in fula4x
                                if line_markers == '':
                                    bool_append_markphas = False
                                    break
                                (pos_phased, alleleA_phased, alleleB_phased
                                 ) = self.parse_marker(line_markers)
                                if pos_phased > pos_phased_prev:
                                    break
                                else:
                                    continue
                                continue
                            pass
                        break
                    ##
                    ## 2) continue loop over genotype likelihoods ("panel 2")
                    ##
                    elif position < pos_phased:
                        bool_append_markphas = False
                        break
                    ##
                    ## 3) continue loop over markers ("panel 0")
                    ##
                    ## elif position > pos_phased:
                    else:
                        ## INDEL
                        if pos_phased == pos_prev:
                            lines_out_markers1 += [line_markers]
                            lines_out_phased1 += [line_phased]
                            if bool_BOF2 == True and pos_phased > pos_term1-edge:
                                lines_out_markers2 += [line_markers]
        ##                            lines_out_phased2 += [line_phased]
                        ## SNP unique to panel 0
        ##                    elif pos_prev < pos_phased:
                        else:
                            lines_out_markers1 += [line_markers]
                            lines_out_phased1 += [line_phased]
                            if bool_BOF2 == True:
                                lines_out_markers2 += [line_markers]
        ##                            lines_out_phased2 += [line_phased]
                        ## append marker not present in genotype likehoods file
                        ## read markers and phased
                        line_phased = fd_phased.readline()
                        line_markers = fd_markers.readline()
                        if line_markers == '':
                            bool_append_markphas = False
                            break
                        pos_phased, alleleA_phased, alleleB_phased = self.parse_marker(
                            line_markers)
                        continue
                    continue


                ##
                ##
                ##
                ## BOF2
                if bool_BOF2 == True:
                    lines_out2 += [line_beagle]
                    ## match
                    if bool_append_markphas == True:
                        lines_out_markers2 += [line_markers]
    ##                    lines_out_phased2 += [line_phased]
                    ## mismatch
                    else: ## ms23/dg11 2013mar12
                        lines_out_markers2 += ['%s:%s %s %s %s\n' %(
                            chrom,position,position,alleleA,alleleB,)]
                if bool_EOF1 == False:
                    lines_out1 += [line_beagle]
                    ## match
                    if bool_append_markphas == True:
                        lines_out_markers1 += [line_markers]
                        lines_out_phased1 += [line_phased]
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
                            lines_out_phased1 += [line_phased]
                        lines_out_markers1 = lines_out_markers2
                    ## append to current lines
                    ## elif bool_append_to_next == True
                    else:
                        lines_out1 += [line]
                        if bool_append_markphas == True:
                            lines_out_markers1 += [line_markers]
                            lines_out_phased1 += [line_phased]
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
                    line_markers = fd_markers.readline()
                    line_phased = fd_phased.readline()
                    pos_phased, alleleA_phased, alleleB_phased = self.parse_marker(
                        line_markers)

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
##        elif l_vcf[4] in ','.join(tup) for tup in itertools.permutations('ACGT',2):
        elif l_vcf[4] in iter(','.join(tup) for tup in itertools.combinations('ACGT',2)):
            return False,False
##        elif l_vcf[4] in ','.join(tup) for tup in itertools.permutations('ACGT',3):
        elif l_vcf[4] in iter(','.join(tup) for tup in itertools.combinations('ACGT',3)):
            return False,False
        ## skip insertions
        else:
            return True,None

        return


    def BEAGLE_divide(self,chrom,):

        ## todo - instead of using "xxx.renamed.markers", just rename on the fly
        ## out: chrom:pos pos alleleA alleleB\n

        ##
        ## size and edge
        ##
        size = self.i_BEAGLE_size*1000000
        edge = self.i_BEAGLE_edge*1000

        fp_in = 'out_ProduceBeagleInput/ProduceBeagleInput.%s.bgl' %(chrom)
        ## genotype likelihood file
        fp_out_prefix = 'in_BEAGLE/%s/%s' %(chrom,chrom,)

        fp_phased = self.fp_BEAGLE_phased.replace('$CHROMOSOME',chrom)
        fp_markers = self.fp_BEAGLE_markers.replace('$CHROMOSOME',chrom)

        ##
        ## parse genotype likelihood file header
        ##
        ## append header
        fd = open('out_ProduceBeagleInput/ProduceBeagleInput.bgl','r')
        header = fd.readline()
        fd.close()

##        ## parse header
##        header_phased = fd_phased.readline()

        ##
        ## open files
        ##
        fd_phased = open(fp_phased,'r')
        fd_markers = open(fp_markers,'r')
        fd_bgl = open(fp_in,'r') ## with header marker allelelA alleleB
##        fd_bgl.readline() ## skip header

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

        d_index2pos = {}
        index = 1 ## LSF does not allow 0 for LSB_JOBINDEX! "Bad job name. Job not submitted." ## e.g. bsub -J"test[0-2]" echo A
##        pos_phased = 0
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
        line_phased = fd_phased.readline()
        line_markers = fd_markers.readline()
        pos_phased, alleleA_phased, alleleB_phased = self.parse_marker(
            line_markers)

        ## loop over BEAGLE genotype likelihood file lines
        for line in fd_bgl:
            l = line.strip().split()
            ## assume markers to be formatted like CHROM:POS
            l_chrom_pos = l[0].split(':')
            chrom_like = chrom = l_chrom_pos[0]
            pos_panel2 = pos_like = pos = position = int(l_chrom_pos[1])
            alleleA_like = l[1]
            alleleB_like = l[2]
            ## ignore INDELs
##            if (
##                alleleA_like not in ['A','C','G','T',] or
##                alleleB_like not in ['A','C','G','T',]
##                ): continue
            if len(alleleA_like) > 1 or len(alleleB_like) > 1:
                continue

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
                if pos_phased == position:
                    ##
                    ## 1a) markers identical
                    ##
                    if (
                        alleleA_phased == alleleA_like
                        and
                        alleleB_phased == alleleB_like
                        ):
                        bool_append_markphas = True
                        pass
                    ##
                    ## 1b) markers different
                    ##
                    else:
                        bool_append_markphas = False
                        pos_phased_prev = pos_phased
                        while True:
                            ## read markers and phased
                            line_phased = fd_phased.readline()
                            line_markers = fd_markers.readline()
                            ## last marker is in the genotype probability file and not the markers file
                            ## e.g. 7:159128574 in fula4x
                            if line_markers == '':
                                bool_append_markphas = False
                                break
                            (pos_phased, alleleA_phased, alleleB_phased
                             ) = self.parse_marker(line_markers)
                            if pos_phased > pos_phased_prev:
                                break
                            else:
                                continue
                            continue
                        pass
                    break
                ##
                ## 2) continue loop over genotype likelihoods ("panel 2")
                ##
                elif position < pos_phased:
                    bool_append_markphas = False
                    break
                ##
                ## 3) continue loop over markers ("panel 0")
                ##
                ## elif position > pos_phased:
                else:
                    ## INDEL
                    if pos_phased == pos_prev:
                        lines_out_markers1 += [line_markers]
                        lines_out_phased1 += [line_phased]
                        if bool_BOF2 == True and pos_phased > pos_term1-edge:
                            lines_out_markers2 += [line_markers]
##                            lines_out_phased2 += [line_phased]
                    ## SNP unique to panel 0
##                    elif pos_prev < pos_phased:
                    else:
                        lines_out_markers1 += [line_markers]
                        lines_out_phased1 += [line_phased]
                        if bool_BOF2 == True:
                            lines_out_markers2 += [line_markers]
##                            lines_out_phased2 += [line_phased]
                    ## append marker not present in genotype likehoods file
                    ## read markers and phased
                    line_phased = fd_phased.readline()
                    line_markers = fd_markers.readline()
                    if line_markers == '':
                        bool_append_markphas = False
                        break
                    pos_phased, alleleA_phased, alleleB_phased = self.parse_marker(
                        line_markers)
                    continue
                continue

            ## BOF2
            if bool_BOF2 == True:
                lines_out2 += [line]
                ## match
                if bool_append_markphas == True:
                    lines_out_markers2 += [line_markers]
##                    lines_out_phased2 += [line_phased]
                ## mismatch
                else: ## ms23/dg11 2013mar12
                    lines_out_markers2 += ['%s:%s %s %s %s\n' %(
                        chrom,position,position,alleleA_like,alleleB_like,)]
            if bool_EOF1 == False:
                lines_out1 += [line]
                ## match
                if bool_append_markphas == True:
                    lines_out_markers1 += [line_markers]
                    lines_out_phased1 += [line_phased]
                ## mismatch
                else: ## ms23/dg11 2013mar12
                    lines_out_markers1 += ['%s:%s %s %s %s\n' %(
                        chrom,position,position,alleleA_like,alleleB_like,)]
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
                        lines_out_phased1 += [line_phased]
                    lines_out_markers1 = lines_out_markers2
                ## append to current lines
                ## elif bool_append_to_next == True
                else:
                    lines_out1 += [line]
                    if bool_append_markphas == True:
                        lines_out_markers1 += [line_markers]
                        lines_out_phased1 += [line_phased]
                    else: ## ms23/dg11 2013mar12
                        lines_out_markers1 += ['%s:%s %s %s %s\n' %(
                            chrom,position,position,alleleA_like,alleleB_like,)]
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
                line_markers = fd_markers.readline()
                line_phased = fd_phased.readline()
                pos_phased, alleleA_phased, alleleB_phased = self.parse_marker(
                    line_markers)

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
        fd_bgl.close()
        fd_markers.close()
        fd_phased.close()

        ##
        ## file i/o checks
        ##
        self.BEAGLE_divide_fileIO_checks(chrom,fp_phased,)

        return d_index2pos


    def BEAGLE_divide_fileIO_checks(self,chrom,fp_phased,):

        cmd = 'cat %s' %(fp_phased)
        cmd += "| awk 'NR>1{print $2}' | sort -u > panel0in%s" %(chrom)
        self.execmd(cmd)

        cmd = "awk 'FNR>1{print $2}' in_BEAGLE/%s/%s.phased" %(chrom,chrom)
        cmd += "| sort -u > panel0out%s" %(chrom)
        self.execmd(cmd)

        cmd = 'cat out_ProduceBeagleInput/ProduceBeagleInput.%s.bgl' %(chrom)
        cmd += " | awk '{if("
        cmd += '($2=="A"||$2=="C"||$2=="G"||$2=="T")'
        cmd += ' && '
        cmd += '($3=="A"||$3=="C"||$3=="G"||$3=="T")'
        cmd += ") print $1}' | sort > panel2in%s" %(chrom)
        self.execmd(cmd)

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


    def BEAGLE(self,l_chroms,):

        ## http://faculty.washington.edu/browning/beagle/beagle_3.3.2_31Oct11.pdf

        ## http://www.broadinstitute.org/gatk/guide/article?id=43
        ## Interface with BEAGLE imputation software - GSA
        ## IMPORTANT: Due to BEAGLE memory restrictions,
        ## it's strongly recommended that BEAGLE be run on a separate chromosome-by-chromosome basis.
        ## In the current use case, BEAGLE uses RAM in a manner approximately proportional to the number of input markers.

        ## dg11: "Then run beagle imputation
        ## (you may have to chunk up for imputation
        ## and use the "known" option with the 1000G reference)"
        ## ask Deepti what the "known" option is...

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
        if version == 2:
            fp_in = self.d_in['BEAGLE']
            bool_exit = self.check_in('ProduceBeagleInput',[fp_in,],)
        elif version == 3:
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
        ## chunk up for imputation due to memory requirements...
        ##
        ## split bgl by chromosome
        ## this should be part of the BEAGLE_divide function
        ## to avoid looping over the same lines twice (a bit stupid at the moment)
        if version == 2:
            cmd = 'cat out_ProduceBeagleInput/ProduceBeagleInput.bgl '
            cmd += ''' | awk 'BEGIN{FS=":"} '''
##        cmd += " | awk -F : '"
            cmd += 'NR>1{print>"out_ProduceBeagleInput/ProduceBeagleInput."$1".bgl"'
            cmd += "}'"
            print(cmd)
            if not os.path.isfile('out_ProduceBeagleInput/ProduceBeagleInput.22.bgl'): ## redundant
                os.system(cmd)
        ## split further into smaller fragments
        d_indexes = {}
        d_chrom_lens = self.parse_chrom_lens()

        if not (
            os.path.isfile('BEAGLE_divide_indexes.txt')
            and
            os.path.isdir('in_BEAGLE')
            ):
            if version == 3:
                d_indexes = self.new_BEAGLE_divide(l_chroms)
            elif version == 2:
                for chrom in l_chroms:
                    d_index2pos = self.BEAGLE_divide(chrom,)
                    d_indexes[chrom] = d_index2pos
                    continue
                    ## clean up
                    os.remove(
                        'out_ProduceBeagleInput/ProduceBeagleInput.%s.bgl' %(
                            chrom))
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
                std_suffix = '%s/%s.%s.%%I' %('BEAGLE','BEAGLE',chrom)
                cmd = self.bsub_cmd(
                    'BEAGLE',J,memMB=memMB,std_suffix=std_suffix,chrom=chrom,
                    queue=queue,)
                os.system(cmd)

        return


    def BEAGLE_write_shell_script(self,memMB,):

        fp_out = self.d_out['BEAGLE']

        ## initiate shell script
        lines = ['#!/bin/bash\n']

        ## parse chromosome from command line
        lines += ['CHROMOSOME=$1\n']

        lines += ['if [ -s %s.gprobs ]; then\nexit\nfi\n' %(fp_out)] ## redundant

        ## init cmd
        lines += ['cmd="']

        ##
        ## initiate BEAGLE
        ##
        lines += ['java -Djava.io.tmpdir=tmp -Xmx%im -jar %s \\' %(
            memMB,self.fp_software_beagle)]

        lines += self.body_BEAGLE(fp_out,)

        ##
        ## terminate BEAGLE
        ##
        lines += [';']

        ##
        ## gunzip the output files after runs to completion
        ##
        l = []
        for suffix in ['dose','phased','gprobs',]:
            fp = '%s' %(fp_out)
            fp += '.%s.gz' %(
                suffix,
                )
            l += ['gunzip %s' %(fp)]
        s = ';'.join(l).replace('$CHROMOSOME','${CHROMOSOME}')
        lines += [s]

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
            'stdout','stderr',
            'shell',
            'out_UnifiedGenotyper',
            'out_VariantRecalibrator',
            'out_ApplyRecalibration',
            'out_ProduceBeagleInput',
            'in_BEAGLE','in_IMPUTE2',
            'out_BEAGLE','out_IMPUTE2',
            ]

        ## create subdirs
        for dn in l_dn:
            if not os.path.isdir(dn):
                os.mkdir(dn)
            if dn != 'touch' and (dn[:4] == 'out_' or dn[:3] == 'in_'):
                if not os.path.isdir(os.path.join('touch',dn)):
                    os.mkdir(os.path.join('touch',dn))
        for std in ['out','err',]:
            for dn in l_dn:
                if dn[:3] != 'out': continue
                if not os.path.isdir('std%s/%s' %(std,dn[4:],)):
                    os.mkdir('std%s/%s' %(std,dn[4:],))

        for chrom in l_chroms:
            for dn in ['in_BEAGLE','out_BEAGLE','out_IMPUTE2',]:
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
                ## append files in chromosomal subdirectories (e.g. BEAGLE and IMPUTE2)
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
                fp_in = 'out_%s/%s.%s.%i.vcf' %(
                    'UnifiedGenotyper','UnifiedGenotyper',
                    chrom,interval,
                    )
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
        memMB = 27500 ## tmp!!! union4x
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
            l_chroms,d_chrom_lens,self.i_UG_size,)
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
            lines = self.init_cmd(analysis_type,memMB,)

            ## GATKwalker, required, in
            lines += ['--input %s \\' %(vcf) for vcf in l_vcfs_in]
            ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_variantrecalibration_VariantRecalibrator.html#--use_annotation
            if mode == 'SNP':
                lines += [
                    '--use_annotation QD \\',
                    '--use_annotation HaplotypeScore \\',
                    '--use_annotation MQRankSum \\',
                    '--use_annotation ReadPosRankSum \\',
                    '--use_annotation MQ \\',
                    '--use_annotation FS \\',
                    '--use_annotation DP \\',
                    ]
            elif mode == 'INDEL':
                lines += ['-an DP -an FS -an ReadPosRankSum -an MQRankSum \\',]

            ##
            ## GATKwalker, required, out
            ##
            lines += ['--recal_file %s \\' %(fp_recal)]
            lines += ['--tranches_file %s \\' %(fp_tranches)]

            ##
            ## GATKwalker, optional, in
            ##
            
            ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_variantrecalibration_VariantRecalibrator.html#--mode
            lines += ['--mode %s \\' %(mode)]

            ## http://gatkforums.broadinstitute.org/discussion/1259/what-vqsr-training-sets-arguments-should-i-use-for-my-specific-project
            if mode == 'INDEL':
                ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_variantrecalibration_VariantRecalibrator.html#--maxGaussians
                lines += ['--maxGaussians 4 \\'] ## default 10
                ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_variantrecalibration_VariantRecalibrator.html#--percentBadVariants
                lines += ['--percentBadVariants 0.01 \\'] ## default 0.03
                ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_variantrecalibration_VariantRecalibrator.html#--minNumBadVariants
                lines += ['--minNumBadVariants 1000 \\'] ## default 2500

            ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_variantrecalibration_VariantRecalibrator.html#--resource
            fd = open(d_resources[mode],'r')
            lines_resources = fd.readlines()
            fd.close()
            lines += ['%s \\' %(line.strip()) for line in lines_resources]

            ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_variantrecalibration_VariantRecalibrator.html#--TStranche
            l_TStranches = []
    ##        l_TStranches += [99.70+i/20. for i in range(6,0,-1,)]
            l_TStranches += [99+i/10. for i in range(10,0,-1,)]
            l_TStranches += [90+i/2. for i in range(18,-1,-1,)]
            l_TStranches += [70+i for i in range(19,-1,-1,)]
            s_TStranches = ''
            for TStranche in l_TStranches:
                s_TStranches += '--TStranche %.1f ' %(TStranche)
            lines += ['%s \\' %(s_TStranches)]

            ##
            ## GATKwalker, optional, out
            ##
            lines += ['--rscript_file out_%s/%s.%s.plots.R \\' %(T,T,mode,)]

            ##
            ## term GATK walker
            ##
            lines += self.term_cmd(
                '%s.%s' %(analysis_type,mode),[fp_tranches,fp_recal,],)

            self.write_shell('shell/%s.%s.sh' %(analysis_type,mode,),lines,)

            J = 'VR'
            cmd = self.bsub_cmd(
                '%s.%s' %(analysis_type,mode),J,memMB=memMB,queue=queue,
                std_suffix='%s/%s.%s' %(analysis_type,analysis_type,mode,),)
            self.execmd(cmd)

        return


    def ApplyRecalibration(self,l_chroms,d_chrom_lens,):

        '''
Validated human SNP data suggests that the Ti/TV should be ~2.1 genome-wide and ~2.8 in exons (ref ???)
http://www.broadinstitute.org/gsa/wiki/index.php/QC_Methods#SNP_callset_metrics
http://www.broadinstitute.org/gsa/wiki/index.php/Variant_quality_score_recalibration#Ti.2FTv-free_recalibration
http://www.broadinstitute.org/gsa/wiki/images/b/b2/TiTv_free_VQSR.pdf
"For new approach: just cut at 99% sensitivity"
http://www.broadinstitute.org/gsa/wiki/images/a/ac/Ngs_tutorial_depristo_1210.pdf
http://www.broadinstitute.org/gsa/wiki/images/b/bc/Variant_Recalibrator_Sanger_June_2010.pdf
##
http://www.broadinstitute.org/gsa/wiki/images/e/eb/FP_TITV.jpg
'''

        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_variantrecalibration_ApplyRecalibration.html

        T = analysis_type = 'ApplyRecalibration'
        ##fula_20120704/stdout/ApplyRecalibration/ApplyRecalibration.out:    Max Memory :      2187 MB
        ##zulu_20121208/stdout/ApplyRecalibration/ApplyRecalibration.out:    Max Memory :      2210 MB
        ##uganda_20130113/stdout/ApplyRecalibration/ApplyRecalibration.out:    Max Memory :      2233 MB
        memMB = 2900
        ##fula_20120704/stdout/ApplyRecalibration/ApplyRecalibration.out:    CPU time   :   2316.11 sec.
        ##zulu_20121208/stdout/ApplyRecalibration/ApplyRecalibration.out:    CPU time   :   3448.45 sec.
        ##uganda_20130113/stdout/ApplyRecalibration/ApplyRecalibration.out:    CPU time   :   3575.38 sec.
        queue = 'normal'

        if os.path.isfile('ApplyRecalibration.touch'): return

        ##
        ## parse list of vcf input files
        ##
        l_vcfs_in = self.get_fps_in(l_chroms,d_chrom_lens,self.i_UG_size,)

        for mode in ['SNP',]:

            ##
            ## 1) check input existence
            ##
            fp_in_recal = 'out_VariantRecalibrator/VariantRecalibrator.%s.recal' %(mode)
            fp_in_tranches = 'out_VariantRecalibrator/VariantRecalibrator.%s.tranches' %(mode)
    ##        T_prev = self.seqsteps[self.seqsteps.index(analysis_type)-1]
    ##        l_fp_in = self.d_out[T_prev]
            bool_exit = self.check_in('VariantRecalibrator',[fp_in_recal,fp_in_tranches,],)

            ##
            ## 2) touch
            ##
            bool_return = self.touch('%s.%s' %(analysis_type,mode))
            if bool_return == True: continue

            ##
            ## send e-mail to user about choice of ts_filter level
            ##
            address = '%s@sanger.ac.uk' %(pwd.getpwuid(os.getuid())[0],)
            cmd = 'echo "Check your tranches file'
            cmd += ' (%s/%s)' %(os.getcwd(),fp_in_tranches,)
            cmd += ' to see if you are satisfied with the chosen TS level of'
            cmd += ' ${ts_filter_level}" '
            cmd += '| mail -s "%s" ' %(self.project)
            cmd += '%s\n' %(address)
            self.execmd(cmd)

            ##
            ## set output prefix
            ##
            fp_out_prefix = 'out_ApplyRecalibration'
            fp_out_prefix += '/ApplyRecalibration.recalibrated.filtered.%s' %(mode)
            fp_out_vcf = '%s.vcf' %(fp_out_prefix)
            fp_out_idx = '%s.vcf.idx' %(fp_out_prefix)
            if os.path.isfile(fp_out_vcf):
                return

            ##
            ## initialize shell script
            ##
            lines = ['#!/bin/bash']

            ##
            ## --ts-filter-level
            ##
            if self.f_ts_filter_level in [None,'None',]:
                cmd = self.determine_TS_level(fp_in_tranches)
                lines += [cmd]
            else:
                lines += [
                    'ts_filter_level=%f' %(
                        float(self.f_ts_filter_level))]

            ##
            ## initialize GATK java command
            ##
            lines += self.init_cmd(analysis_type,memMB,)

            ## GATKwalker, required, in
##            import glob ## tmp!!! for union1x and union4x!!!
##            l_vcfs_in = glob.glob('out_CombineVariants/*.vcf') ## tmp!!!
            lines += ['--input %s \\' %(vcf) for vcf in l_vcfs_in if os.path.isfile(vcf)]
            ## GATKwalker, required, out
            lines += ['--out %s \\' %(fp_out_vcf,)]
            ## GATKwalker, required, in
            lines += ['--recal_file %s \\' %(fp_in_recal)]
            lines += ['--tranches_file %s \\' %(fp_in_tranches)]
            ## GATKwalker, optional, in
            ## ts_filter_level should correspond to targetTruthSensitivity prior to novelTiTv dropping (see tranches plot)
            ## default is 99.00
    ##        lines += ['--ts_filter_level 99.0 \\']
            lines += ['--ts_filter_level $ts_filter_level \\']
            ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_variantrecalibration_ApplyRecalibration.html#--mode
            lines += ['--mode %s \\' %(mode)]

            lines += self.term_cmd(
                '%s.%s' %(analysis_type,mode),[fp_out_vcf,fp_out_idx,],)

            self.write_shell('shell/ApplyRecalibration.%s.sh' %(mode),lines,)

            ##
            ## execute shell script
            ##
            J = 'AR'
            cmd = self.bsub_cmd(
                '%s.%s' %(analysis_type,mode),J,memMB=memMB,
                std_suffix='%s/%s.%s' %(analysis_type,analysis_type,mode,),)
            self.execmd(cmd)

        return


    def init_cmd(self,analysis_type,memMB,):

        s = 'cmd="'
        ## run GATK
        s += 'java '
        ## set maximum heap size
        s += '-Xmx%im ' %(memMB)
        s += '-jar %s \\' %(self.fp_GATK)
        lines = [s]

        ## CommandLineGATK, required, in
        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_CommandLineGATK.html#--analysis_type
        lines += ['--analysis_type %s \\' %(analysis_type)]
        ## CommandLineGATK, optional, in
        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_CommandLineGATK.html#--reference_sequence
        lines += ['--reference_sequence %s \\' %(self.fp_FASTA_reference_sequence)]

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
        else:
            memMB = 2900 ## 2600 if 250 samples at 4x (helic), 2100 if 120 samples at 4x or 8x (ethiopia) ## 1900=2000-default

        ##
        ## 1) touch
        ##
        bool_return = self.touch(analysis_type)
        if bool_return == True: return

        fp_out = self.d_out['UnifiedGenotyper']

        ## initiate shell script
        lines = ['#!/bin/bash']

        ## commands prior to GATK command
        lines += self.init_UnifiedGenotyper(
            l_chroms,d_chrom_lens,)

        ## initiate GATK command
        lines += self.init_cmd('UnifiedGenotyper',memMB,)

        ## append GATK command options
        lines += self.body_UnifiedGenotyper(fp_out,)

        ## terminate shell script
        lines += self.term_cmd(analysis_type,[fp_out],)

        ## write shell script
        self.write_shell('shell/UnifiedGenotyper.sh',lines,)

        ## execute shell script
        for chrom in l_chroms:
            intervals = int(math.ceil(
                d_chrom_lens[chrom]/float(self.i_UG_size)))
            J = '%s%s[%i-%i]' %('UG',chrom,1,intervals,)
            std_suffix = '%s/%s.%s.%%I' %(T,T,chrom)
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
        ):

        if not std_suffix:
            std_suffix = '%s/%s' %(analysis_type,analysis_type,)

        cmd = 'bsub -J"%s%s" -q %s' %(self.name,J,queue,)
        cmd += ' -G %s' %(self.project)
        cmd += " -M%i000 -R'select[mem>%i] rusage[mem=%i]'" %(
            memMB,memMB,memMB,)
        cmd += ' -o stdout/%s.out' %(std_suffix)
        cmd += ' -e stderr/%s.err' %(std_suffix)
        cmd += ' ./shell/%s.sh' %(analysis_type,)
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
        ## append input files
        ##
        l = []
        for fp_bam in self.fp_bams:
            for s in os.listdir(fp_bam):
                if s[-4:] != '.bam':
                    continue
                ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_CommandLineGATK.html#--input_file
                lines += [' --input_file %s/%s \\' %(
                    fp_bam,s,)]

        ##
        ## CommandLineGATK, optional
        ##

        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_CommandLineGATK.html#--intervals
        lines += [
            '--intervals $CHROMOSOME:$(((${LSB_JOBINDEX}-1)*%i+1))-$posmax \\' %(
                self.i_UG_size)]
        if self.fp_intervals:
            lines += ['--intervals %s \\' %(self.fp_intervals)]
            ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_CommandLineGATK.html#--interval_set_rule
            lines += ['--interval_set_rule INTERSECTION \\']

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
        lines += [' -stand_call_conf 4 \\'] ## ask Deepti if OK
        lines += [' -stand_emit_conf 4 \\'] ## ask Deepti if OK
        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_genotyper_UnifiedGenotyper.html#--genotype_likelihoods_model
        lines += [' --genotype_likelihoods_model BOTH \\'] ## default value SNP
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
        if '2.4' in self.fp_GATK:
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

        return lines


    def term_cmd(self,analysis_type,l_fp_out,):

        if type(l_fp_out) != list:
            print(l_fp_out)
            stop

        ## cont cmd
        lines = [';']
        for fp_out in l_fp_out:
            fp_out = fp_out
            lines += ['echo %s >> %s.touch;' %(fp_out,analysis_type,)]
            lines += ['touch touch/%s;' %(fp_out,)]

        ## write continuation shell script
        ## do not continue as part of previous command
        ## as this will influence CPU statistics
        ## and risk job of hitting CPU walltime
        s = "bsub -R 'select[mem>1500] rusage[mem=1500]' -M1500000 \\\n"
        s += ' -o stdout/rerun.out \\\n'
        s += ' -e stderr/rerun.err \\\n'
        s += ' -G %s \\\n' %(self.project)
        s += ' bash ./rerun_python.sh'
        fd = open('rerun.sh','w')
        fd.write(s)
        fd.close()
        self.execmd('chmod +x rerun.sh')

        s = ' /software/bin/python-2.7.3'
        s += ' %s/GATK_pipeline2.py' %(os.path.dirname(sys.argv[0]))
        for k,v in vars(self.namespace_args).items():
            s += ' --%s %s' %(k,str(v).replace('$','\\$'))
        fd = open('rerun_python.sh','w')
        fd.write(s)
        fd.close()
        self.execmd('chmod +x rerun_python.sh')

        ## cont cmd
        lines += ['bash ./rerun.sh']
        ## term cmd
        lines += ['"']

        lines += ['echo $cmd']
        lines += ['eval $cmd']

        return lines


    def init_UnifiedGenotyper(self,l_chroms,d_chrom_lens,):

        lines = []

        lines += ['\n## parse chromosome from command line']
        lines += ['CHROMOSOME=$1']

        lines += self.bash_chrom2len(l_chroms,d_chrom_lens,)

        ##
        ## do not allow interval to exceed the length of the chromosome
        ## otherwise it will raise an I/O error (v. 1.4-15)
        ##
        lines += ['posmax=$(((${LSB_JOBINDEX}+0)*%i))' %(self.i_UG_size)]
        lines += ['if [ $posmax -gt $LENCHROMOSOME ]']
        lines += ['then posmax=$LENCHROMOSOME']
        lines += ['fi\n']

        fn = 'out_UnifiedGenotyper/'
        fn += 'UnifiedGenotyper.${CHROMOSOME}.${LSB_JOBINDEX}.vcf'

        ## job finished
        lines += ['if [ -s %s.idx ]; then' %(fn)]
        lines += ['exit']
        lines += ['fi\n']

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

        ## http://docs.python.org/2/library/argparse.html#module-argparse
        ## New in version 2.7
        ## optparse deprecated since version 2.7
        parser = argparse.ArgumentParser()

        ##
        ## add arguments
        ##

        parser.add_argument(
            '--bam','--bams','--bamdir','--fp_bams',
            dest='fp_bams',
            help='Path to directory containing improved BAMs',
            type=str,default=None,
            nargs='+',
            required = True,
            )

        parser.add_argument(
            '--GATK','--fp_GATK',
            dest='fp_GATK',
            help='File path to GATK (e.g. /software/varinf/releases/GATK/GenomeAnalysisTK-1.4-15-gcd43f01/GenomeAnalysisTK.jar)',
            type=str,default=None,
            required = True,
            )

        ## wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.3/b37/human_g1k_v37.fasta.gz
        parser.add_argument(
            '--FASTA','--reference','--reference-sequence','--reference_sequence','--fp_FASTA_reference_sequence',
            dest='fp_FASTA_reference_sequence',
            help='File path to reference sequence in FASTA format (e.g. /lustre/scratch111/resources/vrpipe/ref/Homo_sapiens/1000Genomes/human_g1k_v37.fasta)',
            type=str,default=None,
            required = True,
            )

        parser.add_argument(
            '--project','-P','-G',
            dest='project',
            help='Project',
            type=str,default=None,
            required = True,
            )

        parser.add_argument(
            '--arguments','--args','--options','--opts','--fp_arguments',
            dest='fp_arguments',
            type=str,default=None,
            required = False,
            )

        parser.add_argument(
            '--name','--dataset',
            dest='name',
            help='Name of dataset',
            type=str,default=None,
            required = False,
            )

        ##
        ## UnifiedGenotyper
        ##

        parser.add_argument(
            '--i_UG_size',
            dest='i_UG_size',
            help='Size (bp) of divided parts.',
            type=str,default=10*10**6,
            required = False,
            )

        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_genotyper_UnifiedGenotyper.html#--genotyping_mode
        parser.add_argument(
            '--genotyping_mode','-gt_mode',
            dest='genotyping_mode',
            help='Specifies how to determine the alternate alleles to use for genotyping.',
            type=str,default='DISCOVERY',
            required = False,
            )

        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_genotyper_UnifiedGenotyper.html#--output_mode
        parser.add_argument(
            '--output_mode','-out_mode',
            dest='output_mode',
            help='Specifies which type of calls we should output.',
            type=str,default='EMIT_VARIANTS_ONLY',
            required = False,
            )

        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_genotyper_UnifiedGenotyper.html#--alleles
        parser.add_argument(
            '--alleles','-alleles','--fp_alleles',
            dest='fp_alleles',
            help='The set of alleles at which to genotype when --genotyping_mode is GENOTYPE_GIVEN_ALLELES.',
            default='',
            type=str,
            required = False,
            )

        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_CommandLineGATK.html#--intervals
        parser.add_argument(
            '--intervals','-L','--fp_intervals',
            dest='fp_intervals',
            help='Additionally, one may specify a rod file to traverse over the positions for which there is a record in the file (e.g. -L file.vcf).',
            type=str,default=None,
            required = False,
            )

        ##
        ## VariantRecalibrator resources
        ##

        parser.add_argument(
            '--resources','--VariantRecalibrator','--fp_resources_SNP','--fp_resources',
            dest='fp_resources_SNP',
            help='File path to a file with -resource lines to append to GATK VariantRecalibrator',
            type=str,default=None,
            required = False,
            )

        parser.add_argument(
            '--resources_INDEL','--VariantRecalibrator_INDEL','--fp_resources_INDEL',
            dest='fp_resources_INDEL',
            help='File path to a file with -resource lines to append to GATK VariantRecalibrator',
            type=str,default=None,
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
            '--dbsnp','--fp_vcf_dbsnp',
            dest='fp_vcf_dbsnp',
            help=s_help,
            type=str,default=None,
            required = True,
            )

        ##
        ## ApplyRecalibration
        ##

        parser.add_argument(
            '--f_ts_filter_level','--ts','--f_ApplyRecalibration_ts_filter_level','--ts_filter_level',
            dest='f_ts_filter_level',
            help='',
            type=float,default=None,
            required = False,
            )

        ##
        ## BEAGLE
        ##
        parser.add_argument(
            '--beagle','--BEAGLE','--BEAGLEjar','--fp_software_beagle',
            dest='fp_software_beagle',
            help='File path to BEAGLE .jar file (e.g. /nfs/team149/Software/usr/share/beagle_3.3.2.jar)',
            type=str,default=None,
            required = True,
            )

        parser.add_argument(
            '--i_BEAGLE_size',
            dest='i_BEAGLE_size',
            help='Size (Mbp) of divided parts.',
            type=str,default=2, ## CPU bound (2Mbp=22hrs,3090MB) with lowmem option...
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
            type=str,default=None,
            required = True,
            )

        s_help = 'Markers file to be divided (e.g.'
        s_help += ' ALL.chr$CHROM.phase1_release_v3.20101123.filt.renamed.markers'
        parser.add_argument(
            '--fp_BEAGLE_markers',
            dest='fp_BEAGLE_markers',
            help=s_help,
            type=str,default=None,
            required = True,
            )

        parser.add_argument(
            '--i_BEAGLE_nsamples',
            dest='i_BEAGLE_nsamples',
            help='',
            type=str,default=20,
            required = False,
            )

        ##
        ## IMPUTE2
        ##
        parser.add_argument(
            '--impute2','--IMPUTE2','--IMPUTE2jar','--fp_software_impute2','--fp_software_IMPUTE2',
            dest='fp_software_IMPUTE2',
            help='File path to IMPUTE2 executable (e.g. /software/hgi/impute2/v2.2.2/bin/impute2)',
            type=str,default=None,
            required = True,
            )

        parser.add_argument(
            ## http://mathgen.stats.ox.ac.uk/impute/impute_v2.html
            '--hap','--impute2-hap','--fp_impute2_hap',
            dest='fp_impute2_hap',
            help='Hap files used by IMPUTE2 (e.g. /nfs/t149_1kg/ALL_1000G_phase1integrated_v3_impute/ALL_1000G_phase1integrated_v3_chr$CHROM_impute.hap.gz)',
            type=str,default=None,
            required = True,
            )

        parser.add_argument(
            ## http://mathgen.stats.ox.ac.uk/impute/impute_v2.html
            '--legend','--impute2-legend','--fp_impute2_legend',
            dest='fp_impute2_legend',
            help='Legend files used by IMPUTE2 (e.g. /nfs/t149_1kg/ALL_1000G_phase1integrated_v3_impute/ALL_1000G_phase1integrated_v3_chr$CHROM_impute.legend.gz)',
            type=str,default=None,
            required = True,
            )

        parser.add_argument(
            ## http://mathgen.stats.ox.ac.uk/impute/impute_v2.html
            '--map','--impute2-map','--fp_impute2_map',
            dest='fp_impute2_map',
            help='Map files used by IMPUTE2 (e.g. /nfs/t149_1kg/ALL_1000G_phase1integrated_v3_impute/genetic_map_chr$CHROM_combined_b37.txt)',
            type=str,default=None,
            required = True,
            )

        parser.add_argument(
            '--i_IMPUTE2_size',
            dest='i_IMPUTE2_size',
            help='Size (bp) of divided parts.',
            type=int,default=5000000, ## memory bound (5Mbp=6GB,11hrs)...
            required = False,
            )

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
            'UnifiedGenotyper':'out_UnifiedGenotyper/UnifiedGenotyper.$CHROMOSOME.${LSB_JOBINDEX}.vcf',
            'ProduceBeagleInput':'out_ProduceBeagleInput/ProduceBeagleInput.bgl',
            'BEAGLE':'out_BEAGLE/$CHROMOSOME/$CHROMOSOME.${LSB_JOBINDEX}.like',
            'IMPUTE2':'out_IMPUTE2/$CHROMOSOME/$CHROMOSOME.${LSB_JOBINDEX}.gen',
            }

        self.d_in = {
            'ProduceBeagleInput':'out_ApplyRecalibration/ApplyRecalibration.recalibrated.filtered.vcf',
            'BEAGLE':'out_ProduceBeagleInput/ProduceBeagleInput.bgl',
            'IMPUTE2':'in_IMPUTE2/$CHROMOSOME.gen',
            }

        self.verbose = True

        return


    def IMPUTE2_without_BEAGLE(self,l_chroms,d_chrom_lens,):

        ## tmp function!!!

        if not os.path.isdir('stderr/IMPUTE2_without_BEAGLE'):
            os.mkdir('stderr/IMPUTE2_without_BEAGLE')
        if not os.path.isdir('stdout/IMPUTE2_without_BEAGLE'):
            os.mkdir('stdout/IMPUTE2_without_BEAGLE')
        if not os.path.isdir('touch/out_IMPUTE2_without_BEAGLE'):
            os.mkdir('touch/out_IMPUTE2_without_BEAGLE')
        if not os.path.isdir('out_IMPUTE2_without_BEAGLE'):
            os.mkdir('out_IMPUTE2_without_BEAGLE')
            for chrom in l_chroms:
                os.mkdir('out_IMPUTE2_without_BEAGLE/%s' %(chrom))

        if not os.path.isdir('in_IMPUTE2_without_BEAGLE'):
            os.mkdir('in_IMPUTE2_without_BEAGLE')

            cmd = 'cat out_ProduceBeagleInput/ProduceBeagleInput.bgl '
            cmd += ''' | awk 'BEGIN{FS=":"} '''
            cmd += 'NR>1{print $1,$2}'
            cmd += "'"
            cmd += " | awk '"
            cmd += '{print $1":"$2,$1":"$2,$0}'
            cmd += "'"
            cmd += ' | cut -d " " -f 3 --complement'
            cmd += ''' | awk 'BEGIN{FS=":"} '''
            cmd += '{print>"in_IMPUTE2_without_BEAGLE/"$1".gen"'
            cmd += "}'"
            print(cmd)
            self.execmd(cmd)

        if True:

            queue = 'normal'
            fp_in = 'in_IMPUTE2_without_BEAGLE/$CHROMOSOME.gen'
            fp_out = 'out_IMPUTE2_without_BEAGLE/$CHROMOSOME/$CHROMOSOME.${LSB_JOBINDEX}.gen'

            ## write shell script to list/memory
            lines = self.IMPUTE2_write_shell_script(fp_in,fp_out)

            ## term cmd
            lines += self.term_cmd('IMPUTE2_without_BEAGLE',[fp_out],)

            ## write shell script to file
            self.write_shell('shell/IMPUTE2_without_BEAGLE.sh',lines,)

            ##
            ## execute shell script
            ##
            for chrom in l_chroms:
                print('bsub IMPUTE2 wo BEAGLE', chrom)
                self.bsub_IMPUTE2('IMPUTE2_without_BEAGLE',chrom,d_chrom_lens)

        return


    def IMPUTE2_unite(self,l_chroms,d_chrom_lens,):

        print('IMPUTE2_unite')

##        ##
##        ## 1) check input existence
##        ##
##        l_fp_in = []
##        for chrom in l_chroms:
##            index_max = int(math.ceil(
##                    d_chrom_lens[chrom]/float(self.i_IMPUTE2_size)))
##            for index in xrange(1,index_max+1):
##                ## use summary instead of genotype file, because of
##                ## --read 0 SNPs in the analysis interval+buffer region
##                l_fp_in += ['out_IMPUTE2/%s/%s.%i.gen_summary' %(
##                    chrom,chrom,index)]
##        self.check_in('IMPUTE2',l_fp_in)

        ##
        ## 2) touch
        ##
        bool_return = self.touch('IMPUTE2_unite')
        if bool_return == True: return

        ##
        ## 3) file I/O checks
        ##
        for chrom in l_chroms:
            print('chrom', chrom)
            self.IMPUTE2_fileIO_checks(chrom)

        ##
        ## 4) IMPUTE output concatenation
        ##
        for chrom in l_chroms:
            index_max = int(math.ceil(
                    d_chrom_lens[chrom]/float(self.i_IMPUTE2_size)))
            cmd = 'cat'
            for index in range(1,index_max+1):
                fp_gen = 'out_IMPUTE2/%s/%s.%i.gen' %(chrom,chrom,index)
                if not os.path.isfile(fp_gen): continue
                cmd += ' %s' %(fp_gen)
            cmd += ' > out_IMPUTE2/%s.gen' %(chrom)
            self.execmd(cmd)

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


## DISK SPACE

## todo20130428 - instead of using "xxx.renamed.markers", just rename on the fly
## out: chrom:pos pos alleleA alleleB\n
