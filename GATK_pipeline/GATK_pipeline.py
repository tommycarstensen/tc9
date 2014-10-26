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
import Bio
from Bio.bgzf import BgzfWriter


## README

## http://www.broadinstitute.org/gatk/about#typical-workflows
## http://www.broadinstitute.org/gatk/guide/topic?name=best-practices

class main():

    def main(self):

        if self.args.AR_input:
            self.run_ApplyRecalibration()

##        self.HaplotypeCaller()
##        self.CombineGVCFs()
##        self.GenotypeGVCFs()
##
##        ## 1000G_phase1.snps.high_confidence.b37.vcf.gz only contains chromosomes 1-22 and X
##        self.VariantRecalibrator()

        self.ApplyRecalibration()

##        ## phased imputation reference panels not available for chrY
##        self.chroms.remove('Y')
##        ## Exception in thread "main" java.lang.ArrayIndexOutOfBoundsException: -2
##        self.chroms.remove('X')
        self.BEAGLE4()
##        self.BEAGLE3()
##        self.unite_BEAGLE3()

        return


    def HaplotypeCaller(self):

        T = analysis_type = 'HaplotypeCaller'
        queue = 'long'  # 4x split per sample
        memMB = 3900
#        queue = 'basement'  # >4x split per sample
        queue = 'normal'  # 4x split per chromosome

        ## Create folders.
        self.mkdir('tmp/')

        ## touch/lock
        if self.touch(analysis_type):
            return

        ## write shell script
        self.shell_HC(analysis_type, memMB)

        ## execute shell script
        for chrom in self.chroms:

            for bam in self.bams:
                basename = os.path.splitext(os.path.basename(bam))[0]
                LSF_affix = '{}/{}/{}'.format(T, chrom, basename)

                ## Skip bam and chromosome if output was generated.
                if os.path.isfile(
                    'out_{}/{}/{}.vcf.gz.tbi'.format(T, chrom, basename)):
                    continue

                ## Skip bam and chromosome if output is being generated.
                if os.path.isfile('LSF/{}.err'.format(LSF_affix)):
                    ## file changed within the past 5 minutes?
                    if time.time()-os.path.getmtime(
                        'LSF/{}.err'.format(LSF_affix)) < 300:
                        continue
                    else:
                        os.remove('LSF/{}.err'.format(LSF_affix))
                if os.path.isfile('LSF/{}.out'.format(LSF_affix)):
                    os.remove('LSF/{}.out'.format(LSF_affix))

                self.mkdir('LSF/{}'.format(LSF_affix))
                self.mkdir('out_{}/{}/'.format(T, chrom))

                J = '{} {}'.format('HC', basename)
                cmd = self.bsub_cmd(
                    analysis_type, J, LSF_affix=LSF_affix,
                    memMB=memMB, queue=queue, num_threads=self.nct,
                    bam=bam, chrom=chrom)
                self.execmd(cmd)

        return


    def CombineGVCFs(self):

        '''Merge gVCFs prior to GenotypeGVCFs'''

        analysis_type = T = 'CombineGVCFs'
        memMB = 9900
        queue = 'long'
        queue = 'basement'

        ## write shell script
        self.shell_CombineGVCFs(T, memMB)

        bool_exit = False
        for chrom in self.chroms:

            ## 1) check input existence / check that previous jobs finished
            l_vcfs_in = [
                'out_HaplotypeCaller/{}/{}.vcf.gz'.format(
                    chrom, os.path.splitext(os.path.basename(bam))[0])
                for bam in sorted(self.bams)]
            if self.check_in(
                'HaplotypeCaller', ['{}.tbi'.format(vcf) for vcf in l_vcfs_in],
                'touch/HaplotypeCaller.touch'):
                bool_exit = True
                continue

            ## 2) check output existence / check that job did not start
            if self.touch('{}.{}'.format(analysis_type, chrom)):
                continue

            self.mkdir('lists')
            l_combined = []
            for fn_list in glob.glob('lists/{}.{}.*.list'.format(T, chrom)):
                with open(fn_list) as f:
                    l_combined += f.read().rstrip().split('\n')

            l_vcfs_in = list(sorted(list(set(l_vcfs_in)-set(l_combined))))

            for i, vcf in enumerate(
                l_vcfs_in,
                self.gVCF_limit*len(glob.glob(
                    'lists/{}.{}.*.list'.format(T, chrom)))):
                if i%self.gVCF_limit == 0:
                    fn_out = 'lists/{T}.{chrom}.{i}.list'.format(
                        T=T, chrom=chrom, i=i//self.gVCF_limit)
                    assert not os.path.isfile(fn_out)
                    fd_out = open(fn_out, 'w')
                fd_out.write('{}\n'.format(vcf))
            fd_out.close()

            self.mkdir('LSF/{}'.format(T))
            for i in range(len(glob.glob(
                'lists/{}.{}.*.list'.format(T, chrom)))):
                ## skip if job initiated
                if os.path.isfile(
                    'out_{}/{}/{}.vcf.gz'.format(T, chrom, i)):
                    continue
                cmd = self.bsub_cmd(
                    T, 'CgVCFs.{}.{}'.format(chrom, i),
                    memMB=memMB, queue=queue,
                    LSF_affix='{}/{}.{}'.format(T, chrom, i)
                    )
                cmd += ' {} {}'.format(chrom, i)
                self.execmd(cmd)

        if bool_exit == True:
            sys.exit()

        return


    def GenotypeGVCFs(self):

        '''Convert gVCFs to VCFs'''

        analysis_type = T = 'GenotypeGVCFs'
        memMB = 4900
        queue = 'basement'

        ## write shell script
        self.shell_GenotypeGVCFs(T, memMB)
        self.mkdir('LSF/{}'.format(T))

        bool_exit = False
        for chrom in self.chroms:

            ## 1) check input existence / check that previous jobs finished
            l_vcfs_in = [
                'out_CombineGVCFs/{}/{}.vcf.gz'.format(chrom, i)
                for i in range(len(glob.glob(
                    'lists/CombineGVCFs.{}.*.list'.format(chrom))))]
            if self.check_in(
                'CombineGVCFs', ['{}.tbi'.format(vcf) for vcf in l_vcfs_in],
                'touch/CombineGVCFs.touch'):
                ## continue loop over chromosomes
                bool_exit = True
                continue

            ## 2) check output existence / check that job did not start
            if self.touch('{}.{}'.format(analysis_type, chrom)):
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
                    LSF_affix='{}/{}'.format(T, chrom),
                    num_threads=4,
                    )
                cmd += ' {}'.format(chrom)
                self.execmd(cmd)

        if bool_exit:
            sys.exit()

        return


    def VariantRecalibrator(self):

        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_variantrecalibration_VariantRecalibrator.html

        T = analysis_type = 'VariantRecalibrator'
        memMB = 19900
        queue = 'yesterday'
        num_threads = 4

        ##
        ## 1) check input existence (vcf)
        ##
        l_vcfs_in = [
            'out_GenotypeGVCFs/{}.vcf.gz'.format(chrom) for chrom in self.chroms]
        if self.check_in(
            'GenotypeGVCFs', ['{}.tbi'.format(vcf) for vcf in l_vcfs_in],
            'touch/GenotypeGVCFs.touch'):
            sys.exit(0)

        d_resources = {'SNP':self.fp_resources_SNP, 'INDEL':self.fp_resources_INDEL,}

        if not os.path.isdir('LSF/VariantRecalibrator'):
            os.mkdir('LSF/VariantRecalibrator')
        if not os.path.isdir('out_VariantRecalibrator'):
            os.mkdir('out_VariantRecalibrator')

        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_variantrecalibration_VariantRecalibrator.html#--mode
        for mode in ['SNP', 'INDEL',]:

            memMB = {'SNP':20900, 'INDEL':8900}[mode]

            ## 2) touch / check output
            bool_continue = self.touch('{}.{}'.format(analysis_type,mode))
            if bool_continue == True:
                continue

            ## Define file paths.
            fp_tranches = 'out_VariantRecalibrator/{}.tranches'.format(mode)
            fp_recal = 'out_VariantRecalibrator/{}.recal'.format(mode)

            lines = ['out={}'.format(fp_tranches)]

            ## Initiate GATK walker.
            lines += self.init_GATK_cmd(analysis_type,memMB,)

            lines += [' --num_threads {} \\'.format(num_threads)]

            ## required, in
            lines += [' --input {} \\'.format(vcf) for vcf in l_vcfs_in]
            ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_variantrecalibration_VariantRecalibrator.html#--use_annotation
            ## http://gatkforums.broadinstitute.org/discussion/2805/howto-recalibrate-variant-quality-scores-run-vqsr
            if mode == 'SNP':
                lines += [
                    ' --use_annotation DP \\',
                    ' --use_annotation QD \\',
                    ' --use_annotation FS \\',
                    ' --use_annotation MQRankSum \\',
                    ' --use_annotation ReadPosRankSum \\',
##                    ' --use_annotation InbreedingCoeff \\',
                    ]
            elif mode == 'INDEL':
##                lines += [' -an InbreedingCoeff ',]
                lines += [' -an QD -an DP -an FS -an ReadPosRankSum -an MQRankSum \\',]

            ##
            ## required, out
            ##
            lines += [' --recal_file {} \\'.format(fp_recal)]
            lines += [' --tranches_file {} \\'.format(fp_tranches)]

            ##
            ## Optional Parameters.
            ##
            ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_variantrecalibration_VariantRecalibrator.html#--mode
            lines += [' --mode {} \\'.format(mode)]

            ## http://gatkforums.broadinstitute.org/discussion/1259/what-vqsr-training-sets-arguments-should-i-use-for-my-specific-project
            if mode == 'INDEL':
                ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_variantrecalibration_VariantRecalibrator.html#--maxGaussians
                lines += [' --maxGaussians 4 \\']  # default 8

            ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_variantrecalibration_VariantRecalibrator.html#--resource
            fd = open(d_resources[mode], 'r')
            lines_resources = fd.readlines()
            fd.close()
            lines += [' {} \\'.format(line.strip()) for line in lines_resources]

            ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_variantrecalibration_VariantRecalibrator.html#--TStranche
            l_TStranches = []
    ##        l_TStranches += [99.70+i/20. for i in range(6,0,-1,)]
            l_TStranches += [99+i/10. for i in range(10,0,-1,)]
            l_TStranches += [90+i/2. for i in range(18,-1,-1,)]
##            l_TStranches += [70+i for i in range(19,-1,-1,)]
            s_TStranches = ''
            for TStranche in l_TStranches:
                s_TStranches += '--TStranche {:.1f} '.format(TStranche)
            lines += [' {} \\'.format(s_TStranches)]

            ##
            ## GATKwalker, optional, out
            ##
            lines += [' --rscript_file out_{}/{}.plots.R \\'.format(T, mode,)]

            ## Terminate command and rerun pipeline.
            lines += self.term_cmd(
                '{}.{}'.format(analysis_type,mode),[fp_tranches,fp_recal,],)

            self.write_shell('shell/{}.{}.sh'.format(analysis_type,mode,), lines,)

            J = 'VR.{}'.format(mode)
            cmd = self.bsub_cmd(
                '{}.{}'.format(analysis_type,mode), J, memMB=memMB, queue=queue,
                LSF_affix='{}/{}'.format(analysis_type,mode,),
                num_threads=num_threads)
            self.execmd(cmd)

        return


    def skip_header(self, fd):

        for line in fd:
            if line[0] == '#':
                continue
            yield line

        return


    def hook_compressed_text(self, filename, mode):

        ext = os.path.splitext(filename)[1]
        if ext == '.gz':
            f = gzip.open(filename, mode + 't')
        else:
            f = open(filename, mode)

        return f


    def parse_recal(self, fd, pattern):

        line = next(self.skip_header(fd))
        l = line.rstrip().split('\t')
        chrom = l[0]
        pos = l[1]
        VQSLod = float(re.match(pattern, l[7]).group(1))

        return chrom, pos, VQSLod


    def parse_sources(self):

        d_sources = {}
        for mode in ('SNP', 'INDEL'):
            with open('out_VariantRecalibrator/{}.recal'.format(mode)) as f:
                for line in f:
                    if line.split('=')[0] != '##GATKCommandLine':
                        continue
                    sources = re.findall(
                        r'\([^\[\]()]*\[\([^\[\]()]+source=([\w./]+)', line)
                    d_sources[mode] = sources
                    break

        return d_sources


    def parse_minVQSLods(self):

        d_ts_filter_level = {
            'SNP':self.ts_filter_level_SNP,
            'INDEL':self.ts_filter_level_INDEL,}
        d_minVQSLod = {}
        for mode in ('SNP','INDEL'):
            with open(
                'out_VariantRecalibrator/{}.tranches'.format(mode)) as f:
                for line in f:
                    if line[0] == '#':
                        continue
                    l = line.split(',')
                    if l[0] == 'targetTruthSensitivity':
                        index = l.index('minVQSLod')
                        continue
                    if float(l[0]) < d_ts_filter_level[mode]:
                        continue
                    d_minVQSLod[mode] = float(l[index])
                    break

        return d_minVQSLod


    def bsub_ApplyRecalibration(self):

        d_sources = self.parse_sources()
        assert d_sources['SNP'] == d_sources['INDEL']

        if not os.path.isdir('LSF/ApplyRecalibration'):
            os.makedirs('LSF/ApplyRecalibration')

        for source_SNP, source_INDEL in zip(
            d_sources['SNP'], d_sources['INDEL']):
            chrom = os.path.basename(source_SNP).split('.')[0]
            assert chrom in [str(i) for i in range(1,23)]+['X', 'Y', 'MT']
            s = ''
            s += 'bsub -G {} '.format(self.project)
            s += ' -o LSF/ApplyRecalibration/{}.out'.format(chrom)
            s += ' -e LSF/ApplyRecalibration/{}.err'.format(chrom)
            s += ' -J AR{}'.format(chrom)
            self.args.AR_input = source_SNP
            s += self.args_to_command_line()
            subprocess.call(s, shell=True)
        sys.exit()

        return


    def ApplyRecalibration(self):

        '''Does the same as GATK ApplyRecalibration,
except does *not* emit INDELs untouched in the output VCF
and requires less than 100MB of memory'''

        analysis_type = T = 'ApplyRecalibration'
##        num_threads = 1

        ## check input existence
        for mode in ('SNP', 'INDEL'):
            if self.check_in(
                'VariantRecalibrator',
                [
                    'out_VariantRecalibrator/{}.recal'.format(mode),
                    'out_VariantRecalibrator/{}.tranches'.format(mode),],
                'touch/VariantRecalibrator.{}.touch'.format(mode)):
                sys.exit()

        ## check output existence
        if self.touch(T):
            return

        self.bsub_ApplyRecalibration()

        return


    def run_ApplyRecalibration(self):

        d_minVQSLod = self.parse_minVQSLods()

        chrom = os.path.basename(self.args.AR_input).split('.')[0]

        pattern = re.compile(r'.*VQSLOD=([-\d.]*)')
        fp_out = 'out_ApplyRecalibration/{}.vcf.gz'.format(chrom)
        with open('out_VariantRecalibrator/SNP.recal') as fd_recal_SNP, \
             open('out_VariantRecalibrator/INDEL.recal') as fd_recal_INDEL, \
             gzip.open(self.args.AR_input, 'rt') as fd_source, \
             BgzfWriter(fp_out, 'wb') as fd_out:
            chrom_SNP = chrom_INDEL = None
            line_VCF = next(self.skip_header(fd_source))
            chrom_VCF, pos_VCF = line_VCF.split('\t', 2)[:2]
            while chrom_SNP != chrom_VCF:
                chrom_SNP, pos_SNP, VQSLod_SNP = self.parse_recal(
                    fd_recal_SNP, pattern)
            while chrom_INDEL != chrom_VCF:
                chrom_INDEL, pos_INDEL, VQSLod_INDEL = self.parse_recal(
                    fd_recal_INDEL, pattern)
            while True:
                print(chrom, pos_VCF, file=sys.stderr)
                if pos_VCF == pos_INDEL:
                    assert chrom_VCF == chrom_INDEL
                    if VQSLod_INDEL >= d_minVQSLod['INDEL']:
                        print(line_VCF, end='', file=fd_out)
                    try:
                        line_VCF = next(self.skip_header(fd_source))
                    except StopIteration:
                        break
                    chrom_VCF, pos_VCF = line_VCF.split('\t', 2)[:2]
                    try:
                        chrom_INDEL, pos_INDEL, VQSLod_INDEL = self.parse_recal(
                            fd_recal_INDEL, pattern)
                    except StopIteration:
                        continue
                    continue
                else:
                    assert pos_VCF == pos_SNP
                    assert chrom_VCF == chrom_SNP
                    if VQSLod_SNP >= d_minVQSLod['SNP']:
                        print(line_VCF, end='', file=fd_out)
                    try:
                        line_VCF = next(self.skip_header(fd_source))
                    except StopIteration:
                        break
                    chrom_VCF, pos_VCF = line_VCF.split('\t', 2)[:2]
                    try:
                        chrom_SNP, pos_SNP, VQSLod_SNP = self.parse_recal(
                            fd_recal_SNP, pattern)
                    except StopIteration:
                        continue
                    continue
                continue
            pass

        ## index bgz output
        subprocess.call('tabix -p vcf {}'.format(fp_out), shell=True)
        ## confirm process has run to completion by writing to file
        with open('touch/ApplyRecalibration.touch', 'a') as f:
            f.write('{}.tbi\n'.format(fp_out))

        return


    def BEAGLE4(self):

        ## http://faculty.washington.edu/browning/beagle

        memMB = 12900  # todo: move to argparse
        window = 50000  # todo: move to argparse
        queue = 'basement'  # todo: move to argparse

        l_chroms = []
        for source_SNP in self.parse_sources()['SNP']:
            l_chroms.append(os.path.basename(source_SNP).split('.')[0])

        ## 1) check input existence
        if self.check_in(
            'ApplyRecalibration',
            [
                'out_ApplyRecalibration/{}.vcf.gz.tbi'.format(chrom)
                for chrom in l_chroms],
            'touch/ApplyRecalibration.touch'):
            sys.exit()

        ## 2) check that process didn't start or end
        if self.touch('BEAGLE'):
            return

        ## initiate shell script
        lines = ['#!/bin/bash\n']
        ## parse chromosome from command line
        lines += ['chrom=$1']
        lines += ['pos1=$2']
        lines += ['pos2=$3']
        lines += ['out=out_BEAGLE/$chrom/${LSB_JOBINDEX}']
        lines += ['mkdir -p $(dirname $out)']
        ## exit if output already exists
        lines += ['if [ -s $out.gprobs.gz ]; then exit; fi']
        ## initiate BEAGLE
        lines += ['{} \\'.format(
            self.init_java(self.fp_software_beagle, memMB))]
        ## Arguments for specifying data
        lines += [' gl=out_ApplyRecalibration/$chrom.vcf.gz \\']
        if self.ped:
            lines += [' ped={} \\'.format(self.ped)]
        lines += [' out=$out \\']
##        lines += [' excludemarkers={} \\'.format(excludemarkers)]
        lines += [' chrom=$chrom:$pos1-$pos2 \\']
        ## Other arguments
        lines += [' nthreads=1 \\']
        lines += [' window={} \\'.format(window)]
        lines += [' overlap=3000 \\']
        lines += [' gprobs=true \\']
        lines += [' usephase=false \\']
        lines += [' seed=-99999 \\']
        lines += [' singlescale=1.0 \\']
        lines += [' duoscale=1.0 \\']
        lines += [' trioscale=1.0 \\']
        lines += [' burnin-its=5 \\']
        lines += [' phase-its=5 \\']
        lines += [' impute-its=5 \\']
        ## Advanced options not recommended for general use
        lines += [' nsamples=4 \\']
        lines += [' buildwindow=1200 \\']
        ## term cmd
        lines += self.term_cmd('BEAGLE', ['$out.gprobs.gz'],)
        ## write shell script
        if not os.path.isfile('shell/BEAGLE.sh'):
            self.write_shell('shell/BEAGLE.sh',lines,)

        if not os.path.isdir('LSF/BEAGLE'):
            os.mkdir('LSF/BEAGLE')

        ##
        ## execute shell script
        ##
        for chrom in l_chroms:
            print('BEAGLE chrom', chrom)
            fd_vcf = gzip.open(
                'out_ApplyRecalibration/{}.vcf.gz'.format(chrom), 'rt')
            cnt = 0
            pos_prev = None
            for line in fd_vcf:
                if line[0] == '#':
                    continue
                l = line.split('\t',2)
                chrom = l[0]
                pos = int(l[1])
                cnt += 1
                if cnt == 1 or pos_prev == None:
                    pos1 = pos
                elif cnt % window == 0:
                    pos2 = pos
                    index = cnt//window
                    self.bsub_BEAGLE(chrom, pos1, pos2, index, memMB, queue)
                    print(chrom, ':', pos1, '-', pos2, index, cnt)
                    pos = pos_prev = None
                else:
                    pass
                pos_prev = pos
                continue

            pos2 = pos
            index = (cnt//window)+1
            stop3
            self.bsub_BEAGLE(chrom, pos1, pos2, index, memMB, queue)

        return


    def bsub_BEAGLE(self, chrom, pos1, pos2, index, memMB, queue):

        ## finished?
        fn_out = 'out_BEAGLE/{}/{}.{}.gprobs.gz'.format(
            chrom,chrom,index)
        if os.path.isfile(fn_out):
            return

        ## started and running?
        fn = 'LSF/BEAGLE/{}.{}.out'.format(chrom,index)
        if os.path.isfile(fn):
            if os.path.getsize(fn):
                ## running?
                with open(fn) as f:
                    if 'iteration' in f.readlines()[-1]:
                        return

        print(chrom,index)

        J = '{}.{}[{}-{}]'.format('BEAGLE',chrom,index,index,)
        LSF_affix = '{}/{}.%%I'.format('BEAGLE',chrom)
        cmd = self.bsub_cmd(
            'BEAGLE', J, memMB=memMB, LSF_affix=LSF_affix,
            chrom=chrom, queue=queue, pos1=pos1, pos2=pos2)
        os.system(cmd)

        return


    def shell_CombineGVCFs(self, T, memMB):

        lines = ['#!/bin/bash\n']
        lines += ['chrom=$1']
        lines += ['i=$2']
        lines += ['out=out_{}/$chrom/$i.vcf.gz'.format(T)]
        lines += ['## exit if job started']
        lines += ['if [ -s $out ]; then exit; fi\n']

        lines += self.init_GATK_cmd(T, memMB,)
        lines += [' -L $chrom \\']
        lines += [' -V lists/{}.$chrom.$i.list \\'.format(T)]
        lines += [' -o $out \\']

        ## terminate shell script
        lines += self.term_cmd(T, ['$out.tbi'],)

        ## write shell script
        self.write_shell('shell/{}.sh'.format(T), lines,)

        return


    def shell_GenotypeGVCFs(self, T, memMB):

        lines = ['#!/bin/bash\n']
        lines += ['chrom=$1']
        lines += ['out=out_{}/$chrom.vcf.gz'.format(T)]
        lines += ['## exit if job started']
        lines += ['if [ -s $out ]; then exit; fi\n']
        lines += ['## exit if job finished']
        lines += ['if [ -s $out.tbi ]; then exit; fi\n']

        lines += self.init_GATK_cmd(T, memMB,)
        lines += [' -L $chrom \\']
        lines += [' -V lists/{}.$chrom.list \\'.format(T)]
        lines += [' -o out_{}/$chrom.vcf.gz \\'.format(T)]
        lines += [' -nt 4 \\'.format(T)]
        lines += [' --annotation InbreedingCoeff \\']  # default
        lines += [' --annotation FisherStrand \\']  # default
        lines += [' --annotation QualByDepth \\']  # default
        lines += [' --annotation ChromosomeCounts \\']  # default
        lines += [' --annotation GenotypeSummaries \\']  # default
        lines += [' --annotation MappingQualityRankSumTest \\']
        lines += [' --annotation ReadPosRankSumTest \\']

        ## terminate shell script
        lines += self.term_cmd(T, ['$out.tbi'])

        ## write shell script
        self.write_shell('shell/{}.sh'.format(T), lines,)

        return



    def execmd(self,cmd):

        print(cmd)
        subprocess.call(cmd,shell=True)

        return


    def BEAGLE3_fileIO_checks(self,chrom):

        ##
        ## check that gprobs is identical to markers
        ##
        cmd = "awk '{print $1}' in_BEAGLE/{}/{}.*.markers" %(chrom,chrom)
        cmd += ' | sort -u > markers{}'.format(chrom)
        self.execmd(cmd)

        if os.path.isfile('gprobs{}'.format(chrom)):
            os.remove('gprobs{}'.format(chrom))
        for file in glob.glob(
            'out_BEAGLE/{}/{}.*.like.gprobs.gz'.format(chrom,chrom)):
            cmd = "zcat {}" %(file)
            cmd += " | awk 'FNR>1{print $1}'"
            cmd += ' >> gprobs{}'.format(chrom)
            self.execmd(cmd)

        self.execmd('sort -u gprobs{} -o gprobs{}'.format(chrom,chrom))

        cmd = 'comm -3 gprobs{} markers{}'.format(chrom,chrom)
        i = int(os.popen('{} | wc -l'.format(cmd)).read())
        if i > 0:
            print(os.popen('{} | head'.format(cmd)).readlines())
            print(cmd)
            sys.exit()

        os.remove('markers{}'.format(chrom))
        os.remove('gprobs{}'.format(chrom))

        return


    def BEAGLE3_unite(self,chrom):

        ##
        ##
        ##
        fd = open('BEAGLE_divide_indexes.txt', 'r')
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
        fp_out = 'out_BEAGLE/{}.gprobs.gz'.format(chrom,)
        for index in range(1,max(d_index2pos.keys())+1,):
            fp_in = 'out_BEAGLE/{}/{}.{}.like.gprobs.gz'.format(chrom,chrom,index,)
            if index == 1:
                if not os.path.isfile(fp_out):
                    cmd = 'zcat {} | head -n1 | gzip > {}'.format(fp_in,fp_out,)
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
            cmd = 'zcat {}'.format(fp_in)
            cmd += " | awk 'NR>1{pos=int(substr($1,%i));" %(len(chrom)+2)
            cmd += " if(pos>%i&&pos<=%i) print $0}'" %(pos1,pos2)
            cmd += ' | gzip >> {}'.format(fp_out,)
            self.execmd(cmd)

        return


    def unite_BEAGLE3(self):

        ##
        ## 1) check input existence
        ##
        l_fp_in = self.BEAGLE3_parse_output_files()
        if self.check_in('BEAGLE', l_fp_in, 'touch/BEAGLE.touch'):
            sys.exit()

        for chrom in self.chroms:
            print(chrom)
            fp = 'out_BEAGLE/{}.gprobs.gz'.format(chrom)
            if os.path.isfile(fp):
                continue
            self.BEAGLE3_fileIO_checks(chrom)
            self.BEAGLE3_unite(chrom)

        return


    def BEAGLE3_parse_output_files(self,):

        ## open file
        fd = open('BEAGLE_divide_indexes.txt', 'r')
        ## read all lines into memory
        lines = fd.readlines()
        ## close file
        fd.close()

        l_fp_in = []
        for line in lines:
            l = line.strip().split()
            chrom = l[0]
            if not chrom in self.chroms:
                continue
            index = int(l[1])
            fp_in = 'out_BEAGLE/{}/'.format(chrom,)
            fp_in += '{}.{}.like.gprobs.gz'.format(chrom,index,)
            l_fp_in += [fp_in]

        return l_fp_in


    def write_shell(self,fp,lines,):

        self.mkdir(fp)

        if type(lines) != list:
            print(type(lines))
            stop

        s = '\n'.join(lines)+'\n\n'
        fd = open(fp, 'w')
        fd.write(s)
        fd.close()
        os.system('chmod +x {}'.format(fp))

        return


    def parse_marker(self,line_m,):

        l_markers = line_m.split()
        pos_ref = int(l_markers[1])
        A_ref = l_markers[2]
        B_ref = l_markers[3]

        return pos_ref, A_ref, B_ref


    def BEAGLE3_divide(self,):

        d_indexes = {}
        for chrom in self.chroms:
            d_indexes[chrom] = {}

        ## parse the input tranches file describing where to cut the data
        fp_tranches = 'out_VariantRecalibrator/SNP.tranches'
        minVQSLOD = self.parse_minVQSLOD(fp_tranches,self.ts_filter_level,)

        ## open the recal file
        fp_recal = 'out_VariantRecalibrator/SNP.recal'

        for vcf in glob.glob('out_UnifiedGenotyper/*.vcf.gz'):
            assert os.path.getmtime(vcf) < os.path.getmtime(fp_recal)

        with open(fp_recal, 'r') as fd_recal:

            ## loop over the raw input variants to be recalibrated
            if os.path.islink('out_UnifiedGenotyper'):
                path = os.readlink('out_UnifiedGenotyper')
            else:
                path = 'out_UnifiedGenotyper'

            for chrom in self.chroms:
                print('chrom {}'.format(chrom), flush=True)
                s = '{}/{}.*.vcf.gz'.format(path,chrom,)
                l_files = glob.glob(s)
                if len(l_files) == 0:
                    print('no files found')
                    print(s)
                    sys.exit()
                l_files_sorted = self.sort_nicely(l_files)
                d_indexes[chrom] = self.loop_UG_out_BEAGLE3(
                    chrom,l_files_sorted,fd_recal,minVQSLOD)

        return d_indexes


    def parse_minVQSLOD(self,fp_tranches,ts_filter_level,):

        with open(fp_tranches) as fd:
            for line in fd:
                if line[0] == '#':
                    continue
                l = line.split(',')
                if l[0] == 'targetTruthSensitivity':
                    index = l.index('minVQSLod')
                    continue
                print(l[0])
                targetTruthSensitivity = float(l[0])
                if targetTruthSensitivity == ts_filter_level:
                    minVQSLOD = float(l[index])
                    break

        return minVQSLOD


    def alphanum_key(self,s):
        ## http://dave.st.germa.in/blog/2007/12/11/exception-handling-slow/
        NUM_RE = re.compile('([0-9]+)')
        return [int(c) if c.isdigit() else c for c in NUM_RE.split(s)]


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


    def vcf2beagle3_header(self,file):

        with gzip.open(file, 'rt') as fd_vcf:
            for line_vcf in fd_vcf:
                ## skip header
                if line_vcf[0] != '#':
                    break
                line_vcf_header = line_vcf
        l_samples = line_vcf_header.strip().split('\t')[9:]
        l = ['{} {} {}'.format(sample,sample,sample) for sample in l_samples]
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


    def loop_UG_out_BEAGLE3(self,chrom,l_files,fd_recal,minVQSLOD):

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
        header = self.vcf2beagle3_header(file)
        d_index2pos = {}

        ##
        ## size and edge
        ##
        size = self.i_BEAGLE_size*1000000
        edge = self.i_BEAGLE_edge*1000

        ## genotype likelihood file
        fp_out_prefix = 'in_BEAGLE/{}/{}'.format(chrom,chrom,)
        fp_markers = self.chrom2file(self.BEAGLE_markers, chrom)
        fp_phased = self.chrom2file(self.BEAGLE_phased, chrom)

        ##
        ## open files
        ##
        fd_phased = open(fp_phased, 'r')
        fd_markers = open(fp_markers, 'r')
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

        if os.path.isfile('{}.phased'.format(fp_out_prefix,)):
            os.remove('{}.phased'.format(fp_out_prefix,))
        fd = open('{}.phased'.format(fp_out_prefix,), 'w')
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
                        lines_out_markers2 += ['{}:{} {} {} {}\n'.format(
                            chrom,position,position,alleleA,alleleB,)]
                if bool_EOF1 == False:
                    lines_out1 += [line_beagle]
                    ## match
                    if bool_append_markphas == True:
                        lines_out_markers1 += [line_m]
                        lines_out_phased1 += [line_p]
                    ## mismatch
                    else: ## ms23/dg11 2013mar12
                        lines_out_markers1 += ['{}:{} {} {} {}\n'.format(
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
                        cmd = 'tail -n1 {}.{}.like | cut -d " " -f1'.format(
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
                            ) = self.remove_duplicate_lines_BEAGLE3(
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
                    print('%2s, pos_curr %9i, pos_init %9i, pos_term %9i, i %3i, n %5i'.format(
                        chrom,position,pos_init1,pos_term1,index,len(lines_out1)))
                    sys.stdout.flush()
                    if bool_append_to_next == False:
                        ## write/append current lines to file
                        fd_out = open('{}.{}.like'.format(fp_out_prefix,index,),mode)
                        fd_out.writelines(lines_out1)
                        fd_out.close()
                        fd = open('{}.{}.markers'.format(fp_out_prefix,index,),mode)
                        fd.writelines(lines_out_markers1)
                        fd.close()
    ##                    fd = open('{}.%i.phased'.format(fp_out_prefix,index,),mode)
                        fd = open('{}.phased'.format(fp_out_prefix,), 'a')
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
                            lines_out_markers1 += ['{}:{} {} {} {}\n'.format(
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
                ) = self.remove_duplicate_lines_BEAGLE3(
                    index, fp_out_prefix,
                    lines_out1, lines_out_markers1)
            mode = 'a'
            index -= 1
            pos_init1 = d_index2pos[index][0]
        else:
            mode = 'w'
        print('%2s, pos_curr %9i, pos_init %9i, pos_term %9i, i %3i, n %5i'.format(
            chrom,position,pos_init1,pos_term1,index,len(lines_out1)))

        ## write remaining few lines to output files
        ## without checking for duplicate positions
        lines_out_markers1 += fd_markers.readlines()
        lines_out_phased1 += fd_phased.readlines()

        ##
        ## write/append lines
        ##
        fd_out = open('{}.{}.like'.format(fp_out_prefix,index,),mode)
        fd_out.writelines(lines_out1)
        fd_out.close()

        fd = open('{}.{}.markers'.format(fp_out_prefix,index,),mode)
        fd.writelines(lines_out_markers1)
        fd.close()

##        fd = open('{}.%i.phased'.format(fp_out_prefix,index,),mode)
        fd = open('{}.phased'.format(fp_out_prefix,), 'a')
        fd.writelines(lines_out_phased1)
        fd.close()

        ## append position to dictionary
        d_index2pos[index] = [pos_init1,pos_term1,]

        ##
        ## close all files after looping over last line
        ##
        fd_markers.close()
        fd_phased.close()

        with open('summaryVR.txt', 'a') as file_summary:
            file_summary.write('{} {}\n'.format(chrom,cnt_variants))
##        print('{} %i\n'.format(chrom,cnt_variants),file='summaryVR.txt')

        ##
        ## file i/o checks
        ##
##        self.BEAGLE3_divide_fileIO_checks(chrom,fp_phased,)

        return d_index2pos


    def vcf2beagle(self,l_vcf):

        ## this function is very slow; especially split, sum and pow

        ## append chrom:pos alleleA alleleB
        line_beagle = '{}:{} {} {}'.format(l_vcf[0],l_vcf[1],l_vcf[3],l_vcf[4],)

        index = l_vcf[8].split(':').index('PL')
        for s_vcf in l_vcf[9:]:
            ## variant not called
            if s_vcf[:3] == './.':
                line_beagle += ' 0.3333 0.3333 0.3333'
                continue
            l_probs = []
            l_log10likelihoods = s_vcf.split(':')[index].split(', ')
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
                line_beagle += ' %6.4f'.format(prob/sum(l_probs))
            if ', ' in l_vcf[4]:
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
        if not l_vcf[3] in ['A', 'C', 'G', 'T',]:
            return True,None
        ## dialleic SNP
        if l_vcf[4] in ['A', 'C', 'G', 'T',]:
            return False,True
        ## triallelic and other non-diallelic SNPs
##        elif l_vcf[4] in [
##            'A,C', 'A,G', 'A,T', 'C,G', 'C,T', 'G,T',
##            'A,C,G', 'A,C,T', 'A,G,T', 'C,G,T',
##            ]:
##            return False,False,False
##        ## UG DISCOVERY
##        elif l_vcf[4] in iter(', '.join(tup) for tup in itertools.combinations('ACGT',2)):
        ## UG GENOTYPE_GIVEN_ALLELES
        elif l_vcf[4] in iter(', '.join(tup) for tup in itertools.permutations('ACGT',2)):
            return False,False
##        ## UG DISCOVERY
##        elif l_vcf[4] in iter(', '.join(tup) for tup in itertools.combinations('ACGT',3)):
        ## UG GENOTYPE_GIVEN_ALLELES
        elif l_vcf[4] in iter(', '.join(tup) for tup in itertools.permutations('ACGT',3)):
            return False,False
        ## skip insertions
        else:
            return True,None

        return


    def BEAGLE3_divide_fileIO_checks(self,chrom,fp_phased,):

        '''Compare VR out and BEAGLE in'''

        cmd = 'cat {}'.format(fp_phased)
        cmd += "| awk 'NR>1{print $2}' | sort -u > panel0in{}" %(chrom)
        self.execmd(cmd)

        cmd = "awk 'FNR>1{print $2}' in_BEAGLE/{}/{}.phased" %(chrom,chrom)
        cmd += "| sort -u > panel0out{}" %(chrom)
        self.execmd(cmd)

        fp_tranches = 'out_VariantRecalibrator/SNP.tranches'
        minVQSLOD = self.parse_minVQSLOD(fp_tranches,self.ts_filter_level,)
        fp_recal = 'out_VariantRecalibrator/SNP.recal'

##        with contextlib.ExitStack() as stack:
##            fd_recal = stack.enter_context(open(fp_recal))
##            fd_out = stack.enter_context(open('panel2in{}'.format(chrom), 'w'))
##            fd_UG = stack.enter_context(fileinput.fileinput(files_sorted))
##            l_files = self.sort_nicely(
##                glob.glob('out_UnifiedGenotyper/{}.*.vcf.gz'.format(chrom)))
##            assert len(l_files) > 0
##                l_files_sorted = self.sort_nicely(l_files)
##            fd_vcf = fileinput.input(
##                files=l_files,openhook=self.hook_compressed_text)
##            for l_vcf in self.generate_line_vcf_PASS_split(
##                fd_vcf,fd_recal,minVQSLOD,))

        with open(fp_recal, 'r') as fd_recal, open('panel2in{}'.format(chrom), 'w') as fd_out:
            for line in fd_recal:
                if line[0] == '#':
                    continue
                VQSLOD, l_recal = self.parse_VQSLOD(line)
                if VQSLOD < minVQSLOD:
                    continue
                ## slow conversion to integer just to make sure it is an integer
                pos = int(l_recal[1])
                fd_out.write('{}:{}\n'.format(chrom,pos))
        self.execmd('sort panel2in{} -o panel2in{}'.format(chrom,chrom))

        cmd = "awk 'FNR>1{print $1}' in_BEAGLE/{}/{}.*.like" %(chrom,chrom)
        cmd += ' | sort -u > panel2out{}'.format(chrom)
        self.execmd(cmd)

        cmd = "awk '{print $1}' in_BEAGLE/{}/{}.*.markers" %(chrom,chrom)
        cmd += ' | sort -u > markers{}'.format(chrom)
        self.execmd(cmd)

        ## check that like (panel2out) is identical to PBI (panel2in)
        cmd = 'comm -3 panel2out{} panel2in{}'.format(chrom,chrom)
        i = int(os.popen('{} | wc -l'.format(cmd)).read())
        if i > 0:
            print(os.popen('{} | head'.format(cmd)).readlines())
            print(cmd)
            sys.exit()

        ## check that panel0out is a subset of panel0in
        cmd = 'comm -23 panel0out{} panel0in{}'.format(chrom,chrom)
        i = int(os.popen('{} | wc -l'.format(cmd)).read())
        if i > 0:
            print(os.popen('{} | head'.format(cmd)).readlines())
            print(cmd)
            sys.exit()

        ## check that like (panel2out) is a true subset of markers
        cmd = 'comm -32 panel2out{} markers{}'.format(chrom,chrom)
        i = int(os.popen('{} | wc -l'.format(cmd)).read())
        if i > 0:
            print(os.popen('{} | head'.format(cmd)).readlines())
            print(cmd)
            sys.exit()

        ## check that phased (panel0out) is a true subset of markers
        cmd = 'comm -32 panel0out{} markers{}'.format(chrom,chrom)
        i = int(os.popen('{} | wc -l'.format(cmd)).read())
        if i > 0:
            print(os.popen('{} | head'.format(cmd)).readlines())
            print(cmd)
            sys.exit()

        ## check that markers is identical to panel2in+panel0in
        cmd = 'cat panel2out{} panel0out{} | sort -u > panelsout'.format(chrom,chrom)
        self.execmd(cmd)
        cmd = 'comm -3 panelsout markers{}'.format(chrom)
        i = int(os.popen('{} | wc -l'.format(cmd)).read())
        if i > 0:
            print(os.popen(cmd).readlines()[:10])
            print(cmd)
            sys.exit()

        for affix in ['panel2out', 'panel2in', 'panel0in', 'panel0out', 'markers',]:
            os.remove('{}{}'.format(affix,chrom,))
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
        cmd = 'tail -n1 {}.{}.like | cut -d " " -f1'.format(
            fp_out_prefix,index-1)
        pos_prev_like = int(os.popen(cmd).read().split(':')[1])

        cmd = 'tail -n1 {}.{}.markers'.format(
            fp_out_prefix,index-1)
        l = os.popen(cmd).read().strip().split()
        pos_prev_markers = int(l[1])
        alleleA_prev = l[2]
        alleleB_prev = l[3]

##        cmd = 'tail -n1 {}.%i.phased | cut -d " " -f2'.format(
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


    def BEAGLE3(self):

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
        fp_in_recal = 'out_VariantRecalibrator/{}.recal'.format(mode)
        fp_in_tranches = 'out_VariantRecalibrator/{}.tranches'.format(mode)
        if self.check_in(
            'VariantRecalibrator',[fp_in_recal,fp_in_tranches,],
            'touch/VariantRecalibrator.SNP.touch'):
            sys.exit()
        stoptmp

        ##
        ## 2) touch
        ##
        if self.touch('BEAGLE'):
            return

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
            d_indexes = self.BEAGLE3_divide()
            s = ''
            for chrom,d_index2pos in d_indexes.items():
                for index,[pos1,pos2,] in d_index2pos.items():
                    s += '{} {} {} {}\n'.format(chrom,index,pos1,pos2,)
            fd = open('BEAGLE_divide_indexes.txt', 'w')
            fd.write(s)
            fd.close()
        else:
            print('BEAGLE_divide_indexes.txt exists')
            d_indexes = {}
            for chrom in self.chroms:
                d_indexes[chrom] = {}
            fd = open('BEAGLE_divide_indexes.txt', 'r')
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

            print('bsub BEAGLE {}'.format(chrom))

##            J = '{}{}[%i-%i]'.format('BEAGLE',chrom,1,max(d_indexes[chrom].keys()),)
            for index in d_indexes[chrom].keys():

                ## finished?
                fn_out = 'out_BEAGLE/{}/{}.{}.like.gprobs.gz'.format(
                    chrom,chrom,index)
                if os.path.isfile(fn_out):
                    continue

                ## started?
                fn = 'LSF/BEAGLE/{}.{}.out'.format(chrom,index)
                if os.path.isfile(fn):
                    if os.path.getsize(fn):
                        with open(fn) as f:
                            if 'iteration' in f.readlines()[-1]:
                                continue

                print(chrom,index)

                J = '{}.{}[{}-{}]'.format('BEAGLE',chrom,index,index,)
                LSF_affix = '{}/{}.%%I'.format('BEAGLE',chrom)
                cmd = self.bsub_cmd(
                    'BEAGLE', J, memMB=memMB, LSF_affix=LSF_affix,
                    chrom=chrom, queue=queue,)
                os.system(cmd)

        return


    def init_java(self, jar, memMB, java='java', bool_checkpoint=False):

        s = '{} -Djava.io.tmpdir={}'.format(java, 'tmp')
        ## set maximum heap size
        s += ' -Xmx{}m'.format(memMB)
        if bool_checkpoint:
            s += ' -XX:-UsePerfData -Xrs '
        s += ' \\\n -jar {}'.format(jar)

        return s


    def BEAGLE3_write_shell_script(self,memMB,):

        fp_out = 'out_BEAGLE/$CHROMOSOME/$CHROMOSOME.${LSB_JOBINDEX}.like'

        ## initiate shell script
        lines = ['#!/bin/bash\n']

        ## parse chromosome from command line
        lines += ['CHROMOSOME=$1\n']

        lines += ['if [ -s {}.gprobs.gz ]; then exit; fi\n'.format(fp_out)] ## redundant

##        ## init cmd
##        lines += ['cmd="']

        ##
        ## initiate BEAGLE
        ##
        s_java = self.init_java(self.fp_software_beagle,memMB)
        lines += ['{} \\'.format(s_java)]

        lines += self.body_BEAGLE(fp_out,)

        ## term cmd
        lines += self.term_cmd('BEAGLE',['{}.gprobs.gz'.format(fp_out)],)

        ## write shell script
        self.write_shell('shell/BEAGLE.sh', lines)

        return


    def body_BEAGLE3(self,fp_out,):

        fp_like = 'in_BEAGLE/$CHROMOSOME/$CHROMOSOME.${LSB_JOBINDEX}.like'
##        fp_phased = 'in_BEAGLE/$CHROMOSOME/$CHROMOSOME.${LSB_JOBINDEX}.phased'
        fp_phased = 'in_BEAGLE/$CHROMOSOME/$CHROMOSOME.phased'
        fp_markers = 'in_BEAGLE/$CHROMOSOME/$CHROMOSOME.${LSB_JOBINDEX}.markers'

        lines = []

##like=<unphased likelihood data file> where <unphased likelihood data file> is
##the name of a genotype likelihoods file for unphased, unrelated data
##(see Section 2.2). You may use multiple like arguments if data from different
##cohorts are in different files.
        lines += [' like={} \\'.format(fp_like)]
####arguments for phasing and imputing data ...
#### Arguments for specifying files
## phased=<phased unrelated file> where <phased unrelated file> is the name of a
## Beagle file containing phased unrelated data (see Section 2.1).
## You may use multiple phased arguments if data from different cohorts are in
## different files.
##        lines += [' phased=in_BEAGLE/ALL.chr$CHROMOSOME.phase1_release_v3.20101123.filt.renamed.bgl \\']
        lines += [' phased={} \\'.format(fp_phased)]
####  unphased=<unphased data file>                     (optional)
####  phased=<phased data file>                         (optional)
####  pairs=<unphased pair data file>                   (optional)
####  trios=<unphased trio data file>                   (optional)
####  like=<unphased likelihood data file>              (optional)
##markers=<markers file> where <markers file> is the name of the markers file containing
## marker identifiers, positions, and alleles described in Section 2.4.
## The markers argument is optional if you specify only one Beagle file,
## and is required if you specify more than one Beagle file.
        lines += [' markers={} \\'.format(fp_markers)]
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
        lines += [' nsamples={} \\'.format(int(self.i_BEAGLE_nsamples))]

        lines += [' niterations=10 \\'] ## default 10

        lines += [' omitprefix=true \\'] ## default false

        lines += [' verbose=false \\'] ## default false

##        lines += [' lowmem=false \\'] ## default false
        lines += [' lowmem=true \\'] ## default false

        ## non-optional output prefix
        lines += [' out={} \\'.format(fp_out)]

        return lines


    def check_in(self, analysis_type, l_fp_in, fp_touch,):

        d_l_fp_out = {}

        with open(fp_touch) as fd:
            s = fd.read()
        l_fp_out = s.split('\n')
        d_l_fp_out['touch'] = l_fp_out

        ## todo: use os.walk here instead...
        for dirname in ['',]:
            d_l_fp_out[dirname] = []
            l = os.listdir(os.path.join(dirname, 'out_{}'.format(analysis_type)))
            for s in l:
                path1 = os.path.join('out_{}'.format(analysis_type),s)
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
                print('{} and possibly {} other files not generated.'.format(
                    list(set(l_fp_in)-set(l_fp_out))[0],
                    len(set(l_fp_in)-set(l_fp_out))-1,))
                print('dirname', dirname)
                print('{} has not run to completion. Exiting.'.format(analysis_type))
                bool_exit = True
#                print(inspect.stack()[1])
##                sys.exit()

        return bool_exit


    def touch(self,analysis_type):

        bool_return = False
        fn_touch = 'touch/{}.touch'.format(analysis_type)
        if os.path.isfile(fn_touch):
            if self.verbose == True:
                print('in progress or completed:', analysis_type)
            bool_return = True
        else:
            if not os.path.isdir(os.path.dirname(fn_touch)):
                os.mkdir(os.path.dirname(fn_touch))
            self.execmd('touch {}'.format(fn_touch))

        return bool_return


    def write_brestart(self,):

        with open('brestart.sh', 'w') as f:
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
        lines += ['CHROM=$1']
        lines += ['BAM=$2']
        lines += ['BAMBASENAME=$(basename $BAM | rev | cut -d "." -f2- | rev)']
        lines += ['out=out_HaplotypeCaller/$CHROM/$BAMBASENAME.vcf.gz']
##        ## exit if job started
##        lines += ['if [ -s $out ]; then exit; fi\n']
        ## exit if job finished
        lines += ['if [ -s $out.tbi ]; then exit; fi\n']

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
        self.write_shell('shell/{}.sh'.format(analysis_type),lines,)

        return


    def body_HaplotypeCaller(self,):

        '''Write walker specific command line arguments.'''

        ## https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_haplotypecaller_HaplotypeCaller.html

        lines = []

        ##
        ## Inherited arguments
        ##

        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_CommandLineGATK.html#--input_file
        lines += [' --input_file $BAM \\']

        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_CommandLineGATK.html#--intervals
        if self.intervals:
            lines += ['--intervals {} \\'.format(self.intervals)]
            ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_CommandLineGATK.html#--interval_set_rule
            lines += ['--interval_set_rule INTERSECTION \\']
            pass
        else:
            lines += ['--intervals $CHROM \\']

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
            lines += [' --alleles {} \\'.format(self.alleles)]

        ## dbSNP file. rsIDs from this file are used to populate the ID column of the output.
        ## https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_haplotypecaller_HaplotypeCaller.html#--dbsnp
        if self.dbsnp:
            lines += [' --dbsnp {} \\'.format(self.dbsnp)]

        ##
        ## Optional Outputs
        ##

        ## https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_haplotypecaller_HaplotypeCaller.html#--out
        lines += [' --out $out \\']

        ##
        ## Optional Parameters
        ##

        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_haplotypecaller_HaplotypeCaller.html#--genotyping_mode
        lines += [' -gt_mode {} \\'.format(self.genotyping_mode)] ## default value DISCOVERY

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
        lines += [' {} \\'.format(s_annotation)]

        lines += [' --emitRefConfidence GVCF \\']

        lines += ['\n']

        return lines


    def ApplyRecalibration_emit_INDELs(self):

        analysis_type = T = 'ApplyRecalibration'
        mode = 'SNP'
##        num_threads = 1
        queue = 'normal'
        memMB = 15900

        for mode in ['SNP']:

            ## check input existence
            fp_recal = 'out_VariantRecalibrator/{}.recal'.format(mode)
            fp_tranches = 'out_VariantRecalibrator/{}.tranches'.format(mode)
            if self.check_in(
                'VariantRecalibrator', [fp_recal, fp_tranches,],
                'touch/VariantRecalibrator.SNP.touch'):
                continue

            ## check if process already started and otherwise lock for this process
            if self.touch('{}.{}'.format(analysis_type, mode)):
                continue

            ## write shell script
            lines = []
            lines += ['chrom=$1']
            lines += ['mode=$2']
            lines += ['out=out_{}/$mode.$chrom.vcf.gz\n'.format(T)]
            lines += self.init_GATK_cmd(analysis_type, memMB,)
            lines += [
                ' --input out_GenotypeGVCFs/$chrom.vcf.gz \\']
            lines += [' --recal_file {} \\'.format(fp_recal)]
            lines += [' --tranches_file {} \\'.format(fp_tranches)]
            lines += [' --out $out \\']
            lines += [' --mode $mode \\']
            lines += [' --excludeFiltered \\']
            lines += [
                ' --ts_filter_level {:.1f} \\'.format(self.ts_filter_level)]
            lines += self.term_cmd(
                '{}.$mode'.format(analysis_type), ['$out'],)
            if not os.path.isfile('shell/{}.sh'.format(T)):
                self.write_shell('shell/{}.sh'.format(T), lines,)

            ## create LSF folder
            if not os.path.isdir('LSF/{}'.format(T)):
                os.mkdir('LSF/{}'.format(T))

            ## execute shell script
            for chrom in self.chroms:
                J = 'AR.{}.{}'.format(mode, chrom)
                cmd = self.bsub_cmd(
                    analysis_type, J, memMB=memMB,
                    LSF_affix='{}/{}.{}'.format(analysis_type, mode, chrom),
                    queue=queue, mode=mode, chrom=chrom,)
                self.execmd(cmd)

        return



    def init_GATK_cmd(self,analysis_type,memMB,bool_checkpoint=False):

        lines = []

        ## exit if output exists
        lines += ['if [ -s $out ]; then exit; fi']

        ## create output folder
        lines += ['mkdir -p $(dirname $out)']

        s = ''
        ## Java version alert: starting with release 2.6, GATK now requires Java 1.7. See Version Highlights for 2.6 for details.
        ## http://www.broadinstitute.org/gatk/guide/article?id=2846
        s_java = self.init_java(
            self.fp_GATK, memMB, java=self.java, bool_checkpoint=bool_checkpoint)
        s += ' {} \\'.format(s_java)
        lines += ['\n{}'.format(s)]

        ## CommandLineGATK, required, in
        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_CommandLineGATK.html#--analysis_type
        lines += [' --analysis_type {} \\'.format(analysis_type)]
        ## CommandLineGATK, optional, in
        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_CommandLineGATK.html#--reference_sequence
        lines += [' --reference_sequence {} \\'.format(self.reference_sequence)]
        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_CommandLineGATK.html#-nct
        if self.nct != 1:
            lines += [' --num_cpu_threads_per_data_thread {} \\'.format(self.nct)]


        return lines


    def bsub_cmd(
        self,
        analysis_type,
        J,
        queue='normal',memMB=4000,
        LSF_affix=None,
        chrom=None, index=None, pos1=None, pos2=None,
        bool_checkpoint=False,
        num_threads=None, bam=None, mode=None,
        ):

        if not LSF_affix:
            LSF_affix = '{}/{}'.format(analysis_type,analysis_type,)

        cmd = 'bsub -J"{}" -q {}'.format(J,queue,)
        cmd += ' -G {}'.format(self.project)
        cmd += " -M%i -R'select[mem>%i] rusage[mem=%i]'" %(
            memMB,memMB,memMB,)
        cmd += ' -o {}/LSF/{}.out'.format(os.getcwd(), LSF_affix)
        cmd += ' -e {}/LSF/{}.err'.format(os.getcwd(), LSF_affix)
        if num_threads:
            cmd += ' -n{} -R"span[hosts=1]"'.format(num_threads)
        if bool_checkpoint:
            cmd += ' -k "{} method=blcr 710"'.format(
                os.path.join(os.getcwd(), 'checkpoint'))
            cmd += ' -r'
        if bool_checkpoint:
            cmd += ' cr_run'
        cmd += ' bash {}/shell/{}.sh'.format(os.getcwd(),analysis_type,)
        for x in (chrom, index, bam, mode, pos1, pos2):
            if x:
                cmd += ' {}'.format(x)

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
            lines += ['if [ ! -s {} ]; then exit; fi'.format(fp_out)]
            lines += ['echo {} >> touch/{}.touch'.format(fp_out, analysis_type,)]

        if extra:
            lines += ['{}\n'.format(extra)]

        lines += ['bash ./rerun.sh']
        lines += ['fi']

        ## write continuation shell script
        ## do not continue as part of previous command
        ## as this will influence CPU statistics
        ## and risk job of hitting CPU walltime
        s = "bsub -R 'select[mem>1500] rusage[mem=1500]' -M1500 \\\n"
        s += ' -o LSF/rerun.out \\\n'
        s += ' -e LSF/rerun.err \\\n'
        s += ' -G {} \\\n'.format(self.project)
        s += ' bash ./rerun_python.sh'
        fd = open('rerun.sh', 'w')
        fd.write(s)
        fd.close()
        self.execmd('chmod +x rerun.sh')

        s = self.args_to_command_line()
        fd = open('rerun_python.sh', 'w')
        fd.write(s)
        fd.close()
        self.execmd('chmod +x rerun_python.sh')

        return lines


    def args_to_command_line(self):

        s = ''
        s += ' {}'.format(sys.executable)
        s += ' {}'.format(sys.argv[0])
        for k, v in vars(self.args).items():
            if v == False:
                continue
            elif v == None:
                continue
            elif v == True and type(v) == bool:
                v = ''
            else:
                pass
            if type(v) == list:
                v = ' '.join(v)
            s += ' --{} {}'.format(k, str(v))

        return s


    def is_file(self, str_):
        if not os.path.isfile(str_) and not os.path.islink(str_):
            msg = '{} is neither a readable file nor a symbolic link' % str_
            raise argparse.ArgumentTypeError(msg)
        return str_


    def is_file_or_dir(self, str_):
        print(str_)
        if not any([
            os.path.isfile(str_),os.path.islink(str_),os.path.isdir(str_)]):
            msg = '{} is neither a readable file nor a directory' % str_
            raise argparse.ArgumentTypeError(msg)
        return str_


    def add_arguments(self,parser):

        ## required arguments

        parser.add_argument(
            '--fp_bams', '--bam', '--bams', '--input',
            help='Path to BAM and/or directory containing BAMs',
            nargs='+', required=True, type=self.is_file_or_dir)

        parser.add_argument('--coverage', required=True, type=float)

        parser.add_argument(
            '--fp_GATK', '--GATK', '--gatk', '--jar', required=True,
            help='File path to GATK',)

        parser.add_argument('--project', required=True)

        parser.add_argument('--arguments', '--args')

        parser.add_argument('--java', required=True)

        ##
        ## CommandLineGATK arguments
        ##

        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_CommandLineGATK.html#--intervals
        parser.add_argument(
            '--intervals', '-L',
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
            '--alleles', help='The set of alleles at which to genotype when --genotyping_mode is GENOTYPE_GIVEN_ALLELES.',)

        ## https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_haplotypecaller_HaplotypeCaller.html#--dbsnp
        parser.add_argument('--dbsnp', '-D')

        ## Optional Parameters

        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_genotyper_UnifiedGenotyper.html#--genotyping_mode
        parser.add_argument(
            '--genotyping_mode', '-gt_mode',
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

        ## What VQSR training sets / arguments should I use for my specific project?
        ## https://www.broadinstitute.org/gatk/guide/article?id=1259
        parser.add_argument(
            '--ts_filter_level_SNP', '--ts_SNP', type=float, required=True,)

        parser.add_argument(
            '--ts_filter_level_INDEL', '--ts_INDEL', type=float, required=True,)

        parser.add_argument('--AR_input')

        ##
        ## BEAGLE
        ##
        parser.add_argument(
            '--fp_software_beagle', '--beagle', '--BEAGLE', '--BEAGLEjar',
            help='File path to BEAGLE .jar file (e.g. beagle_3.3.2.jar)',
            required=True,
            )

        parser.add_argument(
            '--i_BEAGLE_size',  # BEAGLE3 argument
            help='Size (Mbp) of divided parts.',
            type=int,default=2,  # CPU bound (2Mbp=22hrs,3090MB) with lowmem option...
            )

        ## BEAGLE 3
        parser.add_argument(
            '--i_BEAGLE_edge',
            help='Window size (kbp) at either side of the divided part to avoid edge effects.',
            type=int,default=150,
            )

        ## BEAGLE 3
        s_help = 'Phased file to be divided (e.g.'
        s_help += ' ALL.chr$CHROM.phase1_release_v3.20101123.filt.renamed.bgl'
        parser.add_argument(
            '--BEAGLE_phased', '--phased',
            help=s_help,
            required=False,
            )

        ## BEAGLE 3
        s_help = 'Markers file to be divided (e.g.'
        s_help += ' ALL.chr$CHROM.phase1_release_v3.20101123.filt.renamed.markers'
        parser.add_argument(
            '--BEAGLE_markers', '--markers',
            help=s_help, required=False,
            )

        parser.add_argument('--i_BEAGLE_nsamples', default=20)

        ##
        ## optional arguments
        ##
        parser.add_argument(
            '--checkpoint', dest='bool_checkpoint', action='store_true', default=False)

        parser.add_argument(
            '--chroms', type=str, nargs='+',
            default=[str(i+1) for i in range(22)]+['X', 'Y', 'MT',])

        parser.add_argument(
            '--ped', type=self.is_file)

        return parser


    def parse_arguments(self):

        parser = argparse.ArgumentParser()

        parser = self.add_arguments(parser)

        ## parse arguments to argparse NameSpace
        self.args = namespace_args = parser.parse_args()

        ## setatrr
        for k,v in vars(namespace_args).items():
            setattr(self,k,v)

        if self.fp_GATK is None and self.fp_options is None:
            parser.error('--GATK or --arguments')

        s_arguments = ''
        for k,v in vars(namespace_args).items():
            s_arguments += '{} {}\n'.format(k,v)

        if self.arguments == None or self.arguments == 'None':
            self.arguments = '{}.arguments'.format(self.project)
            fd = open(self.arguments, 'w')
            fd.write(s_arguments)
            fd.close()
        else:
            fd = open(self.arguments, 'r')
            lines = fd.readlines()
            fd.close()
            for line in lines:
                l = line.strip().split()
                k = l[0]
                v = l[1]
                setattr(self, k, v)

        self.bams = []
        for fp in self.fp_bams:
            if os.path.isdir(fp):
                self.bams += glob.glob(os.path.join(fp, '*.bam'))
            elif os.path.isfile(fp):
                self.bams += [fp]
            else:
                stop_take_care_of_symlinks

        return


    def __init__(self):

        ## parse command line arguments
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
