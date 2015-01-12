#!/bin/python3

## Tommy Carstensen, Wellcome Trust Sanger Institute, 2012-2015

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
import urllib


## README

## http://www.broadinstitute.org/gatk/about#typical-workflows
## http://www.broadinstitute.org/gatk/guide/topic?name=best-practices

class main():

    def main(self):

        if self.args.AR_input:
            self.run_ApplyRecalibration()

        if self.coverage >= 15:
            self.HaplotypeCaller()
            self.CombineGVCFs()
            self.GenotypeGVCFs()
        elif self.coverage < 15:
            self.UnifiedGenotyper()

##        ## 1000G_phase1.snps.high_confidence.b37.vcf.gz only contains chromosomes 1-22 and X
##        self.chroms.remove('X')  ## tmp!!!
##        self.chroms.remove('Y')  ## tmp!!!
##        self.chroms.remove('MT')  ## tmp!!!
        self.VariantRecalibrator()

        self.ApplyRecalibration()

        ## phased imputation reference panels not available for chrY
        try:
            self.chroms.remove('Y')
        except ValueError:
            pass
##        ## Exception in thread "main" java.lang.ArrayIndexOutOfBoundsException: -2
##        self.chroms.remove('X')
        self.BEAGLE4()

        return


    def UnifiedGenotyper(self):

        T = analysis_type = 'UnifiedGenotyper'
        queue = 'normal'
        size_fragment_bp = 5*10**5  # normal
        queue = 'long'
        size_fragment_bp = 2*10**6  # long
        ## http://gatkforums.broadinstitute.org/discussion/1975/how-can-i-use-parallelism-to-make-gatk-tools-run-faster
        nct = 3
        nt = 8
        ## Each data thread needs to be given the full amount of memory you’d normally give a single run. So if you’re running a tool that normally requires 2 Gb of memory to run, if you use -nt 4, the multithreaded run will use 8 Gb of memory. In contrast, CPU threads will share the memory allocated to their “mother” data thread, so you don’t need to worry about allocating memory based on the number of CPU threads you use.
        memMB = 2900+8*len(self.bams)
        memMB = nt*4000-100
        memMB = 47900

        ## Parse chromosome ranges.
        d_chrom_ranges = self.parse_chrom_ranges()

        ## Write bam lists for
        ## all (autosomes, PAR), male (Y) and female (nonPAR X) samples.
        self.write_bam_and_XL_lists(d_chrom_ranges)

        ## Create folders.
        if not os.path.isdir('tmp'):
            os.mkdir('tmp')

        ## touch/lock
        if self.touch(analysis_type):
            return

        ## write shell script
        self.shell_UG(T, memMB, nct, nt)

        ## execute shell script
        ## loop over chroms
        for chrom in self.chroms:

            d_sex2bamlist, XL = self.get_sex_and_XL(chrom)
            cstart = d_chrom_ranges[chrom][0]
            cstop = d_chrom_ranges[chrom][1]

            for sex in d_sex2bamlist.keys():

                if not os.path.isfile(d_sex2bamlist[sex]):
                    continue
                if not os.path.getsize(d_sex2bamlist[sex]):
                    continue

                sample_ploidy = self.get_ploidy(chrom, sex)

                for i in range(
                    1, 1+math.ceil((cstop-cstart)/size_fragment_bp)):

                    pos1 = 1+cstart//size_fragment_bp+(i-1)*size_fragment_bp
                    pos2 = min(cstop, (pos1-1)+size_fragment_bp)

                    affix = '{}/{}/{}'.format(T, chrom, i)
                    if sex and chrom == 'X':
                        affix += '.{}'.format(sex)

                    ## Skip if output was generated.
                    if os.path.isfile(
                        'out_{}.vcf.gz.tbi'.format(affix)):
                        continue

                    ## Skip if output is being generated.
                    if os.path.isfile('LSF/{}.err'.format(affix)):
                        if time.time()-os.path.getmtime(
                            'LSF/{}.err'.format(affix)) < 300:
                            continue
                        if os.path.isfile('LSF/{}.out'.format(affix)):
                            os.remove('LSF/{}.out'.format(affix))
                        os.remove('LSF/{}.err'.format(affix))

                    os.makedirs(os.path.dirname(
                        'LSF/{}'.format(affix)), exist_ok=True)
                    os.makedirs(os.path.dirname(
                        'out_{}'.format(affix)), exist_ok=True)

                    variables = []
                    if chrom in ('PAR1', 'PAR2'):
                        variables += ['chrom=X']
                    else:
                        variables += ['chrom={}'.format(chrom)]
                    variables += ['pos1={}'.format(pos1)]
                    variables += ['pos2={}'.format(pos2)]
                    if XL:
                        variables += ['XL={}'.format(XL)]
                    variables += ['input_file={}'.format(d_sex2bamlist[sex])]
                    variables += ['out=out_{}.vcf.gz'.format(affix)]
                    variables += ['nct={}'.format(nct)]
                    variables += ['nt={}'.format(nt)]
                    variables += ['sample_ploidy={}'.format(sample_ploidy)]
                    variables = ' '.join(variables)

                    LSB_JOBNAME = '{} {} {} {}'.format('UG', chrom, i, sex)
                    cmd = self.bsub_cmd(
                        T, LSB_JOBNAME, LSF_affix=affix,
                        LSF_memMB=memMB, LSF_queue=queue, LSF_n=nct,
                        variables=variables)
                    self.execmd(cmd)

        return


    def HaplotypeCaller(self):

        T = analysis_type = 'HaplotypeCaller'
        queue = 'long'  # 4x split per sample
        memMB = 3900
#        queue = 'basement'  # >4x split per sample
        queue = 'normal'  # 4x split per chromosome
        nct = 4

        ## Create folders.
        os.makedirs('tmp', exist_ok=True)

        ## touch/lock
        if self.touch(analysis_type):
            return

        ## Parse chromosome ranges.
        d_chrom_ranges = self.parse_chrom_ranges()

        ## Write bam lists for
        ## all (autosomes, PAR), male (Y) and female (nonPAR X) samples.
        self.write_bam_and_XL_lists(d_chrom_ranges)

        ## write shell script
        self.shell_HC(analysis_type, memMB)

        ## execute shell script
        for chrom in self.chroms:

            d_sex2bamlist, XL = self.get_sex_and_XL(chrom)

            for sex in d_sex2bamlist:
                with open(d_sex2bamlist[sex], 'r') as f:
                    list_bams = [line.strip() for line in f]

                sample_ploidy = self.get_ploidy(chrom, sex)

                for bam in self.bams:
                    basename = os.path.splitext(os.path.basename(bam))[0]

                    if not bam in list_bams:
                        continue

                    affix = '{}/{}/{}'.format(T, chrom, i)
                    if sex:
                        affix += '.{}'.format(sex)

                    ## Skip bam and chromosome if output was generated.
                    if os.path.isfile(
                        'out_{}/{}/{}.vcf.gz.tbi'.format(T, chrom, basename)):
                        continue

                    ## Skip bam and chromosome if output is being generated.
                    if os.path.isfile('LSF/{}.err'.format(affix)):
                        ## file changed within the past 5 minutes?
                        if time.time()-os.path.getmtime(
                            'LSF/{}.err'.format(affix)) < 300:
                            continue
                        else:
                            os.remove('LSF/{}.err'.format(affix))
                    if os.path.isfile('LSF/{}.out'.format(affix)):
                        os.remove('LSF/{}.out'.format(affix))

                    os.makedirs(os.path.dirname(
                        'LSF/{}'.format(affix)), exist_ok=True)
                    os.makedirs(os.path.dirname(
                        'out_{}'.format(affix)), exist_ok=True)

                    variables = []
                    if chrom in ('PAR1', 'PAR2'):
                        variables += ['chrom=X']
                    else:
                        variables += ['chrom={}'.format(chrom)]
                    variables += ['input_file={}'.format(bam)]
                    if XL:
                        variables += ['XL={}'.format(XL)]
                    variables += ['out=out_{}.vcf.gz'.format(affix)]
                    variables += ['nct={}'.format(nct)]
                    variables += ['nt={}'.format(nt)]
                    variables += ['sample_ploidy={}'.format(sample_ploidy)]
                    variables = ' '.join(variables)

                    J = '{} {}'.format('HC', basename)
                    cmd = self.bsub_cmd(
                        T, LSB_JOBNAME, LSF_affix=affix,
                        memMB=memMB, queue=queue, num_threads=nct,
                        variables=variables)
                    self.execmd(cmd)

        return


    def shell_HC(self, T, memMB):

        lines = ['#!/bin/bash\n']
        lines += ['basename=$(basename $file_input | rev | cut -d "." -f2- | rev)']
##        ## exit if job started
##        lines += ['if [ -s $out ]; then exit; fi\n']
        ## exit if job finished
        lines += ['if [ -s $out.tbi ]; then exit; fi\n']

        ## initiate GATK command
        lines += self.init_GATK_cmd(T, memMB)

        ## append GATK command options
        lines += self.body_HaplotypeCaller()

        ## terminate shell script
##        lines += self.term_cmd(
##            analysis_type, ['$out.tbi'], extra='tabix -p vcf $out')
        lines += self.term_cmd(T, ['$out.tbi'])

        ## write shell script
        self.write_shell('shell/{}.sh'.format(T),lines,)

        return


    def parse_chrom_ranges(self):

        d_chrom_ranges = {}

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
            d_chrom_ranges[chrom] = (1, chrom_len)
##            ## break, when we reach the Y chromosome
##            if chrom == 'Y':
##                break

        if self.build == 37:
            url = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/GCF_000001405.25.regions.txt'
        elif self.build == 37:
            url = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/GCF_000001405.28.regions.txt'
        else:
            print('url not provided for build', self.build)
            sys.exit()
        with urllib.request.urlopen(url) as f:
            for line in f:
                line = line.decode('utf-8')
                if line[0] == '#':
                    continue
                l = line.split()
                if l[0] in ('PAR#1', 'PAR#2') and l[1] == 'X':
                    d_chrom_ranges[l[0].replace('#','')] = (int(l[2]),int(l[3]))

        return d_chrom_ranges


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

            os.makedirs('lists', exist_ok=True)
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

            os.makedirs('LSF/{}'.format(T), exist_ok=True)
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
        os.makedirs('LSF/{}'.format(T), exist_ok=True)

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

            os.makedirs('lists', exist_ok=True)
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
            lines += self.init_GATK_cmd(analysis_type, memMB,)

            lines += [' --num_threads {} \\'.format(num_threads)]

            ## required, in
            lines += [' --input {} \\'.format(vcf) for vcf in l_vcfs_in]
            ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_variantrecalibration_VariantRecalibrator.html#--use_annotation
            ## http://gatkforums.broadinstitute.org/discussion/2805/howto-recalibrate-variant-quality-scores-run-vqsr
            if mode == 'SNP':
                lines += [
                    ' -an DP \\',
                    ' -an QD \\',
                    ' -an FS \\',
                    ' -an SOR \\',
                    ' -an MQ \\',
                    ' -an MQRankSum \\',
                    ' -an ReadPosRankSum \\',
##                    ' -an InbreedingCoeff \\',
                    ]
            elif mode == 'INDEL':
                lines += [
                    ' -an QD \\',
                    ' -an DP \\',
                    ' -an FS \\',
                    ' -an SOR \\',
                    ' -an MQRankSum \\',
                    ' -an ReadPosRankSum \\',
##                    ' -an InbreedingCoeff \\',
                    ]

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

            lines += ['"\n\neval $cmd']

            ## Terminate command and rerun pipeline.
            lines += self.term_cmd(
                '{}.{}'.format(analysis_type,mode),[fp_tranches,fp_recal,],)

            self.write_shell('shell/{}.{}.sh'.format(analysis_type,mode,), lines,)

            LSB_JOBNAME = 'VR.{}'.format(mode)
            cmd = self.bsub_cmd(
                '{}.{}'.format(analysis_type, mode), LSB_JOBNAME, memMB=memMB, queue=queue,
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
        if os.path.isfile(fp_out):
            sys.exit()
        with open('out_VariantRecalibrator/SNP.recal') as fd_recal_SNP, \
             open('out_VariantRecalibrator/INDEL.recal') as fd_recal_INDEL, \
             gzip.open(self.args.AR_input, 'rt') as fd_source, \
             BgzfWriter(fp_out, 'wb') as fd_out:
            ## write meta-information header
            print('##fileformat=VCFv4.1', file=fd_out)
            print('##filedate={}'.format(
                datetime.datetime.now().strftime("%Y%m%d")), file=fd_out)
            print('##source="GATK_pipeline.py"', file=fd_out)
            chrom_SNP = chrom_INDEL = None
            for line_VCF in fd_source:
                if line_VCF[:2] == '##':
                    continue
                assert line_VCF[:6] == '#CHROM'
                ## write sample IDs to output
                print(line_VCF, end='', file=fd_out)
                break
            while chrom_SNP != chrom:
                chrom_SNP, pos_SNP, VQSLod_SNP = self.parse_recal(
                    fd_recal_SNP, pattern)
            while chrom_INDEL != chrom:
                chrom_INDEL, pos_INDEL, VQSLod_INDEL = self.parse_recal(
                    fd_recal_INDEL, pattern)
            for line_VCF in fd_source:
                chrom_VCF, pos_VCF = line_VCF.split('\t', 2)[:2]
                assert chrom == chrom_VCF
                print(chrom, pos_VCF, file=sys.stderr)
                if pos_VCF == pos_INDEL:
                    assert chrom_VCF == chrom_INDEL
                    if VQSLod_INDEL >= d_minVQSLod['INDEL']:
                        print(line_VCF, end='', file=fd_out)
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

        ## return and continue with BEAGLE if all AR processes completed
        return


    def BEAGLE4(self):

        ## http://faculty.washington.edu/browning/beagle

        memMB = 12900  # todo: move to argparse
        window = 50000  # todo: move to argparse
        if self.checkpoint:
            queue = 'normal'  # todo: move to argparse
            nthreads = 1
        else:
            queue = 'basement'
            nthreads = 1
            queue = 'long'
            nthreads = 4
            nthreads = 12

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
        lines += ['pos1=$3']
        lines += ['pos2=$4']
        lines += ['LSB_JOBINDEX=$2']
        lines += ['out=$out']
        lines += ['mkdir -p $(dirname $out)']
        ## exit if output already exists
        lines += ['if [ -s $out.vcf.gz ]; then exit; fi']
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
        lines += [' nthreads={:d} \\'.format(nthreads)]
        lines += [' window={} \\'.format(window)]  # default 50000 as of r1274
        lines += [' overlap=3000 \\']  # default 3000
        lines += [' gprobs=true \\']
        lines += [' usephase=false \\']
        lines += [' seed=-99999 \\']
        lines += [' singlescale=0.8 \\']  # default 0.8 as of r1389
        lines += [' duoscale=1.0 \\']  # default 1.0
        lines += [' trioscale=1.0 \\']  # default 1.0
        lines += [' burnin-its=5 \\']  # default 5
        lines += [' phase-its=5 \\']  # default 5
        lines += [' impute-its=5 \\']  # default 5
        ## Advanced options not recommended for general use
        lines += [' nsamples=4 \\']  # default 4
        lines += [' buildwindow=1200 \\']  # default 1200 as of r1274
        ## term cmd
        lines += self.term_cmd(
            'BEAGLE', ['$out.vcf.gz'], extra='tabix -p vcf $out.vcf.gz')
        ## write shell script
        if not os.path.isfile('shell/BEAGLE.sh'):
            self.write_shell('shell/BEAGLE.sh',lines,)
        if self.checkpoint:
            self.write_brestart()
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
                    variables = []
                    variables += ['out=out_BEAGLE/{}/{}'.format(chrom, index)]
                    variables += ['chrom={}'.format(chrom)]
                    variables += ['pos1={}'.format(pos1)]
                    variables += ['pos2={}'.format(pos2)]
                    variables = ' '.join(variables)
                    self.bsub_BEAGLE(
                        chrom, pos1, pos2, index, memMB, queue, nthreads)
                    print(chrom, ':', pos1, '-', pos2, index, cnt)
                    pos = pos_prev = None
                else:
                    pass
                pos_prev = pos
                continue

            pos2 = pos
            index = (cnt//window)+1

            variables = []
            variables += ['out=out_BEAGLE/{}/{}'.format(chrom, index)]
            variables += ['chrom={}'.format(chrom)]
            variables += ['pos1={}'.format(pos1)]
            variables += ['pos2={}'.format(pos2)]
            variables = ' '.join(variables)

            self.bsub_BEAGLE(
                LSF_memMB=memMB, LSF_queue=queue, LSF_n=nthreads,
                variables=variables)

        return


    def bsub_BEAGLE(self, chrom, pos1, pos2, index, memMB, queue, nthreads):

        fn_out = 'out_BEAGLE/{}/{}.vcf.gz'.format(
            chrom, index)

        ## finished?
        if os.path.isfile('{}.tbi'.format(fn_out)):
            return

        ## started and running?
        fn = 'LSF/BEAGLE/{}.{}.out'.format(chrom, index)
        if os.path.isfile(fn):
            print('a', fn)
            if os.path.getsize(fn):
                print('b', fn)
                ## running?
                with open(fn) as f:
                    line = f.readlines()[-1]
                    print(line)
                    if 'mean edges' in line or 'iterations' in line:
                        return
            else:
                stop

        ## started and finished? ## tmp!!!
        fn = 'out_BEAGLE/{}/{}.log'.format(chrom, index)
        if os.path.isfile(fn):
            print('a', fn)
            if os.path.getsize(fn):
                print('b', fn)
                ## running?
                with open(fn) as f:
                    line = f.readlines()[-1]
                    print(line)
                    if line.rstrip().split()[-1] == 'finished':
                        stopshouldnothappen ## tmp!!!
##                        subprocess.call(
##                            'tabix -p vcf out_BEAGLE/{}/{}.vcf.gz'.format(
##                                chrom, index), shell=True)
                        return
            else:
                stop

        print(chrom,index)

        LSB_JOBNAME = '{}.{}.{}'.format('BEAGLE', chrom, index,)
        LSF_affix = '{}/{}.{}'.format('BEAGLE', chrom, index)
        cmd_BEAGLE = self.bsub_cmd(
            'BEAGLE', LSB_JOBNAME, memMB=memMB, LSF_affix=LSF_affix,
            chrom=chrom, queue=queue, pos1=pos1, pos2=pos2, index=index,
            num_threads=nthreads)

        if self.checkpoint:
            s = subprocess.check_output(cmd_BEAGLE, shell=True).decode()
            print(s)
            jobID = int(re.match('.*?<(.*?)>',s).group(1))
            print(jobID)
            cmd_brestart = 'bsub -G %s' %(self.project)
            cmd_brestart += ' -o brestart.out -e brestart.err'
            cmd_brestart += ' -q small -w "ended(%i)"' %(jobID)
            cmd_brestart += ' bash shell/brestart.sh %i %s %i' %(
                jobID, self.project, memMB)
            print(cmd_brestart)
            print()
            subprocess.call(cmd_brestart, shell=True)
        else:
            subprocess.call(cmd_BEAGLE, shell=True)

        return


    def bsub_cmd(
        self,
        ## variables outside script
        shell_affix, LSB_JOBNAME,
        LSF_queue='normal', LSF_memMB=4000, LSF_affix=None, LSF_n=1,
        ## variables inside script
        variables='',
        ):

        if not LSF_affix:
            LSF_affix = shell_affix

        cmd = 'bsub -J"{}" -q {}'.format(LSB_JOBNAME, LSF_queue)
        cmd += ' -G {}'.format(self.project)
        cmd += " -M%i -R'select[mem>%i] rusage[mem=%i]'" %(
            LSF_memMB, LSF_memMB, LSF_memMB,)
        cmd += ' -o {}/LSF/{}.out'.format(os.getcwd(), LSF_affix)
        cmd += ' -e {}/LSF/{}.err'.format(os.getcwd(), LSF_affix)
        if LSF_n > 1:
            cmd += ' -n{:d} -R"span[hosts=1]"'.format(LSF_n)
        if self.checkpoint:
            cmd += ' -k "{} method=blcrkill 600"'.format(
                os.path.join(os.getcwd(), 'checkpoint'))
            cmd += ' -r'
        if self.checkpoint:
            cmd += ' cr_run'
##        for k, v in {
##            'chrom':chrom, 'index':index, 'bam':bam, 'mode':mode,
##            'pos1': pos1, 'pos2':pos2}.items():
##            if v:
##                cmd += ' {}={}'.format(k, v)
        cmd += ' {}'.format(variables)
        bash_cmd = ' bash {}/shell/{}.sh'.format(os.getcwd(), shell_affix)
        cmd += bash_cmd
##        for x in (chrom, index, bam, mode, pos1, pos2):
##            if x:
##                cmd += ' {}'.format(x)

        return cmd


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

        lines += ['"\n\neval $cmd']

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
        lines += [' -nt 8 \\'.format(T)]
        lines += [' --annotation InbreedingCoeff \\']  # default
        lines += [' --annotation FisherStrand \\']  # default in 3.3?
        lines += [' --annotation StrandOddsRatio \\']  # default in 3.3?
        lines += [' --annotation QualByDepth \\']  # default
        lines += [' --annotation ChromosomeCounts \\']  # default
        lines += [' --annotation GenotypeSummaries \\']  # default
        lines += [' --annotation MappingQualityRankSumTest \\']
        lines += [' --annotation ReadPosRankSumTest \\']
        lines += [' --standard_min_confidence_threshold_for_calling 30 \\']
        lines += [' --standard_min_confidence_threshold_for_emitting 30 \\']

        lines += ['"\n\neval $cmd']

        ## terminate shell script
        lines += self.term_cmd(T, ['$out.tbi'])

        ## write shell script
        self.write_shell('shell/{}.sh'.format(T), lines,)

        return



    def execmd(self,cmd):

        print(cmd)
        subprocess.call(cmd,shell=True)

        return


    def write_shell(self,fp,lines,):

        os.makedirs(os.path.dirname(fp), exist_ok=True)

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


    def init_java(self, jar, memMB, java='java'):

        s = '{} -Djava.io.tmpdir={}'.format(java, 'tmp')
        ## set maximum heap size
        s += ' -Xmx{}m'.format(memMB)
        if self.checkpoint:
            s += ' -XX:-UsePerfData -Xrs '
        s += ' \\\n -jar {}'.format(jar)

        return s


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

        with open('shell/brestart.sh', 'w') as f:
            f.write('sleep 30\n')
            ## internal field separator
            f.write("IFS=$'\\n'\n")
            f.write('jobID=$1\n')
            f.write('project=$2\n')
            f.write('memMB=$3\necho memMB $memMB\n')
            f.write('pwd=$(pwd)\n')
            ## parse bhist
            f.write('bhist=$(bhist -l $jobID)\n')
            ## Checkpoint succeeded
            f.write('''cpsucc=$(echo $bhist | sed 's/ *//g' | grep Checkpointsucceeded | wc -l)\n''')
            ##  exit code 13, TERM_CHKPNT, Job killed after checkpointing
            f.write('''exit13=$(echo $bhist | sed 's/ *//g' | grep "Exitedwithexitcode13" | wc -l)\n''')  # could also be 13x
            ##  exit code 143, SIGTERM
            f.write('''exit143=$(echo $bhist | sed 's/ *//g' | grep "Exitedwithexitcode143" | wc -l)\n''')
            ## exit code 140, run limit
            f.write('''exit140=$(echo $bhist | grep TERM_RUNLIMIT | wc -l)\n''')
            ## exit code 16, pid taken
            f.write('''exit16=$(echo $bhist | sed 's/ *//g'| grep Exitedwithexitcode16 | wc -l)\n''')
            ## Checkpoint failed
            f.write('''cpfail=$(echo $bhist | sed 's/ *//g'| grep "Checkpointfailed" | wc -l)\n''')
            ## Done successfully
            f.write('''donesuc=$(echo $bhist | sed 's/ *//g'| grep "Donesuccessfully" | wc -l)\n''')
            ## exit if done succesfully
            f.write('if [ $donesuc -eq 1 ]; then echo $bhist >> bhist_success.tmp; exit; fi\n')
            ## exit if not checkpoint succeeded and not PID taken
            f.write('if [ $exit143 -eq 0 -a $cpsucc -eq 0 -a $exit13 -eq 0 -a $exit16 -eq 0 ]; then echo $bhist >> bhist_unexpectederror.tmp; exit; fi\n')
            ## restart job and capture jobID
            f.write('s=$(brestart -G $project -M$memMB $pwd/checkpoint/$jobID)\n')
            f.write('''jobID=$(echo $s | awk -F "[<>]" '{print $2}')\n''')
            ## report if checkpoint failed
            f.write('if [ $cpfail -ne 0 ]; then echo $s >> checkpointfailed_brestartout.tmp; fi\n')
            ## be verbose
            f.write('echo s $s\n')
            f.write('echo jobID $jobID\n')
            f.write('echo memMB $memMB\n')
            ## bsub this chaperone restart script again
            f.write("bsub -R 'select[mem>'$memMB'] rusage[mem='$memMB']' -M$memMB \\\n")
            f.write(' -o brestart.out -e brestart.err \\\n')
            f.write(' -G $project -q normal -w "ended($jobID)" \\\n')
            f.write(' bash shell/brestart.sh $jobID $project $memMB\n')

        return


##    def bash_chrom2len(self, d_chrom_ranges, size_fragment_bp):
##
##        l_chroms = list(sorted(d_chrom_ranges.keys()))
##        lines = []
##        lines += ['\n## define arrays']
##        lines += ['chroms=(%s)' %(' '.join(l_chroms))]
##        lines += ['starts=(%s)' %(' '.join([
##            str(d_chrom_ranges[chrom][0]) for chrom in l_chroms],),)]
##        lines += ['stops=(%s)' %(' '.join([
##            str(d_chrom_ranges[chrom][1]) for chrom in l_chroms],),)]
##
##        lines += ['\n## find chromosome index in array of chromosomes']
##        lines += ['for ((index=0; index<${#chroms[@]}; index++))\ndo']
##        lines += ['if [ "${chroms[$index]}" == "$chrom" ]\nthen']
##        lines += ['break\nfi\ndone\n']
##
##        lines += ['\n## find length of chromosome passed via the command line']
##        lines += ['start=${starts[$index]}\n']
##        lines += ['stop=${stops[$index]}\n']
##
##        ##
##        ## do not allow interval to exceed the length of the chromosome
##        ## otherwise it will raise an I/O error
##        ##
##        lines += ['posmax=$(($start-1+(${LSB_JOBINDEX}+0)*%i))' %(size_fragment_bp)]
##        lines += ['echo $posmax $start $LSB_JOBINDEX; exit']
##        lines += ['if [ $posmax -gt $stop ]']
##        lines += ['then posmax=$lenchrom']
##        lines += ['fi\n']
##
##        return lines


    def shell_UG(self, T, memMB, nct, nt):

        lines = ['#!/bin/bash\n']

##        lines += self.bash_chrom2len(d_chrom_ranges, size_fragment_bp)

##        ## exit if job started
##        lines += ['if [ -s $out ]; then exit; fi\n']
        ## exit if job finished
        lines += ['if [ -s $out.tbi ]; then exit; fi\n']

        ## initiate GATK command
        lines += self.init_GATK_cmd(T, memMB)

        ## append GATK command options
        lines += self.body_UnifiedGenotyper()

        ## terminate shell script
##        lines += self.term_cmd(
##            analysis_type, ['$out.tbi'], extra='tabix -p vcf $out')
        lines += self.term_cmd(T, ['$out.tbi'])

        ## write shell script
        self.write_shell('shell/{}.sh'.format(T),lines)

        return


    def body_UnifiedGenotyper(self):

        '''Write walker specific (non command line) command line arguments.'''

        ## https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_genotyper_UnifiedGenotyper.php
        ## https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_engine_CommandLineGATK.php

        lines = []

        ##
        ## CommandLineGATK
        ##

        ## https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_engine_CommandLineGATK.php#--input_file
        lines += [' --input_file $input_file \\']

        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_CommandLineGATK.html#--intervals
        lines += [' --intervals $chrom:$pos1-$pos2 \\']
        s = '"\nif [ $XL != "" ]; then cmd=$cmd"'
        s += ' --excludeIntervals $XL"; fi\ncmd=$cmd" \\'
        lines += [s]
        if self.intervals:
            lines += ['--intervals {} \\'.format(self.intervals)]
            ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_CommandLineGATK.html#--interval_set_rule
            lines += ['--interval_set_rule INTERSECTION \\']
            pass

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

        ## https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_genotyper_UnifiedGenotyper.php#--annotation
        s_annotation = ''
        ## https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_annotator_Coverage.php
        s_annotation += ' -A Coverage'
        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_annotator_FisherStrand.html
        s_annotation += ' -A FisherStrand'
        ## https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_annotator_StrandOddsRatio.php
        s_annotation += ' -A StrandOddsRatio'
        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_annotator_MappingQualityRankSumTest.html
        s_annotation += ' -A MappingQualityRankSumTest'
        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_annotator_QualByDepth.html
        s_annotation += ' -A QualByDepth'
        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_annotator_RMSMappingQuality.html
        s_annotation += ' -A RMSMappingQuality'
        ## https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_annotator_ReadPosRankSumTest.php
        s_annotation += ' -A ReadPosRankSumTest'
        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_annotator_HaplotypeScore.html
        s_annotation += ' -A HaplotypeScore'
        ## https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_annotator_InbreedingCoeff.php
        s_annotation += ' -A InbreedingCoeff'
        lines += [' {} \\'.format(s_annotation)]

        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_genotyper_UnifiedGenotyper.html#--genotype_likelihoods_model
        lines += [' --genotype_likelihoods_model BOTH \\'] ## default value SNP

        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_haplotypecaller_HaplotypeCaller.html#--genotyping_mode
        lines += [' -gt_mode {} \\'.format(self.genotyping_mode)] ## default value DISCOVERY

        lines += [' --sample_ploidy $sample_ploidy \\']

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

        lines += ['"\n\neval $cmd']

        return lines


    def write_bam_and_XL_lists(self, d_chrom_ranges):

        if not os.path.isdir('lists'):
            os.mkdir('lists')
        with open('lists/bams.list', 'w') as f:
            for bam in self.bams:
                f.write('{}\n'.format(bam))
        if self.sex:
            with open(self.sex) as f:
                d_sample2sex = {
                    line.split()[0]:line.split()[1].lower()[0] for line in f}
            d_sex2bam = {'m':[], 'f':[]}
            for bam in self.bams:
                sample = os.path.splitext(os.path.basename(bam))[0]
                try:
                    sex = d_sample2sex[sample]
                except KeyError:
                    continue
                d_sex2bam[sex].append(bam)
            for sex in ('m', 'f'):
                with open('lists/bams.{}.list'.format(sex), 'w') as f:
                    for bam in d_sex2bam[sex]:
                        f.write('{}\n'.format(bam))

        ## Create region lists. One of them empty.
        ## http://gatkforums.broadinstitute.org/discussion/1204/what-input-files-does-the-gatk-accept-require
        with open('lists/XL.PAR.list', 'w') as f:
            for region in ('PAR1','PAR2'):
                f.write('X:{:d}-{:d}\n'.format(
                    d_chrom_ranges[region][0], d_chrom_ranges[region][1]))

        return


    def body_HaplotypeCaller(self,):

        '''Write walker specific command line arguments.'''

        ## https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_haplotypecaller_HaplotypeCaller.html

        lines = []

        ##
        ## Inherited arguments
        ##

        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_CommandLineGATK.html#--input_file
        lines += [' --input_file $input_file \\']

        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_CommandLineGATK.html#--intervals
        lines += [' --intervals $chrom:$pos1-$pos2 \\']
        s = '"\nif [ $XL != "" ]; then cmd=$cmd"'
        s += ' --excludeIntervals $XL"; fi\ncmd=$cmd" \\'
        lines += [s]
        if self.intervals:
            lines += [' --intervals {} \\'.format(self.intervals)]
            ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_CommandLineGATK.html#--interval_set_rule
            lines += [' --interval_set_rule INTERSECTION \\']
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

        ## https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_haplotypecaller_HaplotypeCaller.html#--annotation
        s_annotation = ''
        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_annotator_Coverage.html
        s_annotation += ' -A Coverage'
        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_annotator_FisherStrand.html
        s_annotation += ' -A FisherStrand'
        ## https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_annotator_StrandOddsRatio.php
        s_annotation += ' -A StrandOddsRatio'
        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_annotator_MappingQualityRankSumTest.html
        s_annotation += ' -A MappingQualityRankSumTest'
        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_annotator_QualByDepth.html
        s_annotation += ' -A QualByDepth'
        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_annotator_RMSMappingQuality.html
        s_annotation += ' -A RMSMappingQuality'
        s_annotation += ' -A ReadPosRankSumTest'
        lines += [' {} \\'.format(s_annotation)]

        lines += [' --sample_ploidy $sample_ploidy \\']

        lines += [' --emitRefConfidence GVCF \\']

        lines += ['"\n\neval $cmd']

        return lines


    def get_ploidy(self, chrom, sex):
        
        if chrom == 'Y':
            sample_ploidy = 1
        elif chrom == 'X' and sex == 'f':
            sample_ploidy = 2
        elif chrom == 'X' and sex == 'm':
            sample_ploidy = 1
        else:
            sample_ploidy = 2

        return sample_ploidy


    def get_sex_and_XL(self, chrom):

        if chrom == 'Y':
            d_sex2bamlist = {'m':'lists/bams.m.list'}
            XL = None
        elif chrom == 'X':
            d_sex2bamlist = {
                'm':'lists/bams.m.list',
                'f':'lists/bams.f.list'}
##                elif chrom in ('PAR1', 'PAR2'):
            XL = 'lists/XL.PAR.list'
        else:
            d_sex2bamlist = {'':'lists/bams.list'}
            XL = None

        return d_sex2bamlist, XL


    def init_GATK_cmd(self, analysis_type, memMB):

        ## https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_engine_CommandLineGATK.php

        lines = []

        ## exit if output exists
        lines += ['if [ -s $out ]; then exit; fi']

        ## create output folder
        lines += ['mkdir -p $(dirname $out)']

        s = ''
        ## Java version alert: starting with release 2.6, GATK now requires Java 1.7. See Version Highlights for 2.6 for details.
        ## http://www.broadinstitute.org/gatk/guide/article?id=2846
        s_java = self.init_java(self.fp_GATK, memMB, java=self.java)
        s += ' cmd="{} \\'.format(s_java)
        lines += ['\n{}'.format(s)]

        ## CommandLineGATK, required, in
        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_CommandLineGATK.html#--analysis_type
        lines += [' --analysis_type {} \\'.format(analysis_type)]
        ## CommandLineGATK, optional, in
        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_CommandLineGATK.html#--reference_sequence
        lines += [' --reference_sequence {} \\'.format(self.reference_sequence)]
        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_CommandLineGATK.html#-nct
        lines += [' --num_cpu_threads_per_data_thread $nct \\']
        ## https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_engine_CommandLineGATK.php#--num_threads
        lines += [' --num_threads $nt \\']


        return lines


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
            help='Path to bam file and/or directory containing bam files',
            nargs='+', required=True, type=self.is_file_or_dir)

        parser.add_argument('--coverage', required=True, type=float)

        parser.add_argument(
            '--build', required=True, type=int, choices=[37,38])

        parser.add_argument('--sex', required=False, type=self.is_file)

        parser.add_argument(
            '--fp_GATK', '--GATK', '--gatk', '--jar', required=True,
            help='File path to GATK', type=self.is_file)

        parser.add_argument('--project', required=True)

        parser.add_argument('--arguments', '--args')

        parser.add_argument('--java', required=True, type=self.is_file)

        ##
        ## CommandLineGATK arguments
        ##

        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_CommandLineGATK.html#--intervals
        parser.add_argument(
            '--intervals', '-L',
            help='Additionally, one may specify a rod file to traverse over the positions for which there is a record in the file (e.g. -L file.vcf).',
            )

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

        parser.add_argument('--i_BEAGLE_nsamples', default=20)

        ##
        ## optional arguments
        ##
        parser.add_argument(
            '--checkpoint', action='store_true', default=False)

        parser.add_argument(
            '--chroms', type=str, nargs='+',
            default=[str(i+1) for i in range(22)]+['X', 'Y', 'MT', 'PAR1', 'PAR2'])

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
