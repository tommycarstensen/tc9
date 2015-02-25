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
import gzip
import subprocess
import contextlib
import Bio
from Bio.bgzf import BgzfWriter
import urllib
import datetime
import pysam
import urllib.parse
import operator
import socket


## README

## http://www.broadinstitute.org/gatk/about#typical-workflows
## http://www.broadinstitute.org/gatk/guide/topic?name=best-practices

class main():

    def main(self):

        print('main')

        if self.args.AR_input:
            self.run_ApplyRecalibration()

        if self.coverage >= 15:
            self.HaplotypeCaller()
            self.CombineGVCFs()
            self.GenotypeGVCFs()
            self.caller = 'HC'
        elif self.coverage < 15:
            self.UnifiedGenotyper()
            self.caller = 'UG'

        self.VariantRecalibrator()

        self.ApplyRecalibration()

        self.beagle4()

        return

    def UnifiedGenotyper(self):

        T = analysis_type = 'UnifiedGenotyper'
        ## fragment sizes determines runtime determines queue
        queue = 'normal'
        size_bp = 5 * 10 ** 5  # normal
        queue = 'long'
        size_bp = 2 * 10 ** 6  # long
        ## http://gatkforums.broadinstitute.org/discussion/1975/how-can-i-use-parallelism-to-make-gatk-tools-run-faster
        nct = 3
        nt = 8
        queue = 'normal'
        ## Each data thread needs to be given the full amount of memory
        ## you’d normally give a single run. So if you’re running a tool
        ## that normally requires 2 Gb of memory to run, if you use -nt 4,
        ## the multithreaded run will use 8 Gb of memory. In contrast,
        ## CPU threads will share the memory allocated to their “mother”
        ## data thread, so you don’t need to worry about allocating memory
        ## based on the number of CPU threads you use.

        ## Use all memory if all cores used anyway.
        if nct * nt == 32:
            memMB = 255900
        ## Memory sample count and nt dependent.
        else:
            memMB = min(255900, 47900 + 2250 * nt)

        nct = 3; nt = 8; memMB = 127900
        memMB = 191900
        nct = 1; nt = 1; memMB = 63900; queue = 'normal'

        ## Parse chromosome ranges.
        d_chrom_ranges = self.parse_chrom_ranges()

##        d_Y = {}
##        url = 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/chrY/chrY_callable_regions.20130802.bed'
##        with urllib.request.urlopen(url) as f:
##            for line in f:
##                line = line.decode('utf-8')
##                pos1, pos2 = [int(x) for x in line.split()[1:]]

        ## touch/lock
        if self.touch(analysis_type):
            return
##            print('touch file exists')
##            ## Start retrying failed jobs, once all stderr files are generated.
##            for chrom in self.chroms:
##                cstart = d_chrom_ranges[chrom][0]
##                cstop = d_chrom_ranges[chrom][1]
##                for i in range(
##                    1, 1 + math.ceil((cstop - cstart) / size_bp)):
##                    affix = '{}/{}/{}'.format(T, chrom, i)
##                    # Return if not all jobs started or completed.
##                    if not os.path.isfile('LSF/{}.err'.format(affix)):
##                        return

        ## Write bam lists for
        ## all (autosomes, PAR), male (Y) and female (nonPAR X) samples.
        if self.sample_genders:
            self.write_bam_and_XL_lists(d_chrom_ranges)

        ## Create folders.
        os.makedirs('tmp', exist_ok=True)

        ## write shell script
        self.shell_UG(T, memMB, nct, nt)

        ## execute shell script
        ## loop over chroms
        s_out = ''
        for chrom in self.chroms:

            print('UG', chrom)

            if chrom in ('MT', 'X', 'Y', 'PAR1', 'PAR2') and self.sample_genders:
                nct = 1; nt = 1; queue = 'long'

            d_sex2bamlist, XL = self.get_sex_and_XL(chrom)
            ## Memory consumption seems to explode for UG,
            ## when doing haploidy for Y and nonPAR X.
            ## Therefore do diploidy.
            if not self.sample_genders:
                d_sex2bamlist = {'': 'lists/bams.list'}
                XL = None

            ## Look up chromosome ranges.
            cstart = d_chrom_ranges[chrom][0]
            cstop = d_chrom_ranges[chrom][1]

            for sex in d_sex2bamlist.keys():

                if not os.path.isfile(d_sex2bamlist[sex]):
                    stop1
                    continue
                if not os.path.getsize(d_sex2bamlist[sex]):
                    stop2
                    continue

                ## See comment above about memory and haploidy.
                sample_ploidy = self.get_ploidy(chrom, sex)
                if not self.sample_genders:
                    sample_ploidy = 2

                for i in range(
                    cstart // size_bp + 1,
                    cstart // size_bp + 1 + math.ceil((cstop - cstart) / size_bp)):

##                    pos1 = 1 + cstart // size_bp + (i - 1) * size_bp
##                    pos2 = min(cstop, (pos1 - 1) + size_bp)

##                    pos1 = max(cstart, cstart - cstart % size_bp + (i - 1) * size_bp)
                    pos1 = max(cstart, (i - 1) * size_bp + 1)
                    pos2 = min(cstop, i * size_bp)
##                    if chrom in [str(i) for i in range(1,23)]:
##                        continue  # tmp!!!

                    affix = '{}/{}/{}'.format(T, chrom, i)
                    if sex and chrom == 'X':
                        affix += '.{}'.format(sex)

                    out = 'out_{}.vcf.gz'.format(affix)
                    s_out += '{}\n'.format(out)

                    ## Skip if output was generated.
                    if os.path.isfile(
                        'out_{}.vcf.gz.tbi'.format(affix)):
                        continue

                    ## Skip if output is being generated.
                    if os.path.isfile('LSF/{}.err'.format(affix)):
                        if time.time() - os.path.getmtime(
                            'LSF/{}.err'.format(affix)) < 300:
                            print('err recently modified', affix)
                            continue
                        ## Otherwise delete prior to run.
                        if os.path.isfile('LSF/{}.out'.format(affix)):
                            os.remove('LSF/{}.out'.format(affix))
                        os.remove('LSF/{}.err'.format(affix))

                    ## Delete old output.
                    if (
                        not os.path.isfile('out_{}.vcf.gz.tbi'.format(affix))
                        and os.path.isfile('out_{}.vcf.gz'.format(affix)) and
                        time.time() - os.path.getmtime(
                            'out_{}.vcf.gz'.format(affix)) > 15 * 60
                        ):
                        print(out)
                        os.remove(out)

                    os.makedirs(os.path.dirname(
                        'LSF/{}'.format(affix)), exist_ok=True)
                    os.makedirs(os.path.dirname(
                        'out_{}'.format(affix)), exist_ok=True)

                    d_arguments = {}
                    if chrom in ('PAR1', 'PAR2'):
                        d_arguments['chrom'] = 'X'
                    else:
                        d_arguments['chrom'] = chrom
                    d_arguments['pos1'] = pos1
                    d_arguments['pos2'] = pos2
                    if XL:
                        d_arguments['XL'] = XL
                    d_arguments['out'] = out
                    d_arguments['input_file'] = d_sex2bamlist[sex]
                    d_arguments['nct'] = nct
                    d_arguments['nt'] = nt
                    d_arguments['sample_ploidy'] = sample_ploidy
                    arguments = self.args_dict2str(d_arguments)

                    LSB_JOBNAME = '{}.{}.{} {}'.format('UG', chrom, i, sex)
                    cmd = self.bsub_cmd(
                        T, LSB_JOBNAME, LSF_affix=affix,
                        LSF_memMB=memMB, LSF_queue=queue, LSF_n=nt * nct,
                        arguments=arguments)
                    self.execmd(cmd)

        with open('lists/{}.list'.format(T), 'w') as f:
            f.write(s_out)

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
                        if time.time() - os.path.getmtime(
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

##        ## exit if job started
##        lines += ['if [ -f $out ]; then exit; fi\n']
        ## exit if job finished
        lines += ['if [ -s $out.tbi ]; then exit; fi\n']

        ## initiate GATK command
        lines += self.init_GATK_cmd(T, memMB)

        ## append GATK command options
        lines += self.body_HaplotypeCaller()

        ## terminate shell script
        lines += self.term_cmd(T, ['$out.tbi'])

        ## write shell script
        self.write_shell('shell/{}.sh'.format(T), lines,)

        return

    def parse_chrom_ranges(self):

        d_chrom_ranges = {}

        ## 1000G chromosome ranges
        fn = '{}.fai'.format(self.reference_sequence)
        fd = open(fn)
        lines = fd.readlines()
        fd.close()
        for line in lines:
            l = line.strip().split()
            chrom = l[0]
            chrom_len = int(l[1])
##            if not chromosome in ('X','Y'):
##                chrom = int(chromosome)
            d_chrom_ranges[chrom] = (1, chrom_len)
##            ## break, when we reach the Y chromosome
##            if chrom == 'Y':
##                break

        if self.build == 37:
            url = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All'
            url += '/GCF_000001405.25.regions.txt'
        elif self.build == 38:
            url = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All'
            url += '/GCF_000001405.28.regions.txt'
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
                    k = l[0].replace('#', '')
                    v = (int(l[2]), int(l[3]))
                    d_chrom_ranges[k] = v

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

            l_vcfs_in = list(sorted(list(set(l_vcfs_in) - set(l_combined))))

            for i, vcf in enumerate(
                l_vcfs_in,
                self.gVCF_limit * len(glob.glob(
                    'lists/{}.{}.*.list'.format(T, chrom)))):
                if i % self.gVCF_limit == 0:
                    fn_out = 'lists/{T}.{chrom}.{i}.list'.format(
                        T=T, chrom=chrom, i=i // self.gVCF_limit)
                    assert not os.path.isfile(fn_out)
                    fd_out = open(fn_out, 'w')
                fd_out.write('{}\n'.format(vcf))
            fd_out.close()

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
                cmd += ' {}.{}'.format(chrom, i)
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

                s_out += 'out_{}/{}.vcf.gz\n'.format(T, chrom)

        with open('lists/{}.list'.format(T), 'w') as f:
            f.write(s_out)

        if bool_exit:
            sys.exit()

        return

    def assert_identical_headers(self, l_vcfs):

        ## Assert that all headers are identical.
        for i, vcf in enumerate(l_vcfs):
            with gzip.open(vcf, 'rt') as fd_source:
                for line_vcf in fd_source:
                    if line_vcf[:2] == '##':
                        continue
                    if i == 0:
                        l = line_vcf.rstrip().split()
                    else:
                        assert l == line_vcf.rstrip().split()
                    break

        return

    def VariantRecalibrator(self):

        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_variantrecalibration_VariantRecalibrator.html

        T = analysis_type = 'VariantRecalibrator'
        num_threads = 4
        queue = 'yesterday'
        num_threads = 4
        num_threads = 24
        queue = 'normal'

        ##
        ## 1) check input existence (vcf)
        ##
        if self.caller == 'HC':
            T_prev = 'GenotypeGVCFs'
##            l_vcfs_in = [
##                'out_GenotypeGVCFs/{}.vcf.gz'.format(chrom) for chrom in self.chroms]
            with open('lists/{}.list'.format(T_prev)) as f:
                l_vcfs_in = f.read().rstrip().split('\n')
        elif self.caller == 'UG':
            T_prev = 'UnifiedGenotyper'
            with open('lists/{}.list'.format(T_prev)) as f:
                l_vcfs_in = f.read().rstrip().split('\n')
        if self.check_in(
            T_prev, ['{}'.format(vcf) for vcf in l_vcfs_in],
            'touch/{}.touch'.format(T_prev)):
            sys.exit(0)

        assert len(l_vcfs_in) > 0
        for file in l_vcfs_in:
            if file[-7:] != '.vcf.gz':
                print(file, 'has wrong extension')
                sys.exit()

##        self.assert_identical_headers(l_vcfs_in)

        d_resources = {
            'SNP': self.resources_SNP, 'INDEL': self.resources_INDEL}

        ## Check that resources are present.
        for mode in ('SNP', 'INDEL'):
            with open(d_resources[mode]) as f:
                for line in f:
                    resource = line.split()[-1]
                    print(mode, resource)
                    if not os.path.isfile(resource):
                        print(resource, 'not found')
                        sys.exit(0)

        if not os.path.isdir('LSF/{}'.format(T)):
            os.mkdir('LSF/{}'.format(T))

        d_arguments = {'nt': num_threads, 'nct': 1}

        bool_exit = False

        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_variantrecalibration_VariantRecalibrator.html#--mode
        for mode in ('SNP', 'INDEL'):

            memMB = {'SNP': 191900, 'INDEL': 127900}[mode]

            ## 2) touch / check output
            bool_continue = self.touch('{}.{}'.format(T, mode))
            if bool_continue == True:
                continue

            ## Define file paths.
            tranches = 'out_{}/{}.tranches'.format(T, mode)
            recal = 'out_{}/{}.recal.gz'.format(T, mode)
            if os.path.isfile(tranches) or os.path.isfile(recal):
                continue

            lines = []

            d_arguments['out'] = tranches
            arguments = self.args_dict2str(d_arguments)

            ## Initiate GATK walker.
            lines += self.init_GATK_cmd(
                analysis_type, memMB, d_arguments.keys())

            ## required, in
##            lines += [' --input {} \\'.format(vcf) for vcf in l_vcfs_in]
            lines += [' --input lists/{}.list'.format(T_prev)]
            ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_variantrecalibration_VariantRecalibrator.html#--use_annotation
            ## http://gatkforums.broadinstitute.org/discussion/2805/howto-recalibrate-variant-quality-scores-run-vqsr
            if mode == 'SNP':
                for an in self.an_SNP:
                    lines += [' -an {} \\'.format(an)]
            elif mode == 'INDEL':
                for an in self.an_indel:
                    lines += [' -an {} \\'.format(an)]

            ##
            ## required, out
            ##
            lines += [' --recal_file {} \\'.format(recal)]
            lines += [' --tranches_file {} \\'.format(tranches)]

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
            lines += [
                ' {} \\'.format(line.strip()) for line in lines_resources]

            ## https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_variantrecalibration_VariantRecalibrator.php#--TStranche
            l_TStranches = [100, 99.9, 99.8, 99.5, 99.0, 98.0, 95.0, 90.0]
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
                '{}.{}'.format(analysis_type, mode), [tranches, recal])

            self.write_shell('shell/{}.{}.sh'.format(T, mode,), lines)

            LSB_JOBNAME = 'VR.{}'.format(mode)
            cmd = self.bsub_cmd(
                '{}.{}'.format(analysis_type, mode), LSB_JOBNAME,
                LSF_queue=queue,
                LSF_affix='{}/{}'.format(T, mode),
                LSF_n=num_threads, LSF_memMB=memMB, arguments=arguments)
            self.execmd(cmd)

            bool_exit = True

        if bool_exit:
            sys.exit()

        return

    def args2getopts(self, l_args):

        s = ''
        s += 'string_arguments=$@\n'
        s += 'array_arguments=($string_arguments)\n'
        s += 'i=0\n'
        s += 'for argument in $string_arguments; do \n'
        s += ' i=$(($i+1))\n'
        s += ' case $argument in\n'
        for arg in l_args:
            s += '  --{}) {}=${{array_arguments[$i]}};;\n'.format(arg, arg)
        s += ' esac\n'
        s += 'done\n'

##        s += 'while getopts ":{}:" o; do\n'.format(':'.join(l_args))
##        s += '    case "${o}" in\n'
##        for arg in l_args:
##            s += '        {})\n'.format(arg)
##            s += '            {}=${{OPTARG}}\n'.format(arg)
##            s += '            ;;\n'
##        s += '    esac\n'
##        s += 'done\n'

        return s

    def hook_compressed_text(self, filename, mode):

        ext = os.path.splitext(filename)[1]
        if ext == '.gz':
            f = gzip.open(filename, mode + 't')
        else:
            f = open(filename, mode)

        return f

    def parse_recal(self, fd, pattern):

        for line in fd:
            if line[0] == '#':
                continue
            l = line.rstrip().split('\t')
            chrom = l[0]
            pos = l[1]
            VQSLod = re.match(pattern, l[7]).group(1)
            if VQSLod == 'Infinity':
                VQSLod = float('inf')
            else:
                VQSLod = float(VQSLod)

            yield chrom, pos, VQSLod

    def parse_sources(self):

        d_sources = {}
        for mode in ('SNP', 'INDEL'):
            with gzip.open(
                'out_VariantRecalibrator/{}.recal.gz'.format(
                    mode), 'rt') as f:
                for line in f:
                    if line.split('=')[0] != '##GATKCommandLine':
                        continue
                    pattern = '''##GATKCommandLine=<.*?,'''
                    pattern += '''CommandLineOptions=".*?'''
                    pattern += ''' input=\[\(RodBindingCollection'''
                    pattern += ''' \[(.*?)(?=])'''
                    RodBindings = re.match(pattern, line).group(1)
                    sources = re.findall('source=([\w\.\/]+)', RodBindings)
                    d_sources[mode] = sources
                    break

        return d_sources

    def parse_minVQSLods(self):

        d_ts_filter_level = {
            'SNP': self.ts_SNP,
            'INDEL': self.ts_INDEL}
        d_minVQSLod = {}
        for mode in ('SNP', 'INDEL'):
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

        d_chrom2sources = {}
        min_mtime = min([
            os.path.getmtime(
                'out_VariantRecalibrator/{}.{}'.format(mode, suffix))
            for suffix in ('recal.gz', 'tranches')
            for mode in ('SNP', 'INDEL')])
        for source_SNP in d_sources['SNP']:
            if self.caller == 'UG':
                chrom = os.path.basename(os.path.dirname(source_SNP))
            else:
                chrom = os.path.basename(source_SNP).split('.')[0]
            try:
                d_chrom2sources[chrom] += [source_SNP]
            except KeyError:
                d_chrom2sources[chrom] = [source_SNP]
            assert min_mtime > os.path.getmtime(source_SNP)
        for chrom in self.sort_nicely(list(d_chrom2sources.keys())):
            print('AR', chrom)
            assert chrom in [str(i) for i in range(1, 23)] + ['X', 'Y', 'MT']
            dirname_LSF = 'LSF/ApplyRecalibration'
            if not os.path.isdir(dirname_LSF):
                os.makedirs(dirname_LSF)
            s = ''
            s += 'bsub -G {} '.format(self.project)
            s += ' -o {}/{}.out'.format(dirname_LSF, chrom)
            s += ' -e {}/{}.err'.format(dirname_LSF, chrom)
            s += ' -J AR{}'.format(chrom)
            self.args.AR_input = ' '.join(d_chrom2sources[chrom])
            s += self.args_rerun()
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
                    'out_VariantRecalibrator/{}.recal.gz'.format(mode),
                    'out_VariantRecalibrator/{}.tranches'.format(mode)],
                'touch/VariantRecalibrator.{}.touch'.format(mode)):
                sys.exit()

        ## check output existence
        if self.touch(T):
            return

        fn_touch = 'touch/{}.touch'.format(analysis_type)
        self.bsub_ApplyRecalibration()

        return

    def run_ApplyRecalibration(self):

        ## todo20150112: tc9: use heapq.merge() on SNP+INDEL line generators
        ## todo20150208: tc9: write VQSLOD to INFO field

##        chrom = os.path.basename(self.args.AR_input).split('.')[0]
        chrom = chrom_VCF = os.path.basename(os.path.dirname(
            self.args.AR_input[0]))
        if not chrom in [str(i) for i in range(1, 23)] + ['X', 'Y', 'MT']:
            chrom = chrom_VCF = os.path.basename(
                self.args.AR_input[0]).split('.')[0]
        assert chrom in [str(i) for i in range(1, 23)] + ['X', 'Y', 'MT']

##        index = int(os.path.basename(self.args.AR_input).split('.')[0])

        out = 'out_ApplyRecalibration/{}.vcf.gz'.format(chrom)
        ## Don't overwrite existing files!
        if os.path.isfile(out):
            sys.exit()
        os.makedirs(os.path.dirname(out), exist_ok=True)

        d_minVQSLod = self.parse_minVQSLods()

        pattern = re.compile(r'.*VQSLOD=([-\d\.\w]+);')

        self.assert_identical_headers(self.args.AR_input)

        ## Open input files and output file.
        with gzip.open(
            'out_VariantRecalibrator/SNP.recal.gz', 'rt') as fd_recal_SNP, \
            gzip.open(
                'out_VariantRecalibrator/INDEL.recal.gz',
                'rt') as fd_recal_INDEL, \
                BgzfWriter(out, 'wb') as fd_out:
            ## write meta-information header
            print('##fileformat=VCFv4.2', file=fd_out)
            print('##fileDate={}'.format(
                datetime.datetime.now().strftime("%Y%m%d")), file=fd_out)
            print('##source={}'.format(
                sys.argv[0]), file=fd_out)
            ## Skip lines in recal files preceding first position.
            chrom_SNP = None
            while chrom_SNP != chrom_VCF:
                chrom_SNP, pos_SNP, VQSLod_SNP = next(self.parse_recal(
                    fd_recal_SNP, pattern))
            chrom_INDEL = None
            while chrom_INDEL != chrom_VCF:
                chrom_INDEL, pos_INDEL, VQSLod_INDEL = next(self.parse_recal(
                    fd_recal_INDEL, pattern))
            ## Loop over input files.
            for i, source in enumerate(self.sort_nicely(self.args.AR_input)):
                ## Open input file.
                with gzip.open(source, 'rt') as fd_source:
                    ## todo: 2015jan28: do heapq.merge() on SNPs and INDELs
                    ## when Python3.5 is released
                    ## Skip metainformation lines and header line
                    ## and write it to the output from the first input file.
                    for line_VCF in fd_source:
                        if line_VCF[:2] == '##':
                            ## Copy metadata from first file.
                            if i == 0:
                                if re.match('^##INFO=', line_VCF):
                                    print(line_VCF, end='', file=fd_out)
                                if re.match('^##FILTER=', line_VCF):
                                    print(line_VCF, end='', file=fd_out)
                                if re.match('^##FORMAT=', line_VCF):
                                    print(line_VCF, end='', file=fd_out)
                                if re.match('^##contig=', line_VCF):
                                    print(line_VCF, end='', file=fd_out)
                                if re.match('^##reference=', line_VCF):
                                    print(line_VCF, end='', file=fd_out)
        ##                    elif not chrom_VCF and re.match('^##contig', line_VCF):
        ##                        assert chrom_VCF = re.match('^##contig=<ID=(\w+),', line_VCF).group(1)
                            ## Print INFO lines, so bcftools does not throw warnings.
        ##                    elif re.match('^##GATKCommandLine=<', line_VCF):
        ##                        pos_VCF = re.match(
        ##                            ''.join([
        ##                                '^##GATKCommandLine=',
        ##                                '<.+,CommandLineOptions=',
        ##                                '".+intervals=\[\w+:(\d+)-\d+\]']),
        ##                            line_VCF).group(1)
        ##                        ## assert that integer
        ##                        pos_VCF = str(int(pos_VCF))
                                pass
                            continue
                        assert line_VCF[:6] == '#CHROM'
                        if i == 0:
                            ## write sample IDs to output
                            print(line_VCF, end='', file=fd_out)
                        break
                    ## Loop over data lines.
                    for line_VCF in fd_source:
                        chrom_VCF, pos_VCF = line_VCF.split('\t', 2)[:2]
                        assert chrom == chrom_VCF
##                        print(chrom_VCF, pos_VCF, file=sys.stderr)
                        if pos_VCF == pos_INDEL:
                            assert chrom_VCF == chrom_INDEL
                            if VQSLod_INDEL >= d_minVQSLod['INDEL']:
                                print(line_VCF, end='', file=fd_out)
                            try:
                                (
                                    chrom_INDEL, pos_INDEL, VQSLod_INDEL
                                    ) = next(self.parse_recal(
                                        fd_recal_INDEL, pattern))
                            except StopIteration:
                                pass
        ##                        continue
        ##                    continue
                        else:
                            assert pos_VCF == pos_SNP
                            assert chrom_VCF == chrom_SNP
                            if VQSLod_SNP >= d_minVQSLod['SNP']:
                                print(line_VCF, end='', file=fd_out)
                            chrom_VCF, pos_VCF = line_VCF.split('\t', 2)[:2]
                            try:
                                chrom_SNP, pos_SNP, VQSLod_SNP = next(
                                    self.parse_recal(
                                        fd_recal_SNP, pattern))
                            except StopIteration:
                                pass
        ##                        continue
        ##                    continue
                        ## Continue loop over source input lines.
                        continue
                    ## Close source input file.
                    pass
                ## Continue loop over source input files.
                continue
            ## Close output and recal input files.
            pass

        ## index bgz output
        subprocess.call('{} -p vcf {}'.format(self.path_tabix, out), shell=True)
        ## confirm process has run to completion by writing to file
        with open('touch/ApplyRecalibration.touch', 'a') as f:
            f.write('{}\n'.format(out))
            f.write('{}.tbi\n'.format(out))

        ## return and continue with beagle if all AR processes completed
        return

    def beagle4(self):

        ## http://faculty.washington.edu/browning/beagle

        memMB = 16900  # todo: move to argparse
        memMB = 30900
        window = 50000  # todo: move to argparse
        if self.checkpoint:
            queue = 'normal'  # todo: move to argparse
            nthreads = 12
        else:
            queue = 'basement'
            nthreads = 1
            queue = 'long'
            nthreads = 12  # Beagle4 does not seem to scale well beyond 10.

        T_prev = 'ApplyRecalibration'

        ## write shell script
        if not os.path.isfile('shell/beagle.sh'):
            self.write_beagle_wrapper_script(memMB, nthreads, window)
        if self.checkpoint:
            self.write_brestart()

        ## Parse actual chromosome ranges after filtering.
        ## Execute shell script.
        d_pos_max = {}
        d_arguments = {'nthreads':nthreads}
        pattern = re.compile(r'.*VQSLOD=([-\d\.\w]+);')
        d_minVQSLod = self.parse_minVQSLods()
        for chrom in self.chroms:
            ## 1) Check input existence
            if self.check_in(
                T_prev, ['out_{}/{}.vcf.gz'.format(T_prev, chrom)],
                'touch/{}.touch'.format(T_prev)):
                continue
            ## 2) Check that process didn't start or end
            if self.touch('beagle.{}'.format(chrom)):
                continue
            ## Parse chromosome range.
            for mode in ('INDEL', 'SNP'):
                print(chrom, mode)
                tbx = pysam.TabixFile(
                    'out_VariantRecalibrator/{}.recal.gz'.format(mode))
                for line in tbx.fetch(chrom):
                    l = line.rstrip().split('\t')
                    VQSLod = re.match(pattern, l[7]).group(1)
                    if VQSLod == 'Infinity':
                        VQSLod = float('inf')
                    else:
                        VQSLod = float(VQSLod)
                    if VQSLod < d_minVQSLod[mode]:
                        continue
                    pos = int(l[1])
                    try:
                        d_pos_max[chrom] = max(int(pos), d_pos_max[chrom])
                    except KeyError:
                        d_pos_max[chrom] = int(pos)
                        print(mode, chrom)

            ## Execute shell script.
            print('beagle chrom', chrom)
            with gzip.open(
                'out_{}/{}.vcf.gz'.format(T_prev, chrom), 'rt') as fd_vcf:
                cnt = 0
                pos2 = 0
                ## Previous position in current fragment.
                pos_max = d_pos_max[chrom]
                for line in fd_vcf:
                    ## Skip metainformation lines and header line.
                    if line[:2] == '##':
                        continue
                    ## Break after reading header line.
                    break
                for line in fd_vcf:
                    l = line.split('\t', 2)
                    chrom = l[0]
                    pos = int(l[1])
                    ## Skip if position already covered by previous fragment.
                    if cnt == 0 and pos == pos2:
                        continue
                    cnt += 1
                    if cnt == 1:
                        pos1 = pos
                    if cnt % window == 0 or pos == pos_max:
                        pos2 = pos
                        index = cnt // window
                        if pos == pos_max and cnt % window > 0:
                            index += 1
                        d_arguments['out'] = 'out_beagle/{}/{}'.format(
                            chrom, index)
                        d_arguments['chrom'] = chrom
                        d_arguments['pos1'] = pos1
                        d_arguments['pos2'] = pos2
                        d_arguments['nthreads'] = nthreads
                        d_arguments['memMB'] = nthreads
                        arguments = self.args_dict2str(d_arguments)
                        ## Generate optional output with Beagle window ranges.
                        with open('lists/beagle.coords', 'a') as f:
                            f.write('{}\t{}\t{}\n'.format(chrom, pos1, pos2))
                        self.bsub_beagle(
                            chrom, pos1, pos2, index, memMB, queue, nthreads,
                            arguments=arguments)
                        print(chrom, ':', pos1, '-', pos2, index, cnt)
                        cnt = 0

##            pos2 = pos
##            index = (cnt // window) + 1
##
##            d_arguments['out'] = 'out_beagle/{}/{}'.format(chrom, index)
##            d_arguments['chrom'] = chrom
##            d_arguments['pos1'] = pos1
##            d_arguments['pos2'] = pos2
##            arguments = self.args_dict2str(d_arguments)
##
##            with open('lists/beagle.coords', 'a') as f:
##                f.write('{}\t{}\t{}\n'.format(chrom, pos1, pos2))
##            self.bsub_beagle(
##                LSF_memMB=memMB, LSF_queue=queue, LSF_n=nthreads,
##                variables=variables)

        return

    def write_beagle_wrapper_script(self, memMB, nthreads, window):

        ## initiate shell script
        lines = ['#!/bin/bash\n']
        ## Parse arguments from command line.
        lines += [self.args2getopts((
            'chrom', 'pos1', 'pos2', 'out', 'nthreads'))]
        ## Make parent directories of output.
        lines += ['mkdir -p $(dirname $out)']
        ## exit if output already exists
        lines += ['if [ -f $out.vcf.gz ]; then exit; fi']
        ## initiate beagle
        lines += ['{} \\'.format(self.init_java(self.jar_beagle, memMB))]
        ## Arguments for specifying data
        lines += [' gl=out_ApplyRecalibration/$chrom.vcf.gz \\']
        if self.ped:
            lines += [' ped={} \\'.format(self.ped)]
        lines += [' out=$out \\']
        if self.beagle4_excludesamples:
            lines += [' excludesamples={} \\'.format(
                self.beagle4_excludesamples)]
##        lines += [' excludemarkers={} \\'.format(excludemarkers)]
        lines += [' chrom=$chrom:$pos1-$pos2 \\']
        ## Other arguments
        lines += [' nthreads=$nthreads \\']
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
            'beagle', ['$out.vcf.gz'],
            extra='{} -p vcf $out.vcf.gz'.format(self.path_tabix))

        self.write_shell('shell/beagle.sh', lines,)

        return

    def bsub_beagle(
        self, chrom, pos1, pos2, index, memMB, queue, nthreads,
        arguments=''):

        ## Finished? Output tabix indexed.
        if os.path.isfile('out_beagle/{}/{}.vcf.gz.tbi'.format(chrom, index)):
            return

        ## Stalled? Log file exists.
        if os.path.isfile('out_beagle/{}/{}.log'.format(chrom, index)):
            ## Log file not updated for 24 hours.
            if time.time() - os.path.getmtime(
                'out_beagle/{}/{}.log'.format(chrom, index)) > 60*60*24:
##            if (time.time() - os.path.getmtime(
##                'out_beagle/{}/{}.log'.format(chrom, index)) > 60*60*24 or not
##                os.path.isfile('LSF/beagle/{}/{}.out'.format(chrom, index))):
                ## Output exists.
                if os.path.isfile(
                    'out_beagle/{}/{}.vcf.gz'.format(chrom, index)):
                    ## Output is empty.
                    if not os.path.getsize(
                        'out_beagle/{}/{}.vcf.gz'.format(chrom, index)):
                        print('Deleting vcf, log and LSF.out')
                        os.remove(
                            'out_beagle/{}/{}.vcf.gz'.format(chrom, index))
                        os.remove(
                            'out_beagle/{}/{}.log'.format(chrom, index))
                        if os.path.isfile(
                            'LSF/beagle/{}/{}.out'.format(chrom, index)):
                            os.remove(
                                'LSF/beagle/{}/{}.out'.format(chrom, index))

        ## started and running?
        fn = 'LSF/beagle/{}.{}.out'.format(chrom, index)
        if os.path.isfile(fn):
            print('a', fn)
            if os.path.getsize(fn):
                print('b', fn)
                ## running?
                with open(fn) as f:
                    line = f.readlines()[-1]
                    print(line)
                    ## Return if running.
                    if any([
                        'trios' in line,
                        'target markers' in line,
                        ## Starting burn-in iterations
                        ## Starting phasing iterations
                        'iterations' in line,
                        'mean edges' in line,  # mean edges/node
                        ]):
                        return
            else:
                stop

        ## started and finished? ## tmp!!!
        fn = 'out_beagle/{}/{}.log'.format(chrom, index)
        if os.path.isfile(fn):
            print('a', fn)
            if os.path.getsize(fn):
                print('b', fn)
                ## running?
                with open(fn) as f:
                    line = f.readlines()[-1]
                    print(line)
                    if line.rstrip().split()[-1] == 'finished':
                        stopshouldnothappen  # tmp!!!
##                        subprocess.call(
##                            'tabix -p vcf out_beagle/{}/{}.vcf.gz'.format(
##                                chrom, index), shell=True)
                        return
            else:
                stop

        print(chrom, index)

        LSB_JOBNAME = '{}.{}.{}'.format('beagle', chrom, index,)
        LSF_affix = '{}/{}/{}'.format('beagle', chrom, index)
        cmd_beagle = self.bsub_cmd(
            'beagle', LSB_JOBNAME, LSF_memMB=memMB, LSF_affix=LSF_affix,
            LSF_queue=queue, LSF_n=nthreads, arguments=arguments)

        if self.checkpoint:
            s = subprocess.check_output(cmd_beagle, shell=True).decode()
            print(s)
            jobID = int(re.match('.*?<(.*?)>', s).group(1))
            print(jobID)
            cmd_brestart = 'bsub -G {}'.format(self.project)
            cmd_brestart += ' -o brestart.out -e brestart.err'
            cmd_brestart += ' -q small -w "ended({})"'.format(jobID)
            cmd_brestart += ' bash shell/brestart.sh {:d} {} {:d}'.format(
                jobID, self.project, memMB)
            print(cmd_brestart)
            print()
            subprocess.call(cmd_brestart, shell=True)
        else:
            subprocess.call(cmd_beagle, shell=True)

        return

    def bsub_cmd(
        self,
        ## bsub command line arguments
        shell_affix, LSB_JOBNAME,
        LSF_queue='normal', LSF_memMB=4000, LSF_affix=None, LSF_n=1,
        ## bash command line arguments
        arguments='',
        ):

        if not LSF_affix:
            LSF_affix = shell_affix

        os.makedirs(
            os.path.dirname('LSF/{}.out'.format(LSF_affix)),
            exist_ok=True)

        cmd = 'bsub -J"{}" -q {}'.format(LSB_JOBNAME, LSF_queue)
        ## Do not add -G if run on hgs4.
        if all([
            any([
                'bc' in socket.gethostname(),
                'farm3' in socket.gethostname()]),
            'hgs4' not in socket.gethostname()]):
            cmd += ' -G {}'.format(self.project)
        cmd += " -M{:d} -R'select[mem>{:d}] rusage[mem={:d}]'".format(
            LSF_memMB, LSF_memMB, LSF_memMB)
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
##        cmd += ' {}'.format(variables)
        bash_cmd = ' bash {}/shell/{}.sh'.format(os.getcwd(), shell_affix)
        cmd += bash_cmd
        cmd += ' {}'.format(arguments)
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

        lines += self.init_GATK_cmd(T, memMB, tuple_args)
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

        lines += self.init_GATK_cmd(T, memMB, tuple_args)
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
        lines += [' -A StrandBiasBySample \\']
        lines += [' -A VariantType \\']
        lines += [' --standard_min_confidence_threshold_for_calling 30 \\']
        lines += [' --standard_min_confidence_threshold_for_emitting 30 \\']

        lines += ['"\n\neval $cmd']

        ## terminate shell script
        lines += self.term_cmd(T, ['$out.tbi'])

        ## write shell script
        self.write_shell('shell/{}.sh'.format(T), lines,)

        return

    def execmd(self, cmd):

        print(cmd)
        subprocess.call(cmd, shell=True)

        return

    def write_shell(self, fp, lines,):

        os.makedirs(os.path.dirname(fp), exist_ok=True)

        if type(lines) != list:
            print(type(lines))
            stop

        s = '\n'.join(lines) + '\n\n'
        fd = open(fp, 'w')
        fd.write(s)
        fd.close()
        os.system('chmod +x {}'.format(fp))

        return

    def alphanum_key(self, s):
        ## http://dave.st.germa.in/blog/2007/12/11/exception-handling-slow/
        NUM_RE = re.compile('([0-9]+)')
        return [int(c) if c.isdigit() else c for c in NUM_RE.split(s)]

    def sort_nicely(self, l):
        ## http://nedbatchelder.com/blog/200712/human_sorting.html
        """ Sort the given list in the way that humans expect.
        """
        l.sort(key=self.alphanum_key)
        return l

    def init_java(self, jar, memMB):

        s = '{} -Djava.io.tmpdir={}'.format(self.path_java, 'tmp')
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
        for dirname in ['']:
            d_l_fp_out[dirname] = []
            l = os.listdir(
                os.path.join(dirname, 'out_{}'.format(analysis_type)))
            for s in l:
                path1 = os.path.join('out_{}'.format(analysis_type), s)
                path2 = os.path.join(dirname, path1)
                ## append files in chromosomal subdirectories
                if os.path.isdir(path2):
                    l = os.listdir(path2)
                    for fn in l:
                        d_l_fp_out[dirname] += [os.path.join(path1, fn)]
                ## append files in main dir
                elif os.path.isfile(path2):
                    d_l_fp_out[dirname] += [path1]
                ## symlink to die
                elif os.path.islink(path2):
                    l = os.listdir(path2)
                    print(l)
                    stop
                    for fn in l:
                        d_l_fp_out[dirname] += [os.path.join(path1, fn)]                    
                else:
                    print(path2)
                    print(os.path.realpath(path2))
                    print(os.path.isfile(os.path.realpath(path2)))
                    stop_not_expected

        bool_exit = False
        for dirname, l_fp_out in d_l_fp_out.items():
            if len(set(l_fp_in) - set(l_fp_out)) > 0:
                print('{} and possibly {} other files not generated.'.format(
                    list(set(l_fp_in) - set(l_fp_out))[0],
                    len(set(l_fp_in) - set(l_fp_out)) - 1,))
                print('dirname', dirname)
                print(
                    '{} has not run to completion. Goodbye.'.format(
                        analysis_type))
                bool_exit = True
#                print(inspect.stack()[1])
##                sys.exit()

        return bool_exit

    def touch(self, analysis_type, delta=60*60*24*365):

        bool_return = False
        fn_touch = 'touch/{}.touch'.format(analysis_type)
        if os.path.isfile(fn_touch):
##            if time.time() - os.path.getmtime(fn_touch) < delta:
                if self.verbose == True:
                    print('in progress or completed:{}'.format(analysis_type))
                bool_return = True
        else:
            if not os.path.isdir(os.path.dirname(fn_touch)):
                os.mkdir(os.path.dirname(fn_touch))
            self.execmd('touch {}'.format(fn_touch))

        return bool_return

    def write_brestart(self,):

        ## clean up this ugly function!!!

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
            f.write(
                '''exit140=$(echo $bhist | grep TERM_RUNLIMIT | wc -l)\n''')
            ## exit code 16, pid taken
            f.write('''exit16=$(echo $bhist | sed 's/ *//g'| grep Exitedwithexitcode16 | wc -l)\n''')
            ## Checkpoint failed
            f.write('''cpfail=$(echo $bhist | sed 's/ *//g'|''')
            f.write(''' grep "Checkpointfailed" | wc -l)\n''')
            ## Done successfully
            f.write('''donesuc=$(echo $bhist | sed 's/ *//g'| grep "Donesuccessfully" | wc -l)\n''')
            ## exit if done succesfully
            f.write('if [ $donesuc -eq 1 ]; then echo $bhist >> bhist_success.tmp; exit; fi\n')
            ## exit if not checkpoint succeeded and not PID taken
            f.write('if [ $exit143 -eq 0 -a $cpsucc -eq 0 -a $exit13 -eq 0')
            f.write(' -a $exit16 -eq 0 ]; then echo $bhist')
            f.write(' >> bhist_unexpectederror.tmp; exit; fi\n')
            ## restart job and capture jobID
            f.write(
                's=$(brestart -G $project -M$memMB $pwd/checkpoint/$jobID)\n')
            f.write('''jobID=$(echo $s | awk -F "[<>]" '{print $2}')\n''')
            ## report if checkpoint failed
            f.write('if [ $cpfail -ne 0 ]; then echo $s')
            f.write(' >> checkpointfailed_brestartout.tmp; fi\n')
            ## be verbose
            f.write('echo s $s\n')
            f.write('echo jobID $jobID\n')
            f.write('echo memMB $memMB\n')
            ## bsub this chaperone restart script again
            f.write('bsub')
            f.write(" -R 'select[mem>'$memMB'] rusage[mem='$memMB']'")
            f.write(" -M$memMB \\\n")
            f.write(' -o brestart.out -e brestart.err \\\n')
            f.write(' -G $project -q normal -w "ended($jobID)" \\\n')
            f.write(' bash shell/brestart.sh $jobID $project $memMB\n')

        return

    def shell_UG(self, T, memMB, nct, nt):

        lines = ['#!/bin/bash\n']

##        ## exit if job started
##        lines += ['if [ -f $out ]; then exit; fi\n']
        ## exit if job finished
        lines += ['if [ -s $out.tbi ]; then exit; fi\n']

        ## initiate GATK command
        lines += self.init_GATK_cmd(T, memMB, (
            'out', 'chrom', 'pos1', 'pos2', 'input_file',
            'nct', 'nt', 'XL', 'sample_ploidy'))

        ## append GATK command options
        lines += self.body_UnifiedGenotyper()

        ## terminate shell script
        lines += self.term_cmd(T, ['$out.tbi'])

        ## write shell script
        self.write_shell('shell/{}.sh'.format(T), lines)

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
        s = '"\nif [ "$XL" != "" ]; then cmd=$cmd"'
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
        s_annotation += ' -A StrandBiasBySample'
        s_annotation += ' -A VariantType'
        lines += [' {} \\'.format(s_annotation)]

        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_genotyper_UnifiedGenotyper.html#--genotype_likelihoods_model
        lines += [' --genotype_likelihoods_model BOTH \\']  # default value SNP

        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_haplotypecaller_HaplotypeCaller.html#--genotyping_mode
        lines += [' -gt_mode {} \\'.format(self.genotyping_mode)]

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
            print(self.bams)
            print(self.project)
            stop
            lines += [' -stand_call_conf 4 \\']
            lines += [' -stand_emit_conf 4 \\']

        lines += ['"\n\neval $cmd']

        return lines

    def write_bam_and_XL_lists(self, d_chrom_ranges):

        os.makedirs('lists', exist_ok=True)

        ## Create region lists. One of them empty.
        ## http://gatkforums.broadinstitute.org/discussion/1204/what-input-files-does-the-gatk-accept-require
        with open('lists/XL.PAR.list', 'w') as f:
            for region in ('PAR1', 'PAR2'):
                f.write('X:{:d}-{:d}\n'.format(
                    d_chrom_ranges[region][0], d_chrom_ranges[region][1]))

        ## write bam lists
        with open('lists/bams.list', 'w') as f:
            for bam in self.bams:
                f.write('{}\n'.format(bam))

        ## Write gender bam lists.

        ## Convert samples to sex.
        d_sample2sex = {}
        for uri in self.sample_genders:
            print(urllib.parse.urlparse(uri).scheme)
            if urllib.parse.urlparse(uri).scheme == 'ftp':
                if not os.path.isfile(os.path.basename(uri)):
                    urllib.request.urlretrieve(
                        uri, filename=os.path.basename(uri))
                uri = os.path.basename(uri)
            with open(uri) as f:
                for line in f:
                    d_sample2sex[re.split(
                        '[, \t]', line)[0]] = re.split(
                            '[, \t]', line)[-1].lower()[0]

        ## Convert sex to list of bams.
        d_sex2bams = {'m': [], 'f': []}
        for bam in self.bams:
            sample = os.path.splitext(os.path.basename(bam))[0]
##                sample = pysam.Samfile(bam).header['RG'][0]['SM']  # slow
            sample = os.path.basename(bam).split('.')[0]  # fast
            try:
                sex = d_sample2sex[sample]
            except KeyError:
                print('sex missing', bam, sample)
                sys.exit()
            d_sex2bams[sex].append(bam)
        ## bam lists for males and females
        for sex in ('m', 'f'):
            with open('lists/bams.{}.list'.format(sex), 'w') as f:
                for bam in d_sex2bams[sex]:
                    f.write('{}\n'.format(bam))

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
        lines += [' -gt_mode {} \\'.format(self.genotyping_mode)]

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
            d_sex2bamlist = {'m': 'lists/bams.m.list'}
            XL = None
        elif chrom == 'X':
            d_sex2bamlist = {
                'm': 'lists/bams.m.list',
                'f': 'lists/bams.f.list'}
##                elif chrom in ('PAR1', 'PAR2'):
            XL = 'lists/XL.PAR.list'
        else:
            d_sex2bamlist = {'': 'lists/bams.list'}
            XL = None

        return d_sex2bamlist, XL

    def init_GATK_cmd(self, analysis_type, memMB, args):

        ## https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_engine_CommandLineGATK.php

        lines = []

        lines += [self.args2getopts(args)]

        lines += ['if [ "$nct" == "" ]; then nct=1; fi']

        ## exit if output exists
        lines += ['if [ -f $out ]; then exit; fi']

        ## create output folder
        lines += ['mkdir -p $(dirname $out)']

        ## Touch output to make other jobs exit instead of overwriting.
        lines += ['touch $out']

        s = ''
        s_java = self.init_java(self.jar_GATK, memMB)
        s += ' cmd="{} \\'.format(s_java)
        lines += ['\n{}'.format(s)]

        ## CommandLineGATK, required, in
        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_CommandLineGATK.html#--analysis_type
        lines += [' --analysis_type {} \\'.format(analysis_type)]
        ## CommandLineGATK, optional, in
        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_CommandLineGATK.html#--reference_sequence
        lines += [
            ' --reference_sequence {} \\'.format(self.reference_sequence)]
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
            lines += [
                'echo {} >> touch/{}.touch'.format(fp_out, analysis_type,)]

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

        s = self.args_rerun()
        fd = open('rerun_python.sh', 'w')
        fd.write(s)
        fd.close()
        self.execmd('chmod +x rerun_python.sh')

        return lines

    def args_rerun(self):

        s = ''
        s += ' {}'.format(sys.executable)
        s += ' {}'.format(sys.argv[0])
        s += self.args_dict2str(vars(self.args))

        return s

    def args_dict2str(self, _dict):

        s = ''
        for k, v in _dict.items():
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

    def is_file_or_ftp(self, str_):
        if not os.path.isfile(str_) and not os.path.islink(str_):
            if not re.match('^ftp://', str_):
                msg = '{} not found'.format(str_)
                raise argparse.ArgumentTypeError(msg)
        return str_

    def is_file_or_dir(self, str_):
        print(str_)
        if not any([
            os.path.isfile(str_), os.path.islink(str_), os.path.isdir(str_)]):
            msg = '{} is neither a readable file nor a directory' % str_
            raise argparse.ArgumentTypeError(msg)
        return str_

    def add_arguments(self, parser):

        ## required arguments

        parser.add_argument(
            '--path_bams', '--bam', '--input',
            help='Path to bam file and/or directory containing bam files',
            nargs='+', required=True, type=self.is_file_or_dir)

        parser.add_argument('--coverage', required=True, type=float)

        parser.add_argument(
            '--build', required=True, type=int, choices=[37, 38])

        parser.add_argument(
            '--sample_genders', '--sex', required=False, default=[],
            type=self.is_file_or_ftp, nargs='*')

        parser.add_argument(
            '--jar_GATK', '--path_GATK', '--GATK', '--gatk', required=True,
            help='File path to GATK', type=self.is_file)

        parser.add_argument('--project', required=True)

        parser.add_argument('--arguments', '--args')

        parser.add_argument(
            '--path_java', '--java', required=True, type=self.is_file)

        parser.add_argument(
            '--path_tabix', '--tabix', required=True, type=self.is_file)

        ##
        ## CommandLineGATK arguments
        ##

        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_CommandLineGATK.html#--intervals
        parser.add_argument('--intervals', '-L')

        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_CommandLineGATK.html#--reference_sequence
        parser.add_argument(
            '--reference_sequence', '-R', required=True, type=self.is_file)

        ##
        ## HaplotypeCaller specific arguments
        ##

        ## Optional Inputs

        ## https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_haplotypecaller_HaplotypeCaller.html#--alleles
        parser.add_argument('--alleles')

        ## https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_haplotypecaller_HaplotypeCaller.html#--dbsnp
        parser.add_argument('--dbsnp', '-D')

        ## Optional Parameters

        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_genotyper_UnifiedGenotyper.html#--genotyping_mode
        parser.add_argument(
            '--genotyping_mode', '-gt_mode', default='DISCOVERY')

        ##
        ## CombineGVCFs/GenotypeGVCFs related arguments
        ##

        ## http://gatkforums.broadinstitute.org/discussion/4074/file-number-limit-for-genotypegvcfs
        parser.add_argument('--gVCF_limit', default=200, type=int)

        ##
        ## VariantRecalibrator resources
        ##

        parser.add_argument(
            '--resources_SNP', '--VR_snp',
            help='Path to a file with -resource lines to append to GATK VR',)

        parser.add_argument(
            '--resources_INDEL', '--VR_indel',
            help='Path to a file with -resource lines to append to GATK VR',)

        parser.add_argument(
            '--an_SNP', nargs='+',
            choices=[
                'DP', 'QD', 'FS', 'SOR', 'MQ', 'MQRankSum', 'ReadPosRankSum',
                'InbreedingCoeff'],
            default='DP QD FS SOR MQ MQRankSum ReadPosRankSum'.split(),
                )

        parser.add_argument(
            '--an_indel', nargs='+',
            choices=[
                'DP', 'QD', 'FS', 'SOR', 'MQRankSum', 'ReadPosRankSum',
                'InbreedingCoeff'],
            default='DP QD FS SOR MQRankSum ReadPosRankSum'.split(),
                )

        ##
        ## ApplyRecalibration
        ##

        ## https://www.broadinstitute.org/gatk/guide/article?id=1259
        parser.add_argument(
            '--ts_SNP', type=float, required=True,)

        parser.add_argument(
            '--ts_INDEL', type=float, required=True,)

        parser.add_argument('--AR_input', nargs='+')

        ##
        ## beagle
        ##
        parser.add_argument(
            '--jar_beagle', '--beagle',
            help='File path to beagle.jar file (e.g. beagle_3.3.2.jar)',
            required=True,
            )

        parser.add_argument(
            '--beagle4_excludesamples', '--excludesamples',
            required=False, type=self.is_file)

        ##
        ## optional arguments
        ##
        parser.add_argument(
            '--checkpoint', action='store_true', default=False)

        parser.add_argument(
            '--chroms', type=str, nargs='+',
            default=[
                str(i + 1) for i in range(22)] + [
                    'X', 'Y', 'MT', 'PAR1', 'PAR2'])

        parser.add_argument(
            '--ped', type=self.is_file)

        parser.add_argument(
            '--path_bams_exclusion', '--bamXL', '--exclude_bams',
            required=False, type=self.is_file_or_dir, nargs='*')

        return parser

    def parse_bams(self, path_bams):

        list_bams = []
        for path_bam in path_bams:
            if os.path.isdir(path_bam):
                list_bams += glob.glob(os.path.join(path_bam, '*.bam'))
            elif os.path.isfile(path_bam):
                if os.path.splitext(path_bam)[1] == '.bam':
                    list_bams += [path_bam]
                else:
                    with open(path_bam) as f:
                        list_bams += f.read().rstrip().split('\n')
            else:
                print(path_bam, path_bams)
                stop_take_care_of_symlinks

        for bam in list_bams:
##            print('checking for existence of', bam)
            if not os.path.isfile(bam):
                print('bam does not exist', bam)
                sys.exit()
                
        return list_bams

    def parse_arguments(self):

        parser = argparse.ArgumentParser()

        parser = self.add_arguments(parser)

        ## parse arguments to argparse NameSpace
        self.args = namespace_args = parser.parse_args()

        ## setatrr
        for k, v in vars(namespace_args).items():
            setattr(self, k, v)

        if self.jar_GATK is None and self.options is None:
            parser.error('--GATK or --arguments')

        s_arguments = ''
        for k, v in vars(namespace_args).items():
            s_arguments += '{} {}\n'.format(k, v)

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

        self.bams = self.parse_bams(self.path_bams)
        if self.path_bams_exclusion:
            self.bams = list(set(self.bams) - set(
                self.parse_bams(self.path_bams_exclusion)))

        return

    def __init__(self):

        ## parse command line arguments
        self.parse_arguments()
        self.verbose = True

        return

if __name__ == '__main__':
    self = main()
    self.main()
