#!/bin/python

## T. Carstensen (tc9), M.S. Sandhu (ms23), D. Gurdasani (dg11)
## Wellcome Trust Sanger Institute, 2012

## http://www.broadinstitute.org/gatk/about#typical-workflows
## http://www.broadinstitute.org/gatk/guide/topic?name=best-practices

## built-ins
import math, os, sys, time
## New in version 2.7
import argparse

##
## todo20120724: tc9: split chromosomes into smaller parts before BEAGLE step
## to avoid requesting 16gb of memory as those nodes are limited
##

##
## todo20120816: tc9: automate split bgl file step
## between producebeagleinput and BEAGLE
##

##
## todo20120809: tc9/dg11: add a ReduceReads step before UnifiedGenotyper for speed purposes
## and for better variant calling? cf. slide 14 of https://www.dropbox.com/sh/e31kvbg5v63s51t/ajQmlTL6YH/ReduceReads.pdf
## "VQSR Filters are highly empowered by calling all samples together"
## "Reduced BAMs provides better results for large scale analysis projects ( > 100 samples) because it doesn't require batching."
## http://www.broadinstitute.org/gatk/guide/topic?name=best-practices
## "Even for single samples ReduceReads cuts the memory requirements, IO burden, and CPU costs of downstream tools significantly (10x or more) and so we recommend you preprocess analysis-ready BAM files with ReducedReads."
##

## todo20130204: tc9: make memory sample size dependent... only tested on 3 datasets with 100 samples each...

class main():

    def main(self):

        self.init()

        self.check_logs()

        ## parse chromsome lengths from reference sequence
        d_chrom_lens = self.parse_chrom_lens()

        ## define list of chromosomes
        l_chroms = [str(i) for i in xrange(1,22+1,)]+['X','Y',]

        ##
        ## write shell scripts
        ##
        self.main_UnifiedGenotyper(l_chroms,d_chrom_lens,)

        self.main_VariantRecalibrator(l_chroms,d_chrom_lens,)
        sys.exit()

        self.main_determine_TS_level()

        self.main_ApplyRecalibration()

        self.main_ProduceBeagleInput()

        self.BEAGLE()

        self.IMPUTE2(d_chrom_lens,)

        s = '#!/bin/bash\n'
        s += 'for chrom in {1..22} X Y\ndo\n'
        s += './shell/UnifiedGenotyper $CHROMOSOME\n'
        s += 'done'
        self.write_shell('shell/init.sh',s,)

        print self.generate_bsub_line(
            'UnifiedGenotyper',
            job_array_intervals='$maxindex',
            )
        stop

##        self.bsub_all_chroms_all_steps(d_array,d_shell,)

        print 'Execute shell/init.sh to start.'

        return


    def check_logs(self,):

        bool_error = False
        l_fn = os.listdir('stderr')
        for fn in l_fn:
            if os.path.getsize('stderr/%s' %(fn)) > 0:
                fd = open('stderr/%s' %(fn),'r')
                s = fd.read()
                fd.close()
                print s
                print 'stderr/%s' %(fn)
                bool_error = True
                break

        if bool_error == True:
            sys.exit(0)

        return


    def init(self,):

        self.mkdirs()

        return


    def determine_TS_level(self,):

        make_this_a_shell_script_instead

        ## http://www.broadinstitute.org/gsa/wiki/index.php/Variant_quality_score_recalibration
        ## "I took called variants until I found 99% of my known variable sites"

        lines = ['#!/bin/bash']

        ## 
        lines += ['python ~/github/ms23/GATK_pipeline/determine_TS_level.py \\']
        ## input
        lines += ['--input out_GATK/VariantRecalibrator.tranches \\']
        ## output
        lines += ['--output out_Tommy/determine_TS_level.txt \\']

        self.write_shell('shell/determine_TS_level.sh',lines,)

        return s


    def IMPUTE2(self,d_chrom_lens,):

##        intervals_IMPUTE2 = int(math.ceil(
##            d_chrom_lens[chrom]
##            /self.bps_per_interval_IMPUTE2
##            ))

        ## http://mathgen.stats.ox.ac.uk/impute/impute_v2.html
        ## /lustre/scratch107/projects/agv/imputation/code_FINAL/imputechunks.sh
        ## http://mathgen.stats.ox.ac.uk/impute/impute_v2_instructions.pdf
        ## http://mathgen.stats.ox.ac.uk/impute/example_one_phased_panel.html
        ## http://www.stats.ox.ac.uk/~marchini/software/gwas/file_format.html

        lines = []

        lines += ['#!/bin/bash']

        lines += ['CHROMOSOME=$1']
        lines += ['LENCHROMOSOME=$2']

        ##
        ## define files
        ##
        ## Input file options
        ## http://mathgen.stats.ox.ac.uk/impute/ALL_1000G_phase1integrated_v3_impute.tgz
        s_haps = self.fp_impute2_hap
        s_legend = self.fp_impute2_legend

        lines += ['if [ $CHROMOSOME=="X" ]\nthen']
        lines += ['map=%s' %(self.fp_impute2_map.replace('$CHROMOSOME','$CHROMOSOME_PAR1',),)]
        lines += ['else']
        lines += ['map=%s' %(self.fp_impute2_map)]
        lines += ['fi']
        
        ## http://mathgen.stats.ox.ac.uk/impute/input_file_options.html#-g
        fn_genotype    = 'out_BEAGLE/BeagleOutput.$CHROMOSOME.gen'
##        ## Strand alignment options
##        s_strand       = '/lustre/scratch107/projects/agv/imputation/trial3/mock_strand_file/strand_file_chr%s.txt' %(chrom)

        ## IMPUTE2 requires that you specify an analysis interval
        ## in order to prevent accidental whole-chromosome analyses.
        ## If you want to impute a region larger than 7 Mb
        ## (which is not generally recommended),
        ## you must activate the -allow_large_regions flag.
        lines += ['posmax=$((($LSB_JOBINDEX+0)*5000000))']
        lines += ['if test $posmax -ge $LENCHROMOSOME']
        lines += ['then posmax=$LENCHROMOSOME']
        lines += ['fi\n']
        s = '\n'.join(lines)

        s += '/software/hgi/impute2/v2.2.2/bin/impute2 \\'

        ##
        ## Required arguments
        ## http://mathgen.stats.ox.ac.uk/impute/required_arguments.html
        ##
        ## http://mathgen.stats.ox.ac.uk/impute/input_file_options.html#-g
        s += '-g %s \\' %(fn_genotype,)
        ## http://mathgen.stats.ox.ac.uk/impute/required_arguments.html#-m
        ## http://mathgen.stats.ox.ac.uk/impute/input_file_options.html#-m
        ## Ask Deepti whether a different recombination map file should be used
        s += '-m $map \\'
        ## http://mathgen.stats.ox.ac.uk/impute/basic_options.html#-int
        s += '-int $((($LSB_JOBINDEX-1)*%i+1)) $posmax \\' %(5000000)

        ##
        ## Basic options
        ## http://mathgen.stats.ox.ac.uk/impute/basic_options.html
        ##

        ##
        ## Output file options
        ## http://mathgen.stats.ox.ac.uk/impute/output_file_options.html
        ##
        s += '-o out_IMPUTE2/impute2.gen.$CHROMOSOME.$LSB_JOBINDEX \\'

        ##
        ## Input file options
        ##
        ## http://mathgen.stats.ox.ac.uk/impute/input_file_options.html#-l
        s += '-l %s \\' %(s_legend,)
        ## http://mathgen.stats.ox.ac.uk/impute/input_file_options.html#-h
        s += '-h %s \\' %(s_haps,)

        self.write_shell('shell/IMPUTE2.sh',s,)

        return s


    def write_shell(self,fp,lines,):

        if type(lines) != list:
            print type(lines)
            stop

        s = '\n'.join(lines)+'\n\n'
        fd = open(fp,'w')
        fd.write(s)
        fd.close()
        os.system('chmod +x %s' %(fp))

        return


    def BEAGLE(self,):

        '''this step took 3.4hrs on chromosome 22 of Uganda dataset'''

        ## http://faculty.washington.edu/browning/beagle/beagle_3.3.2_31Oct11.pdf

        ## http://www.broadinstitute.org/gsa/wiki/index.php/Interface_with_BEAGLE_imputation_software
        ## IMPORTANT: Due to BEAGLE memory restrictions,
        ## it's strongly recommended that BEAGLE be run on a separate chromosome-by-chromosome basis.
        ## In the current use case, BEAGLE uses RAM in a manner approximately proportional to the number of input markers.

        ## Then run beagle imputation (you may have to chunk up for imputation and use the "known" option with the 1000G reference)
        ## ask Deepti what the "known" option is...

        ## initiate command
        cmd = ''

        ## 
        cmd += 'CHROMOSOME=$1\n'

        ## split bgl
        cmd += 'if [ ! -s in_BEAGLE/BeagleInput.$CHROMOSOME.bgl.00 ]\nthen\n'
        cmd += 'python ~/github/ms23/conversion/BEAGLE_divide.py '
        cmd += '-i out_GATK/ProduceBeagleInput.$CHROMOSOME.excluded.bgl '
        cmd += '-m ../in_BEAGLE/ALL.chr$CHROMOSOME.phase1_release_v3.20101123.filt.renamed.markers '
        cmd += '-p ../in_BEAGLE/ALL.chr$CHROMOSOME.phase1_release_v3.20101123.filt.renamed.bgl '
        cmd += '-o in_BEAGLE/BeagleInput.$CHROMOSOME.bgl'
        cmd += '\nfi\n'

        ## do not loop over fragments, but submit job array instead

        cmd += 'if [ ! -s out_BEAGLE/BeagleOutput.$CHROMOSOME.bgl.00\nthen\n'
        
##        s += 'java -Xmx4000m -jar %s ' %(self.fp_software_beagle)
        cmd += 'java -Xmx16000m Djava.io.tmpdir=out_BEAGLE -jar %s ' %(self.fp_software_beagle)

##like=<unphased likelihood data file> where <unphased likelihood data file> is the name of 
##a genotype likelihoods file for unphased, unrelated data (see Section 2.2).   You may use 
##multiple like arguments if data from different cohorts are in different files.
        cmd += ' like=out_GATK/ProduceBeagleInput.$CHROMOSOME.bgl '
####arguments for phasing and imputing data ...
#### Arguments for specifying files
## phased=<phased unrelated file> where <phased unrelated file> is the name of a Beagle file
## containing phased unrelated data (see Section 2.1).
## You may use multiple phased arguments if data from different cohorts are in different files.
        cmd += ' phased=/lustre/scratch107/projects/uganda/users/tc9/in_BEAGLE/ALL.chr$CHROMOSOME.phase1_release_v3.20101123.filt.renamed.bgl '
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
        cmd += ' markers=/lustre/scratch107/projects/uganda/users/tc9/in_BEAGLE/ALL.chr$CHROMOSOME.phase1_release_v3.20101123.filt.renamed.markers '
####missing=<missing code> where <missing code> is the character or sequence of characters used to represent a missing allele (e.g. missing=-1 or missing=?).
#### The missing argument is required.
##        s += ' missing=? '

##nsamples=<number of samples> where <number of samples> is positive integer giving 
##the number of haplotype pairs to sample for each individual during each iteration of the 
##phasing algorithm. The nsamples argument is optional. The default value is nsamples=4.  
##If you are phasing an extremely large sample (say > 4000 individuals), you may want to 
##use a smaller nsamples parameter (e.g. 1 or 2) to reduce computation time.  If you are 
##phasing a small sample (say < 200 individuals), you may want to use a larger nsamples
##parameter (say 10 or 20) to increase accuracy.
        cmd += ' nsamples=20 '

## lowmem=<true/false> where <true/false> is true if a memory-efficient,
## but slower, implementation of the sampling algorithm should be used.
## If lowmem=true the running time will increase by a factor <= 2,
## and the memory usage will be essentially independent of the number of markers. The lowmem argument is optional. The default value is lowmem=false.
## For haplotype phase inference and imputation of missing data with default BEAGLE options, memory usage increases with the number of markers. It you need to reduce the amount of memory BEAGLE is using, try one or more of the following techniques:
## 1. Use the lowmem=true command line argument (see Section 3.2.2). The lowmem option makes memory requirements essentially independent of the number of markers.
        cmd += ' lowmem=true '

        ## non-optional output prefix
        cmd += ' out=out_BEAGLE/BeagleOutput.$CHROMOSOME.bgl'

        ## terminate BEAGLE
        cmd += '\n'

        cmd += 'fi\n'

        ##
        ## gunzip the output files after runs to completion
        ##
        cmd += '\n\n'
        for suffix in ['dose','phased','gprobs',]:
            fp = 'out_BEAGLE/BeagleOutput.$CHROMOSOME.bgl'
            fp += '.ProduceBeagleInput.$CHROMOSOME.bgl.%s.gz.$LSBJOBINDEX' %(
                suffix,
                )
            s_gunzip = 'gunzip %s\n' %(fp)
            cmd += s_gunzip
            if suffix == 'gprobs':
                fp_gprobs = fp

        ##
        ## IMPUTE2 takes gen files...
        ##
        ##
        ## gprobs > gen
        ##
        ## skip header
        s_more = 'more +2 %s' %(fp_gprobs,)
        ## print columns 1 and 2 when using field separator : (i.e. a replacement of field operator)
        s_awk1 = "awk -F : '{print $1, $2}'"
        ## print columns 1 and 2 with default field separator (i.e. append columns)
        s_awk2 = '''awk '{print $1":"$2, $1":"$2, $0}' '''
        ## cut out column 3 and do *not* print this column (i.e. remove "marker" from header and chromosomeID from subsequent rows)
        s_cut = 'cut -d " " -f 3 --complement > out_BEAGLE/BeagleOutput.$CHROMOSOME.gen'
        s_convert = ' | '.join([s_more,s_awk1,s_awk2,s_cut])
        cmd += '\n'+s_convert+'\n\n'
##        print s_convert
##        os.system(s_convert)
##        stop

        ##
        ## gprobs > sample
        ##
        fp_out = 'out_BEAGLE/chromosome$CHROMOSOME.sample'
        cmd += '''echo "ID_1 ID_2 missing\n0 0 0" > %s ; \
        head -1 %s \
        | sed 's/ /\\n/g' \
        | fgrep -v "marker" | fgrep -v "allele" \
        | uniq \
        | awk '{print $1, $1, "NA"}' \
        >> %s''' %(
            fp_out, fp_gprobs, fp_out,
            )
        cmd += '\n\n'

##        ## gen2ped (gen and sample to ped and map)
##        ## http://www.well.ox.ac.uk/~cfreeman/software/gwas/gtool.html
##        out_prefix = 'out_GTOOL/gtool.%s' %(chrom)
##        cmd += 'gtool=/nfs/team149/Software/usr/share/gtool/gtool'
##        cmd += '$gtool -G '
##        ## in
##        cmd += '--g out_BEAGLE/BeagleOutput.%s.gen ' %(chrom)
##        cmd += '--s out_BEAGLE/chromosome%s.sample ' %(chrom)
##        ## out
##        cmd += '--ped %s.ped ' %(out_prefix)
##        cmd += '--map %s.map ' %(out_prefix)
##        ## parameters
##        cmd += '--phenotype phenotype_1 --threshold 0.9 --snp'
##        cmd += '\n'

        self.write_shell('shell/BEAGLE.sh',cmd,)

        return cmd


    def bsub_all_chroms_all_steps(self,d_array,d_shell,):

        self.write_nonarray_scripts(d_shell,)

        self.write_array_scripts_chromosome(d_shell,)

        self.write_array_scripts_subchromosome(d_array,)

        if self.bool_sequential == True:
            self.write_initial_script()

        return


    def write_nonarray_scripts(self,d_shell,):

        ##
        ## write shell scripts for each non-array step
        ##
##        d_bsub_lines = {}
        for i_step in range(len(self.l_steps)):

            step = self.l_steps[i_step]

            if self.bool_sequential == True:

                ## no additional steps
                if i_step+1 == len(self.l_steps):
                    bsub_line_next = None
                ## this is not the final step; i.e. next step exists
                else:
                    step_next = self.l_steps[i_step+1]
                    if step_next == 'BEAGLE':
                        if step_next in self.d_memory.keys():
                            s_queue = self.d_queues[step_next]
                            memory_MB = self.d_memory[step_next]
                        else:
                            s_queue = 'normal'
                            memory_MB = 4000
                        bsub_line_next = self.generate_bsub_line(
                            step_next,
                            job_array_intervals = 24,
                            job_array_title = '%s' %(step_next[:6],),
                            bool_wait = False,
                            s_queue = s_queue,
                            memory_MB = memory_MB,
                            )
                    else:
                        if step_next == 'VariantRecalibrator':
                            fn_out = self.d_file_outputs[step_next][0].replace('.$CHROMOSOME','')
                            dn_out = self.d_dir_outputs[step_next]
                            fp_out = os.path.join(dn_out,fn_out,)
                            s_queue = 'normal'
                            if os.path.isfile(fp_out):
                                memory_MB = 4000 ## temporary ugly solution...
                            else:
                                memory_MB = 12000
                        elif step_next == 'ApplyRecalibration':
##                            ## 32GB not enough
##                            s_queue = 'long'
##                            memory_MB = 32000
                            s_queue = 'hugemem'
                            memory_MB = 64000
                        else:
                            s_queue = 'normal'
                            memory_MB = 4000
                        bsub_line_next = self.generate_bsub_line(
                            step_next,
                            s_queue = s_queue,
                            memory_MB = memory_MB,
                            )
##                    d_bsub_lines[step_next] = bsub_line_next
            else:
                bsub_line_next = None

            ## skip array steps
            if step in self.l_steps_array_subchromosome+self.l_steps_array_chromosome:
                continue

            s = self.check_input_output_shell(step)
            s += d_shell[step]
            s += '\n\nfi\n\n'
            self.write_shell_script(
                step,
                s,
                line_additional = bsub_line_next,
                )
##            print d_shell[step]
##            stop

        return


    def write_array_scripts_chromosome(self,d_shell,):

        self = GATK()

        chromosome_prefix = ''
        
        ##
        ## write shell script for UnifiedGenotyper (job array and first step)
        ##
        for step in self.l_steps_array_chromosome:
               
            if step == 'ApplyRecalibration':
                s = '#!/bin/bash\n'
            else:
                s = '''#!/bin/bash
LSB_JOBINDEX=$LSB_JOBINDEX
if [ ${LSB_JOBINDEX} -eq 23 ]
then
CHROMOSOME="X"
elif [ ${LSB_JOBINDEX} -eq 24 ]
then
CHROMOSOME="Y"
else
CHROMOSOME=$LSB_JOBINDEX
fi
echo $CHROMOSOME
'''

            ## do not execute if file output exists...
            s += '\n\n# do not execute if file output exists and if not previous file output exists\n'
            ## check that file exists and is not empty (-s)
            s_condition = self.check_input_output_shell(step)
            if step == 'ApplyRecalibration':
                s_condition = s_condition.replace('.$CHROMOSOME','')
            s += s_condition
            if step == 'ApplyRecalibration':
                s += self.ApplyRecalibration(chromosome_prefix)
            elif step == 'BEAGLE':
                s += self.BEAGLE('$CHROMOSOME')
            else:
                print step
                stop
            s += '\n\nfi\n'
            d_shell[step] = s

            ## append step following UnifiedGenotyper; e.g. CombineVariants
            if self.bool_sequential == True:
                step_next = self.l_steps[self.l_steps.index(step)+1]

                bsub_line_next = ''
                ## is LSB_JOBINDEX always equal to 1, if job is not submitted as an array?
                bool_wrap = False
                if not (
                    step in ['ApplyRecalibration',]
                    ):
                    bool_wrap = True
                if bool_wrap == True:
                    ## only call next script once...
                    bsub_line_next += '# only call next script once (e.g. when chromosome 1 finishes)\n'
                    bsub_line_next += '# but will still wait for all chromosomes associated with JOBID to finish\n'
                    bsub_line_next += 'if [ ${LSB_JOBINDEX} -eq 1 ]\n'
                    bsub_line_next += 'then\n'
                bsub_line_next += self.generate_bsub_line(
                    step_next,
                    s_queue = 'normal',
                    memory_MB = 4000,
                    )
                if bool_wrap == True:
                    bsub_line_next += 'fi\n'

            self.write_shell_script(
                step,
                d_shell[step],
                line_additional = bsub_line_next,
                )

        return


    def write_array_scripts_subchromosome(self,d_array,):
        
        ##
        ## write shell script for UnifiedGenotyper (job array and first step)
        ##
        for step in self.l_steps_array_subchromosome:

            s = '#!/bin/sh\n\n'
            for chrom in ['X','Y',]+range(1,22+1,):
##                if chromosome != 22: continue ## tmp
                chrom = str(chromosome)
##                ## replace if statements with dictionary...
##                if step == 'UnifiedGenotyper':
##                    fp = 'out_GATK/UnifiedGenotyper.%s.vcf.1' %(chrom)
##                elif step == 'IMPUTE2':
##                    fp = 'out_IMPUTE2/chromosome%s.impute2.1' %(chrom)
##                else:
##                    print step
##                    print self.d_file_outputs[step]
##                    stop

##                dn_out = self.d_dir_outputs[step]
##                fn_out = self.d_file_outputs[step][0].replace('$CHROMOSOME',chromosome).replace('.$LSB_JOBINDEX','')
##                fp_out = os.path.join(dn_out,fn_out,)
##                if (
##                    self.bool_skip_if_output_exists == True
##                    and
##                    os.path.isfile(fp_out)
##                    and
##                    os.path.getsize(fp_out) > 0
##                    ):
##                    continue

####                if self.l_steps.index(step) != 0:
####                    s += 'if [ ${LSB_JOBINDEX} -eq 1 ]\n'
####                    s += 'then\n'
##                if step != 'IMPUTE2': ## tmpxxx
##                    s += self.check_input_output_shell(step,bool_input=False,chromosome=chromosome,)
                s += self.check_input_output_shell(step,bool_input=False,chromosome=chromosome,)
                s += '%s\n' %(d_array[step][chrom])
####                if self.l_steps.index(step) != 0:
####                    s += 'fi\n'
                s += 'fi\n'

    ##        ## append step following UnifiedGenotyper; e.g. CombineVariants
    ##        if self.bool_sequential == True:
    ##            s_UG += d_bsub_lines[
    ##                self.l_steps[self.l_steps.index('UnifiedGenotyper')+1]
    ##                ]

            self.write_shell_script(step,s,)

        return


    def check_input_output_shell(self,step,bool_input=True,chromosome=None,):

        s = ''
        s += 'if [ '
        s_output = ' ! -s %s/%s ' %(
            self.d_dir_outputs[step],
            self.d_file_outputs[step][0].replace('$LSB_JOBINDEX','1'),
            )
        if (
            step in [
                'CombineVariants',
                'VariantRecalibrator',
                'determine_TS_level',
                'ProduceBeagleInput',
                ]
            ):
            s_output = s_output.replace('.$CHROMOSOME','')
        s += s_output
        ## not first step, so also check for input
        if self.l_steps.index(step) != 0:
            s_input = ' -a -s %s/%s ' %(
                self.d_dir_outputs[self.l_steps[self.l_steps.index(step)-1]],
                self.d_file_outputs[self.l_steps[self.l_steps.index(step)-1]][0].replace('$LSB_JOBINDEX','1'),
                )
            ## ApplyRecalibration (VQSR steps) produces *one* output, when run on all chromosomes simultaneously
            ## use this single output file to check for input for GATKs ProduceBeagleInput
            if (
                step in [
                    'VariantRecalibrator',
                    'ProduceBeagleInput',
                    'determine_TS_level',
                    ]
                ):
                s_input = s_input.replace('.$CHROMOSOME','')
            s += s_input
        s += ' ]\n'
##        if self.l_steps.index(step) == 0:
##            s += 'if [ ! -s %s/%s ]\n' %(
##                self.d_dir_outputs[step],
##                self.d_file_outputs[step][0].replace('$LSB_JOBINDEX','1'),
##                )
##        else:
##            ## if in and not out exists
##            s += 'if [ ! -s %s/%s -a -s %s/%s ]\n' %(
##                self.d_dir_outputs[step],
##                self.d_file_outputs[step][0].replace('$LSB_JOBINDEX','1'),
##                self.d_dir_outputs[self.l_steps[self.l_steps.index(step)-1]],
##                self.d_file_outputs[self.l_steps[self.l_steps.index(step)-1]][0].replace('$LSB_JOBINDEX','1'),
##                )
        s += 'then\n'

        if chromosome:
            s = s.replace('$CHROMOSOME',chromosome)

        return s


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


    def mkdirs(self):

        l_dn = [
            'stdout','stderr',
            'shell',
            'out_BEAGLE','out_Tommy','out_IMPUTE2',
            'out_UnifiedGenotyper',
            'out_VariantRecalibrator',
            ]

        ## create subdirs
        for dn in l_dn:
            if not os.path.isdir(dn):
                os.mkdir(dn)

        return


    def parse_arguments(self,):

        ## http://docs.python.org/2/library/argparse.html#module-argparse
        ## New in version 2.7
        ## optparse deprecated since version 2.7
        parser = argparse.ArgumentParser()

        ##
        ## add arguments
        ##

        parser.add_argument(
            '--bam','--bams','--bamdir',
            dest='fp_bams',
            help='Path to directory containing improved BAMs',
            metavar='FILE',default=None,
            required = True,
            )

        parser.add_argument(
            '--GATK',
            dest='fp_GATK',
            help='File path to GATK (e.g. /software/varinf/releases/GATK/GenomeAnalysisTK-1.4-15-gcd43f01/GenomeAnalysisTK.jar)',
            metavar='FILE',default=None,
            required = True,
            )

        parser.add_argument(
            '--FASTA','--reference','--reference-sequence','--reference_sequence',
            dest='fp_FASTA_reference_sequence',
            help='File path to reference sequence in FASTA format (e.g. /lustre/scratch111/resources/vrpipe/ref/Homo_sapiens/1000Genomes/human_g1k_v37.fasta)',
            metavar='FILE',default=None,
            required = True,
            )

        parser.add_argument(
            '--project','-P','-G',
            dest='project',
            help='Project',
            metavar='STRING',default=None,
            required = True,
            )

        parser.add_argument(
            '--arguments','--args','--options','--opts',
            dest='fp_arguments',
            metavar='FILE',default=None,
            required = False,
            )

        ##
        ## VariantRecalibrator resources
        ##

        parser.add_argument(
            '--resources','--VariantRecalibrator',
            dest='fp_resources',
            help='File path to a file with -resource lines to append to GATK VariantRecalibrator',
            metavar='FILE',default=None,
            required = False,
            )

##        parser.add_argument(
##            '--hapmap',
##            dest='fp_resource_hapmap',
##            help='File path to hapmap vcf to be used by VariantCalibrator (e.g. /lustre/scratch107/projects/uganda/users/tc9/in_GATK/hapmap_3.3.b37.sites.vcf)',
##            metavar='FILE',default=None,
##            required = False,
##            )
##
##        parser.add_argument(
##            '--omni',
##            dest='fp_resource_omni',
##            help='File path to omni vcf to be used by VariantCalibrator (e.g. /lustre/scratch107/projects/uganda/users/tc9/in_GATK/1000G_omni2.5.b37.sites.vcf)',
##            metavar='FILE',default=None,
##            required = False,
##            )

        ## dbSNP file. rsIDs from this file are used to populate the ID column of the output.
        s_help = 'File path to dbsnp vcf to be used by VariantCalibrator'
        s_help += ' (e.g. /lustre/scratch107/projects/uganda/users/tc9/in_GATK/dbsnp_135.b37.vcf)'
##        s_help += '\nYou can get the vcf file with this command:'
##        s_help += '\ncurl -u gsapubftp-anonymous: ftp.broadinstitute.org/bundle/1.5/b37/dbsnp_135.b37.vcf.gz -o dbsnp_135.b37.vcf.gz; gunzip dbsnp_135.b37.vcf.gz'
        parser.add_argument(
            '--dbsnp',
            dest='fp_vcf_dbsnp',
            help=s_help,
            metavar='FILE',default=None,
            required = True,
            )

        ##
        ## BEAGLE
        ##
        parser.add_argument(
            '--beagle','--BEAGLE','--BEAGLEjar',
            dest='fp_software_beagle',
            help='File path to BEAGLE .jar file (e.g. /nfs/team149/Software/usr/share/beagle_3.3.2.jar)',
            metavar='FILE',default=None,
            required = True,
            )

        ##
        ## IMPUTE2
        ##
        parser.add_argument(
            '--impute2','--IMPUTE2','--IMPUTE2jar',
            dest='fp_software_impute2',
            help='File path to BEAGLE .jar file (e.g. /nfs/team149/Software/usr/share/beagle_3.3.2.jar)',
            metavar='FILE',default=None,
            required = True,
            )

        parser.add_argument(
            '--hap','--impute2-hap',
            dest='fp_impute2_hap',
            help='File path to directory containing hap files used by IMPUTE2 (e.g. /nfs/t149_1kg/ALL_1000G_phase1integrated_v3_impute/ALL_1000G_phase1integrated_v3_chr$CHROMOSOME_impute.hap.gz)',
            ## Ask Deepti where this is downloaded from
            ## and whether it is always split by chromosome
            metavar='FILE',default=None,
            required = True,
            )

        parser.add_argument(
            '--legend','--impute2-legend',
            dest='fp_impute2_legend',
            help='File path to directory containing legend files used by IMPUTE2 (e.g. /nfs/t149_1kg/ALL_1000G_phase1integrated_v3_impute/ALL_1000G_phase1integrated_v3_chr$CHROMOSOME_impute.legend.gz)',
            ## Ask Deepti where this is downloaded from
            ## and whether it is always split by chromosome
            metavar='FILE',default=None,
            required = True,
            )

        parser.add_argument(
            '--map','--impute2-map',
            dest='fp_impute2_map',
            help='File path to directory containing map files used by IMPUTE2 (e.g. /nfs/t149_1kg/ALL_1000G_phase1integrated_v3_impute/genetic_map_chr$CHROMOSOME_combined_b37.txt)',
            ## Ask Deepti where this is downloaded from
            ## and whether it is always split by chromosome
            metavar='FILE',default=None,
            required = True,
            )

##        ## IMPUTE2 input files
##        self.fp_impute2_map = '/nfs/t149_1kg/ALL_1000G_phase1integrated_v3_impute/genetic_map_chr$CHROMOSOME_combined_b37.txt'
##        error_input(self.fp_impute2_map)

        ##
        ##
        ##
        
        ## http://docs.python.org/2/library/argparse.html#argparse.ArgumentParser.parse_args
        ## parse arguments to argparse NameSpace
        self.namespace_args = namespace_args = parser.parse_args()        

        ## http://docs.python.org/2/library/functions.html#vars
        for k,v in vars(namespace_args).items():
            setattr(self,k,v)
            continue


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
            if fp == None or fp == 'None': continue
            if fp == self.fp_bams:
                f = os.path.isdir
            else:
                f = os.path.isfile
            if not f(fp):
                print 'file path does not exist:', fp
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

        print 'dirname', os.path.dirname(self.fp_bams)
        print 'basename', os.path.basename(self.fp_bams)

        return


    def ProduceBeagleInput(self,):

        '''this walker takes approximately 5-10 minutes per chromosome to run'''

        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_beagle_ProduceBeagleInput.html
        ## "After variants are called and possibly filtered,
        ## the GATK walker ProduceBeagleInputWalker will take the resulting VCF as input,
        ## and will produce a likelihood file in BEAGLE format.

        lines = ['#!/bin/bash']

        lines += self.GATK_initiate('ProduceBeagleInput',)

        fp_out = 'out_ProduceBeagleInput/ProduceBeagleInput.bgl'

        ## GATKwalker, required, out
        lines += ['--out %s \\' %(fp_out,)]
        ## GATKwalker, required, in
        lines += [
            '--variant out_ApplyRecalibration/ApplyRecalibration.recalibrated.filtered.vcf \\'
            ]

        ## split bgl by chromosome
        cmd = 'cat %s ' %(fp_out)
        cmd += '''awk -v FS=":" 'NR>1'''
        cmd += '''{print>"out_GATK/ProduceBeagleInput."$1".bgl"}' '''
        lines += [cmd]

        s = '\n'.join(lines)+'\n\n'

        self.write_shell('shell/ProduceBeagleInput.sh',s,)

        return s


    def ApplyRecalibration(self,):

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

        fp_in = 'out_Tommy/determine_TS_level.txt'

        lines = ['ts_filter_level=($(cat %s))\n' %(fp_in)]

        lines += self.GATK_initiate('ApplyRecalibration',)

        fp_out = 'out_ApplyRecalibration/ApplyRecalibration.recalibrated.filtered'

        ## GATKwalker, required, in
        lines += ['--input out_CombineVariants/CombineVariants%s.vcf \\']
        ## GATKwalker, required, out
        lines += ['--out %s \\' %(fp_out,)]
        ## GATKwalker, required, in
        lines += ['--recal_file out_VariantRecalibrator/VariantRecalibrator%s.recal \\']
        lines += ['--tranches_file out_VariantRecalibrator/VariantRecalibrator.tranches \\']
        ## GATKwalker, optional, in
        ## ts_filter_level should correspond to targetTruthSensitivity prior to novelTiTv dropping (see tranches plot)
        ## default is 99.00
##        lines += ['--ts_filter_level 99.0 \\']
        lines += ['--ts_filter_level $ts_filter_level \\']

        s = '\n'.join(lines)+'\n\n'

        self.write_shell('shell/ApplyRecalibration.sh',s,)

        return s


    def main_VariantRecalibrator(self,l_chroms,d_chrom_lens):

        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_variantrecalibration_VariantRecalibrator.html

        ##
        ## check input existence
        ##
        fd = open('UnifiedGenotyper.touch','r')
        s = fd.read()
        fd.close()
        l_vcfs_out = s.strip().split('\n')
        
        l_vcfs_in = []
        for chrom in l_chroms:
            intervals = int(math.ceil(
                d_chrom_lens[chrom]/self.bps_per_interval))
            for interval in range(1,intervals+1,):
                fp_in = 'out_%s/%s.%s.%s.%i.vcf' %(
                    'UnifiedGenotyper',self.project,'UnifiedGenotyper',
                    chrom,interval,
                    )
                l_vcfs_in += [fp_in]

        ##
        ##
        ##
        T = analysis_type = 'VariantRecalibrator'
        fp_out = 'xxx'

        bool_return = self.touch(analysis_type)
##        if bool_return == True: return

        lines = self.GATK_initiate('VariantRecalibrator',)

        ## GATKwalker, required, in
        lines += ['--input %s \\' %(vcf) for vcf in l_vcfs_in]
        ## GATKwalker, required, out
        lines += ['--recal_file out_VariantRecalibrator/VariantRecalibrator.recal \\']
        lines += ['--tranches_file out_VariantRecalibrator/VariantRecalibrator.tranches \\']
        ## GATKwalker, required, in (ask Deepti about this...)
        lines += [
            '--use_annotation QD \\',
            '--use_annotation HaplotypeScore \\',
            '--use_annotation MQRankSum \\',
            '--use_annotation ReadPosRankSum \\',
            '--use_annotation MQ \\',
            '--use_annotation FS \\',
            '--use_annotation DP \\',
            ]
        ## GATKwalker, optional, in
        ## Ugandan QCed to be added..!
        fd = open('%s' %(self.fp_resources),'r')
        lines_resources = fd.readlines()
        fd.close()
        lines += ['%s \\' %(line.strip()) for line in lines_resources]

        l_TStranches = [100,]
        l_TStranches += [99+i/10. for i in range(9,0,-1,)]
        l_TStranches += [90+i/2. for i in range(18,-1,-1,)]
        s_TStranches = ''
        for TStranche in l_TStranches:
            s_TStranches += '--TStranche %.1f ' %(TStranche)
        lines += ['%s \\' %(s_TStranches)]

        lines += self.term_cmd(analysis_type,fp_out,)

        self.write_shell('shell/VariantRecalibrator.sh',lines,)

        J = 'VR'
##        std_suffix = '%s.%s.out' %(self.project,analysis_type,)
        cmd = self.bsub_cmd(analysis_type,J,memMB=13000,)
##        os.system(cmd)
        print cmd

        return


    def GATK_initiate(self,analysis_type,):

        s = 'cmd="'
        ## run GATK
        s += 'java '
        ## set maximum heap size
        if analysis_type == 'VariantRecalibrator':
            s += '-Xmx12g ' ## Max Memory :     11856 MB
        elif analysis_type == 'ApplyRecalibration':
            s += '-Xmx62g ' ## 28GB not enough
        else:
            s += '-Xmx4g '
##        s += '-Xmx4g '
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
##            if time.time()-os.path.getmtime(fn_touch) < 60*60:
                print 'in progress:', analysis_type
                bool_return = True
        os.system('touch %s' %(fn_touch))

        return bool_return


    def main_UnifiedGenotyper(self,l_chroms,d_chrom_lens,):

        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_genotyper_UnifiedGenotyper.html

        T = analysis_type = 'UnifiedGenotyper'

        bool_return = self.touch(analysis_type)
##        if bool_return == True: return

        fp_out = 'out_%s/%s.%s.$CHROMOSOME.$LSB_JOBINDEX.vcf' %(
            analysis_type,self.project,analysis_type,)

        ## initiate shell script
        lines = ['#!/bin/bash']

        ## commands prior to GATK command
        lines += self.init_UnifiedGenotyper(
            l_chroms,d_chrom_lens,)

        ## initiate GATK command
        lines += self.GATK_initiate('UnifiedGenotyper',)

        ## append GATK command options
        lines += self.UnifiedGenotyper(fp_out,)

        ## terminate shell script
        lines += self.term_cmd(analysis_type,fp_out)

        ## write shell script
        self.write_shell('shell/UnifiedGenotyper.sh',lines,)

        ## execute shell script
        for chrom in l_chroms:
            fp_stdout = 'stdout/%s.UnifiedGenotyper.%s.1.out' %(
                self.project,chrom,)
            if os.path.isfile(fp_stdout):
                continue
            fp_out_chrom = fp_out.replace('$LSB_JOBINDEX','1')
            fp_out_chrom = fp_out_chrom.replace('$CHROMOSOME',chrom)
            if os.path.isfile(fp_out_chrom):
                continue
            intervals = int(math.ceil(
                d_chrom_lens[chrom]/self.bps_per_interval))

            J = '%s%s[%i-%i]' %('UG',chrom,1,intervals,)
            std_suffix = '%s.%s.%s.%%I' %(self.project,analysis_type,chrom)
            cmd = self.bsub_cmd(
                analysis_type,J,std_suffix=std_suffix,memMB=2000,)
            os.system(cmd)

        return


    def bsub_cmd(
        self,
        analysis_type,
        J,
        queue='normal',memMB=4000,
        std_suffix=None,
        chrom=None,
        ):

        if not std_suffix:
            std_suffix = '%s.%s' %(self.project,analysis_type,)

        cmd = 'bsub -J"%s" -q %s' %(J,queue,)
        cmd += ' -G %s' %(self.project)
        cmd += " -M%i000 -R'select[mem>%i] rusage[mem=%i]'" %(
            memMB,memMB,memMB,)
        cmd += ' -o stdout/%s.out' %(std_suffix)
        cmd += ' -e stderr/%s.err' %(std_suffix)
        cmd += ' ./shell/%s.sh %s' %(analysis_type,chrom,)

        return cmd


    def UnifiedGenotyper(self,fp_out,):

        lines = []

        ##
        ## required
        ##
        ## File to which variants should be written
        lines += [' --out %s \\' %(fp_out)]

        ##
        ## append input files
        ##
        l = os.listdir(self.fp_bams)
        for s in l:
            if s[-4:] != '.bam':
                continue
            ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_CommandLineGATK.html#--input_file
            lines += [' --input_file %s/%s \\' %(
                self.fp_bams,s,)]

        ##
        ## CommandLineGATK, optional
        ##

        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_CommandLineGATK.html#--intervals
        lines += [
            '--intervals $CHROMOSOME:$((($LSB_JOBINDEX-1)*%i+1))-$posmax \\' %(
                int(self.bps_per_interval),)]

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
        lines += [' --output_mode EMIT_VARIANTS_ONLY \\'] ## default value EMIT_VARIANTS_ONLY

        ## http://www.broadinstitute.org/gatk/gatkdocs/#VariantAnnotatorannotations
        s_annotation = ''
        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_annotator_DepthOfCoverage.html
        s_annotation += ' --annotation DepthOfCoverage'
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


    def term_cmd(self,analysis_type,fp_out,):

        ## cont cmd
        lines = [';echo %s >> %s.touch;' %(fp_out,analysis_type,)]
        s = '/software/bin/python-2.7.3'
        s += ' %s/GATK_pipeline2.py' %(os.path.dirname(sys.argv[0]))
        for k,v in vars(self.namespace_args).items():
            s += ' --%s %s' %(k,v)
        lines += [s]
        ## term cmd
        lines += ['"']

        lines += ['echo $cmd']
        lines += ['eval $cmd']

        return lines


    def init_UnifiedGenotyper(self,l_chroms,d_chrom_lens,):

        lines = []

        lines += ['\n## parse chromosome from command line']
        lines += ['CHROMOSOME=$1']

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

        ##
        ## do not allow interval to exceed the length of the chromosome
        ## otherwise it will raise an I/O error (v. 1.4-15)
        ##
        lines += ['posmax=$((($LSB_JOBINDEX+0)*%i))' %(int(
            self.bps_per_interval
            ))]
        lines += ['if test $posmax -gt $LENCHROMOSOME']
        lines += ['then posmax=$LENCHROMOSOME']
        lines += ['fi\n']

        return lines


    def __init__(self,):

        ##
        ## parse command line arguments
        ##
        self.parse_arguments()

        ## A node on the cluster does approximately 4Mbp in 8 hours;
        ## i.e. UnifiedGenotyper, ugandan dataset, chromosome1=249Mbp=3.1weeks=504hours,
        ## so do intervals of 4Mbp to avoid cluster time out after 12 hours
        ## this should be set as a flag with optparse as it differs from dataset to dataset...
        ## this one varies from chromosome to chromosome...
        ## this variable should be dependent on the number of samples (and the coverage / number of variants... ???)
        ## todo: collect stats...
        ##
        ## the same variable is also used for IMPUTE2, so I decided to lower it from 10Mbp to 5Mbp
        self.bps_per_interval = 10.*10**6
        self.bps_per_interval_IMPUTE2 = 5.*10**6

        ## sequential order in which to run GATK walkers and other steps
        self.l_steps = [
            'UnifiedGenotyper',
            'CombineVariants',
            'VariantRecalibrator',
            'determine_TS_level', ## not a GATK step
            'ApplyRecalibration',
            ## Imputation, BEAGLE
            'ProduceBeagleInput',
            'BEAGLE', ## not a GATK step
            ## Imputation, IMPUTE2
            'IMPUTE2', ## not a GATK step
            ]

        ## which jobs need to be run as arrays because of special memory/processor requirements?
        self.l_steps_array_subchromosome = [
            'UnifiedGenotyper',
            'IMPUTE2',
            ]
        self.l_steps_array_chromosome = [
            'ApplyRecalibration',
            'BEAGLE',
            ]

        self.d_file_outputs = {
            'UnifiedGenotyper':['UnifiedGenotyper.$CHROMOSOME.vcf.$LSB_JOBINDEX',],
            'CombineVariants':['CombineVariants.$CHROMOSOME.vcf',],
            'VariantRecalibrator':[
                'VariantRecalibrator.$CHROMOSOME.tranches',
                'VariantRecalibrator.$CHROMOSOME.recal',
                ],
            'determine_TS_level':['determine_TS_level.$CHROMOSOME.txt',], ## not a GATK step
            'ApplyRecalibration':[
                'ApplyRecalibration.recalibrated.filtered.$CHROMOSOME.vcf',
                'ApplyRecalibration.recalibrated.filtered.$CHROMOSOME.vcf.idx',
                ],
            'ProduceBeagleInput':['ProduceBeagleInput.$CHROMOSOME.bgl',],
            'BEAGLE':[
                'BeagleOutput.$CHROMOSOME.bgl.ProduceBeagleInput.$CHROMOSOME.bgl.phased',
                'BeagleOutput.$CHROMOSOME.bgl.ProduceBeagleInput.$CHROMOSOME.bgl.r2',
                'BeagleOutput.$CHROMOSOME.bgl.ProduceBeagleInput.$CHROMOSOME.bgl.gprobs',
                'BeagleOutput.$CHROMOSOME.bgl.ProduceBeagleInput.$CHROMOSOME.bgl.dose',
                'BeagleOutput.$CHROMOSOME.bgl.ProduceBeagleInput.$CHROMOSOME.bgl.r2.idx',
                'BeagleOutput.$CHROMOSOME.bgl.ProduceBeagleInput.$CHROMOSOME.bgl.gprobs.idx',
                'BeagleOutput.$CHROMOSOME.bgl.ProduceBeagleInput.$CHROMOSOME.bgl.phased.idx',
                ],
            'IMPUTE2':[
                'chromosome$CHROMOSOME.impute2.gen.$LSB_JOBINDEX_summary',
                'chromosome$CHROMOSOME.impute2.gen.$LSB_JOBINDEX',
                ],
            }

        self.d_queues = {
            'BEAGLE':'long'
            }

        self.d_memory = {
            'BEAGLE':6000,
            'determine_TS_level':100,
            'CombineVariants':4000,
            'UnifiedGenotyper':2000,
            'VariantRecalibrator':12000,
            }

        ##
        ## binary options
        ##

        ## run everything in sequential;
        ## i.e. wait for previous steps to finish
        ## and then auto submit next step
        self.bool_sequential = True

        ## create empty shell script if file already exists
        ## should really skip the step entirely to avoid submitting 0 second job to cluster and get unpopular and acquire poor cluster stats on behalf of the ms23 group...
        self.bool_skip_if_output_exists = True

        ## not currently used, but will be implemented...
        self.bool_send_mail_upon_job_completion = True

        return


    def BeagleOutputToVCF(self,chrom,):

        ## this function is currently not called

        lines = self.GATK_initiate(analysis_type,)

        ## GATKwalker, required, in
        lines += [
            '--beaglePhased:BEAGLE out_BEAGLE/BeagleOutput.%s.bgl.ProduceBeagleInput.%s.bgl.phased \\' %(
                chromosome,chromosome,
                ),
            ]
        lines += [
            '--beagleProbs:BEAGLE out_BEAGLE/BeagleOutput.%s.bgl.ProduceBeagleInput.%s.bgl.gprobs \\' %(
                chromosome,chromosome,
                ),
            ]
        lines += [
            '--beagleR2:BEAGLE out_BEAGLE/BeagleOutput.%s.bgl.ProduceBeagleInput.%s.bgl.r2 \\' %(
                chromosome,chromosome,
                ),
            ]
        ## GATKwalker, required, out
        fp_out = 'out_GATK/%s.%s.vcf' %(analysis_type,chromosome,)
        lines += ['--out %s \\' %(fp_out,)]
        ## GATKwalker, required, in
        lines += [
            '--variant out_GATK/ApplyRecalibration.recalibrated.filtered.%s.vcf \\' %(
                chromosome
                )
            ]

        s = '\n'.join(lines)+'\n\n'

        return s


if __name__ == '__main__':
    self = main()
    self.main()
