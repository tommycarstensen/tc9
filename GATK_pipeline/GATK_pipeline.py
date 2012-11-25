#!/bin/python
# -*- coding: utf-8 -*-

## T. Carstensen (tc9), M.S. Sandhu (ms23), D. Gurdasani (dg11)
## Wellcome Trust Sanger Institute, 2012

## http://www.broadinstitute.org/gatk/about#typical-workflows

##
## todo: instead of running step by step, run chromosome by chromosome to ensure all steps are finished before proceeding
## no instead create a job array for each chromosome!!!
## multiple bsubs in one shell script will be submitted simultaneously... hmmm..
## if submitted to job small jobs will be run sequentially... too risky a way of queuing...
## see this great page referenced on the farm e-mail list:
## http://apps.sanger.ac.uk/ext-docs/lsf-7.0.6/admin/jobdependencies.html
##

##
## todo20120628: tc9: add some check that the output was generated before proceeding
##

##
## todo20120627: dg11: For WGS data, consider another refinement step with IMPUTE2
##

##
## todo20120628: tc9: add the option of reading parameters from parameter file and/or command line
## make sure that chosen parameters are written to a parameter.out file
##

##
## todo20120717: tc9: give option to run IMPUTE2 before Beagle and vice versa... need to change name of input/output files...
##

##
## todo20120717: tc9: option to clean up files from previous steps?
##

##
## todo20120724: tc9: instead of writing multiple shell scripts when doing sub chromosomal operations
## instead pass the size of each chromosome to a generalized shell script
##

##
## todo20120724: tc9: split chromosomes into smaller parts before BEAGLE step to avoid requesting 16gb of memory as those nodes are limited
##

##
## todo20120726: tc9: make sure that scripts that simple check for file inputs/outputs are submitted to the "small" queue to ensure immediate processing instead of sitting around pending in the fairshare queing system
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


## built-ins
import math, os, sys


class main():

    ## http://www.broadinstitute.org/gsa/gatkdocs/release/index.html

    def main(self):

        if self.bool_discontinue == True:
            return

        d_chromosome_lengths = self.parse_chromosome_ranges()

        d_array, d_shell = self.loop_chromosomes(d_chromosome_lengths,)

        self.bsub_all_chromosomes_all_steps(d_array,d_shell,)

        print 'Run bsub.sh to start.'
        print 'grep "Max Memory" stdout/*'
        print 'grep "CPU time" stdout/*'

        return


    def write_shell_script(
        self,
        prefix,s_out,
        line_additional = None,
        ):

        s = '#!/bin/sh\n\n'+s_out
        if line_additional:
            s += '\n\n'+line_additional
    
        fp = 'shell/%s.sh' %(prefix)
        fd = open(fp,'w')
        fd.write(s)
        fd.close()
        os.system('chmod +x %s' %(fp)) ## only works in unix environment...

        return


    def determine_TS_level(self,chromosome,):

        ## http://www.broadinstitute.org/gsa/wiki/index.php/Variant_quality_score_recalibration
        ## "I took called variants until I found 99% of my known variable sites"

        if chromosome == '':
            chromosome_suffix = ''
        else:
            chromosome_suffix = '.%s' %(chromosome)

        analysis_type = ''

        ## initiate lines
        lines = ['python ~/github/ms23/GATK_pipeline/determine_TS_level.py \\']
        ## input
        lines += ['--input out_GATK/VariantRecalibrator%s.tranches \\' %(chromosome_suffix,)]
        ## input
        lines += ['--output out_Tommy/determine_TS_level%s.txt \\' %(chromosome_suffix,)]

        s = '\n'.join(lines)+'\n\n'

        return s


    def loop_chromosomes(
        self,
        d_chromosome_lengths,
        ):

        instance_GATK = GATK()

        ## initiate dictionaries
        ## shell scripts that need to be individualized for each chromosome
        d_array = {}
        for step in self.l_steps_array_subchromosome:
            d_array[step] = {}
        ## shell scripts that do not need to have the $LSB_JOBINDEX parameter passed
        d_shell = {}
        for step in self.l_steps:
            if self.bool_join_chromosomes == True:
                if step in [
                    'CombineVariants',
                    'VariantRecalibrator',
                    ]:
                    d_shell[step] = '\n'.join(instance_GATK.GATK_initiate(step,))+'\n'
                else:
                    d_shell[step] = ''
            else:
                d_shell[step] = ''

        l_chromosomes = ['X','Y',]+range(1,22+1,)
        for chromosome in l_chromosomes:
##            if chromosome != 22: continue ## tmp
            chromosome = str(chromosome)
            intervals = int(math.ceil(d_chromosome_lengths[chromosome]/self.bps_per_interval))
            intervals_IMPUTE2 = int(math.ceil(d_chromosome_lengths[chromosome]/self.bps_per_interval_IMPUTE2))

            print 'looping over chromosome', chromosome, d_chromosome_lengths[chromosome]

            ##
            ## write shell scripts if array
            ## generate commands if not array; i.e. append to d_shell
            ##

            ##
            ## 1) Call variants (multisample variant calling)
            ##
            instance_GATK.UnifiedGenotyper(chromosome,d_chromosome_lengths,)

##        if self.bool_join_chromosomes == True:
##            l_chromosomes = ['']
##
##        for chromosome in l_chromosomes:

            ## Combine VCF records from different sources
            ## i.e. from different jobs of the job array
##            d_shell['CombineVariants'] += instance_GATK.GATK_CombineVariants(
            d_shell['CombineVariants'] += instance_GATK.CombineVariants(
                chromosome,d_chromosome_lengths,
                )

            ##
            ## 2) Filter input vcf file (variant quality score recalibration)
            ##

            if self.bool_join_chromosomes == False:

                ## Create a Gaussian mixture model
                ## by looking at the annotations values
                ## over a high quality subset of the input call set
                ## and then evaluate all input variants
                d_shell['VariantRecalibrator'] += instance_GATK.VariantRecalibrator(
                    chromosome,
                    )

                ## read out_GATK/VariantRecalibractor.tranches
                ## and decide on a Truth Sensitivity filter level
                ## i.e. when Ti/Tv ratio is 2.1
                d_shell['determine_TS_level'] += instance_main.determine_TS_level(
                    chromosome,
                    )
     
                ## Applies cuts to the input vcf file (by adding filter lines)
                ## to achieve the desired novel truth sensitivity levels
                ## which were specified during VariantRecalibration
                d_shell['ApplyRecalibration'] += instance_GATK.ApplyRecalibration(
                    chromosome,
                    )

##            if self.bool_join_chromosomes == False:
##                d_shell['PrintReads'] += instance_GATK.PrintReads(
##                    chromosome,
##                    )

            ##
            ## 3) Imputation / Genotype refinement
            ##

            ## Convert the input VCF into a format
            ## accepted by the Beagle imputation/analysis program.
            d_shell['ProduceBeagleInput'] += instance_GATK.ProduceBeagleInput(
                chromosome,
                )

##            ## vcf (GATK) --> ped --> gen (IMPUTE2)
##            d_shell['ProduceBeagleInput'] += instance_GATK.ProduceBeagleInput(
##                chromosome,
##                )

            ## Imputation / Genotype refinement
            d_shell['BEAGLE'] += self.BEAGLE(chromosome,)

            ## Take files produced by Beagle imputation engine
            ## and create a vcf with modified annotations.
            if 'BeagleOutputToVCF' in self.l_steps:
                d_shell['BeagleOutputToVCF'] += instance_GATK.BeagleOutputToVCF(
                    chromosome,
                    )

            ## IMPUTE2
            self.IMPUTE2(chromosome,d_chromosome_lengths,)

            ##
            ## end of GATK steps
            ##
            
            ##
            ## append steps that are run as job arrays because of memory/processor requirements
            ## split into sub chromosomal chunks
            ##
            for step in self.l_steps_array_subchromosome:
                if step == 'IMPUTE2':
                    job_array_intervals = intervals_IMPUTE2
                else:
                    job_array_intervals = intervals
                d_array[step][chromosome] = self.generate_bsub_line(
                    '%s.%s' %(step,chromosome),
                    job_array_intervals = job_array_intervals,
                    job_array_title = '%s(%i)' %(chromosome,job_array_intervals,),
                    bool_wait = False,
                    )

##        ##
##        ## append steps that are run as job arrays because of memory/processor requirements
##        ## split into chromosomes
##        ##
##        for step in self.l_steps_array_chromosome:
##            d_array[step] = self.generate_bsub_line(
##                '%s' %(step,),
##                job_array_intervals = 22,
##                job_array_title = '',
##                bool_wait = False,
##                )
##            print d_array[step]
##            print step
##            stop
##            stop

        ##
        ## temporary solution for running in parallel
        ## this is so poorly written and makes so many stupid assumptions
        ## I need to rewrite this... the code is REALLY ugly!!!
        ##
        if self.bool_join_chromosomes == True:
            for step in self.l_steps:
##                ## A VCF file is not always the output...
##                fp_out = 'out_GATK/%s.vcf' %(step,)
                dn_out = instance_main.d_dir_outputs[step]
                fn_out = instance_main.d_file_outputs[step][0].replace('.$CHROMOSOME','').replace('.$LSB_JOBINDEX','')
                fp_out = os.path.join(dn_out,fn_out,)
##                if (
##                    self.bool_skip_if_output_exists == True
##                    and
##                    os.path.isfile(fp_out)
##                    and
##                    os.path.getsize(fp_out) > 0
##                    ):
##                    d_shell[step] = '' 
##                    continue
                if step == 'UnifiedGenotyper':
                    continue
                elif step in ['CombineVariants',]:
                    ## convert string to list
                    lines = d_shell[step].split('\n')
                    ## remove double line breaks
                    lines = [line for line in lines if line != '']
                    ## add output file, which wasn't done previously
                    lines += [' --out %s \\' %(fp_out)]
                    ## convert lines back to string and check for output (which should be done in a shell script instead btw...)
                    s = instance_GATK.GATK_terminate(lines,fp_out,)
                    d_shell[step] = s
                ##
                ## VQSR on whole chromosome
                ##
                elif step in ['VariantRecalibrator',]:
                    s = instance_GATK.VariantRecalibrator('')
                    d_shell[step] = s
                elif step in ['determine_TS_level',]:
                    s = instance_main.determine_TS_level('')
                    d_shell[step] = s
                elif step == 'ApplyRecalibration':
                    s = instance_GATK.VariantRecalibrator('')
                    d_shell[step] = s
                elif step == 'ProduceBeagleInput':
                    s = instance_GATK.ProduceBeagleInput('')
                    d_shell[step] = s
##                else:
##                    print d_shell[step]
##                    print step
##                    stop
                    
        return d_array, d_shell


##    def PrintReads(self,chromosome,):
##
##        analysis_type = 'PrintReads'
##
##        lines += self.GATK_initiate(analysis_type,)
##
####        fp_out = 'out_GATK/ApplyRecalibration.recalibrated.filtered.%s.vcf' %(chromosome)
##        dn_out = instance_main.d_dir_outputs[analysis_type]
##        fn_out = instance_main.d_file_outputs[analysis_type][0]
##        fp_out = os.path.join(dn_out,fn_out,)
##
##        ## GATKwalker, required, in
##        lines += ['--intervals %s \\' %(chromosome,)]
##        ## GATKwalker, required, out
##        lines += ['--out %s \\' %(fp_out,)]
##        ## GATKwalker, required, in
##        lines += ['--recal_file out_GATK/VariantRecalibrator%s.recal \\' %(chromosome_suffix)]
##        lines += ['--tranches_file out_GATK/VariantRecalibrator%s.tranches \\' %(chromosome_suffix)]
##        ## GATKwalker, optional, in
##        ## ts_filter_level should correspond to targetTruthSensitivity prior to novelTiTv dropping (see tranches plot)
##        ## default is 99.00
####        lines += ['--ts_filter_level 99.0 \\']
##        lines += ['--ts_filter_level $ts_filter_level \\']
##
##        s = self.GATK_terminate(lines,fp_out,)
##
##        return s


    def IMPUTE2(self,chromosome,d_chromosome_lengths,):

        ## http://mathgen.stats.ox.ac.uk/impute/impute_v2.html
        ## /lustre/scratch107/projects/agv/imputation/code_FINAL/imputechunks.sh
        ## http://mathgen.stats.ox.ac.uk/impute/impute_v2_instructions.pdf
        ## http://mathgen.stats.ox.ac.uk/impute/example_one_phased_panel.html
        ## http://www.stats.ox.ac.uk/~marchini/software/gwas/file_format.html

        ##
        ## define files (should be in init function...)
        ##
        ## Input file options
##        s_haps         = '/lustre/scratch107/resources/1000g/release/20120117/impute/chr%s.hap' %(chromosome,)
##        s_legend       = '/lustre/scratch107/resources/1000g/release/20120117/impute/chr%s.legend' %(chromosome)
        ## http://mathgen.stats.ox.ac.uk/impute/ALL_1000G_phase1integrated_v3_impute.tgz
        s_haps = self.hap.replace('$CHROMOSOME',chromosome,)
        s_legend = self.legend.replace('$CHROMOSOME',chromosome,)
##        ## Ask Deepti about these input file options...
##        s_known_haps_g = '/lustre/scratch107/projects/agv/shapeit/trial2/prephasing_output/Muganda_autosomes_QCed_excambigsnps_chr%s.haps' %(chromosome)
##        s_map          = '/lustre/scratch107/projects/agv/shapeit/genetic_maps_b37/genetic_map_chr%s_combined_b37.txt' %(chromosome)
        if chromosome == 'X':
            ## ask Deepti if this is the right map file to use for X chromosome
            ## the RSIDs seem to correspond to those in the out_BEAGLE/BeagleOutput.X.gen file
            ## alternatively /lustre/scratch111/resources/variation/grch37/ALL_1000G_phase1integrated_v3_impute
            s_map      = '/nfs/t149_1kg/ALL_1000G_phase1integrated_v3_impute/genetic_map_chr%s_combined_b37.txt' %(chromosome+'_PAR1')
        else:
            s_map      = '/nfs/t149_1kg/ALL_1000G_phase1integrated_v3_impute/genetic_map_chr%s_combined_b37.txt' %(chromosome)
        ## http://mathgen.stats.ox.ac.uk/impute/input_file_options.html#-g
        fn_genotype    = 'out_BEAGLE/BeagleOutput.%s.gen' %(
            chromosome,
            )
##        ## Strand alignment options
##        s_strand       = '/lustre/scratch107/projects/agv/imputation/trial3/mock_strand_file/strand_file_chr%s.txt' %(chromosome)

        ## IMPUTE2 requires that you specify an analysis interval
        ## in order to prevent accidental whole-chromosome analyses.
        ## If you want to impute a region larger than 7 Mb
        ## (which is not generally recommended),
        ## you must activate the -allow_large_regions flag.
        lines = []
        lines += ['max=$((($LSB_JOBINDEX+0)*5000000))']
        lines += ['if test $max -ge %i' %(d_chromosome_lengths[chromosome])]
        lines += ['then max=%i' %(d_chromosome_lengths[chromosome])]
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
        s += '-m %s \\' %(s_map)
        ## http://mathgen.stats.ox.ac.uk/impute/basic_options.html#-int
        s += '-int $((($LSB_JOBINDEX-1)*%i+1)) $max \\' %(5000000)

        ##
        ## Basic options
        ## http://mathgen.stats.ox.ac.uk/impute/basic_options.html
        ##

        ##
        ## Output file options
        ## http://mathgen.stats.ox.ac.uk/impute/output_file_options.html
        ##
        s += '-o out_IMPUTE2/chromosome%s.impute2.gen.$LSB_JOBINDEX \\' %(
            chromosome,
            )

        ##
        ## Input file options
        ##
        ## http://mathgen.stats.ox.ac.uk/impute/input_file_options.html#-l
        s += '-l %s \\' %(s_legend,)
        ## http://mathgen.stats.ox.ac.uk/impute/input_file_options.html#-h
        s += '-h %s \\' %(s_haps,)

        ##
        ## write shell script
        ##
##        if os.path.isfile('out_IMPUTE2/chromosome%s.impute2.$LSB_JOBINDEX'):
##            s = ''

        ## requires individual shell scripts,
        ## because job arrays are used for each chromosome
        fp = 'shell/IMPUTE2.%s.sh' %(chromosome)
        fd = open(fp,'w')
        fd.write(s)
        fd.close()
        os.system('chmod +x %s' %(fp)) ## only works in unix environment...

        return s


    def BEAGLE(self,chromosome,):

        '''this step took 3.4hrs on chromosome 22 of Uganda dataset'''

        ## http://faculty.washington.edu/browning/beagle/beagle_3.3.2_31Oct11.pdf

        ## http://www.broadinstitute.org/gsa/wiki/index.php/Interface_with_BEAGLE_imputation_software
        ## IMPORTANT: Due to BEAGLE memory restrictions,
        ## it's strongly recommended that BEAGLE be run on a separate chromosome-by-chromosome basis.
        ## In the current use case, BEAGLE uses RAM in a manner approximately proportional to the number of input markers.

        ## Then run beagle imputation (you may have to chunk up for imputation and use the "known" option with the 1000G reference)
        ## ask Deepti what the "known" option is...

##        if chromosome != '':
##            chromosome_suffix = '.%s' %(chromosome)
##        else:
##            chromosome_suffix = ''

        ## chromosome 1 requires more than 8GB, will need to request a machine with sufficient memory as well... to be fixed...
##        if chromosome in [1,2,]:
##            s = 'java -Xmx16000m -jar /nfs/team149/Software/usr/share/beagle_3.3.2.jar '
##        else:
##            s = 'java -Xmx4000m -jar /nfs/team149/Software/usr/share/beagle_3.3.2.jar '
        s = ''

##        s += 'if [ ! -s %s/%s ' %(
##        s += self.check_input_output_shell('BEAGLE').replace('$CHROMOSOME',chromosome,)
        
##        s += 'java -Xmx4000m -jar /nfs/team149/Software/usr/share/beagle_3.3.2.jar '
        s += 'java -Xmx16000m Djava.io.tmpdir=out_BEAGLE -jar /nfs/team149/Software/usr/share/beagle_3.3.2.jar '

##like=<unphased likelihood data file> where <unphased likelihood data file> is the name of 
##a genotype likelihoods file for unphased, unrelated data (see Section 2.2).   You may use 
##multiple like arguments if data from different cohorts are in different files.
        s += ' like=out_GATK/ProduceBeagleInput.%s.bgl ' %(chromosome)
        ## ask Deepti why "like" (input file is genotype likelihoods)
        ## 04sep2012 - Deepti noticed the missing phased argument :(
####arguments for phasing and imputing data ...
#### Arguments for specifying files
## phased=<phased unrelated file> where <phased unrelated file> is the name of a Beagle file
## containing phased unrelated data (see Section 2.1).
## You may use multiple phased arguments if data from different cohorts are in different files.
        s += ' phased=/lustre/scratch107/projects/uganda/users/tc9/in_BEAGLE/ALL.chr%s.phase1_release_v3.20101123.filt.renamed.bgl ' %(chromosome)
####  unphased=<unphased data file>                     (optional)
####  phased=<phased data file>                         (optional)
####  pairs=<unphased pair data file>                   (optional)
####  trios=<unphased trio data file>                   (optional)
####  like=<unphased likelihood data file>              (optional)
##markers=<markers file> where <markers file> is the name of the markers file containing
## marker identifiers, positions, and alleles described in Section 2.4.
## The markers argument is optional if you specify only one Beagle file,
## and is required if you specify more than one Beagle file.
##        s += ' markers=/lustre/scratch107/projects/uganda/users/tc9/in_BEAGLE/ALL.chr%s.phase1_release_v3.20101123.filt.markers ' %(chromosome)
        s += ' markers=/lustre/scratch107/projects/uganda/users/tc9/in_BEAGLE/ALL.chr%s.phase1_release_v3.20101123.filt.renamed.markers ' %(chromosome)
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
        s += ' nsamples=20 '

## lowmem=<true/false> where <true/false> is true if a memory-efficient, but slower, implementation of the sampling algorithm should be used. If lowmem=true the running time will increase by a factor ≤ 2, and the memory usage will be essentially independent of the number of markers. The lowmem argument is optional. The default value is lowmem=false.
## For haplotype phase inference and imputation of missing data with default BEAGLE options, memory usage increases with the number of markers. It you need to reduce the amount of memory BEAGLE is using, try one or more of the following techniques:
## 1. Use the lowmem=true command line argument (see Section 3.2.2). The lowmem option makes memory requirements essentially independent of the number of markers.
        s += ' lowmem=true '

        ## optional output prefix --- NOT OPTIONAL!!! took me bloody 15 minutes to debug :(
        s += ' out=out_BEAGLE/BeagleOutput.%s.bgl' %(chromosome)

##        s += '\nfi\n'

        ##
        ## gunzip the output files after runs to completion
        ##
        s += '\n\n'
        for suffix in ['dose','gprobs','phased',]:
            fp = 'out_BEAGLE/BeagleOutput.%s.bgl.ProduceBeagleInput.%s.bgl.%s.gz' %(
                chromosome,chromosome,suffix,
                )
##            ## should check that file exists to avoid error...
##            if not os.path.isfile(fp):
##                continue
            s_gunzip = 'gunzip %s\n' %(fp)
            s += s_gunzip

        ##
        ## IMPUTE2 takes gen files...
        ## Deepti one-liner for file conversion
        ## instead write a script to avoid looping over all lines 4 times
        ## on second thought: not that time consuming, so I can't be bothered
        ## two consecutive greps however is a waste... replace with "text1|text2"...
        ##
        ## gprobs > gen
        ## --revert-match (displays all lines *not* matching a pattern)
##        s_fgrep = 'fgrep -v STDY BeagleOutput.%s.bgl.ProduceBeagleInput.%s.bgl.gprobs' %(chromosome,chromosome,)
        s_more = 'more +2 out_BEAGLE/BeagleOutput.%s.bgl.ProduceBeagleInput.%s.bgl.gprobs' %(chromosome,chromosome,)
        ## print columns 1 and 2 when using field separator : (i.e. a replacement of field operator)
        s_awk1 = "awk -F : '{print $1, $2}'"
        ## print columns 1 and 2 with default field separator (i.e. append columns)
        s_awk2 = '''awk '{print $1":"$2, $1":"$2, $0}' '''
        ## cut out column 3 and do *not* print this column (i.e. remove "marker" from header and chromosomeID from subsequent rows)
        s_cut = 'cut -d " " -f 3 --complement > out_BEAGLE/BeagleOutput.%s.gen' %(chromosome)
        s_convert = ' | '.join([s_more,s_awk1,s_awk2,s_cut])
        s += '\n'+s_convert+'\n\n'
##        print s_convert
##        os.system(s_convert)
##        stop

        ## gprobs > sample
        fp_out = 'out_BEAGLE/chromosome%s.sample' %(chromosome)
        s += '''echo "ID_1 ID_2 missing\n0 0 0" > %s ; \
        head -1 out_BEAGLE/BeagleOutput.%s.bgl.ProduceBeagleInput.%s.bgl.gprobs \
        | sed 's/ /\\n/g' \
        | fgrep -v "marker" | fgrep -v "allele" \
        | uniq \
        | awk '{print $1, $1, "NA"}' \
        >> %s''' %(
            fp_out, chromosome, chromosome, fp_out,
            )
        s += '\n\n'

        ## gen2ped (gen and sample to ped and map)
        ## http://www.well.ox.ac.uk/~cfreeman/software/gwas/gtool.html
        out_prefix = 'out_GTOOL/gtool.%s' %(chromosome)
        s += 'gtool=/nfs/team149/Software/usr/share/gtool/gtool'
        s += '$gtool -G '
        ## in
        s += '--g out_BEAGLE/BeagleOutput.%s.gen ' %(chromosome)
        s += '--s out_BEAGLE/chromosome%s.sample ' %(chromosome)
        ## out
        s += '--ped %s.ped ' %(out_prefix)
        s += '--map %s.map ' %(out_prefix)
        ## parameters
        s += '--phenotype phenotype_1 --threshold 0.9 --snp'
        s += '\n'

        ## plink

        return s


    def generate_bsub_line(
        self,
        prefix_out,
        prefix_shell = None, ## also boolean for write to file or not
        ## job array options
        job_array_intervals = None,
        job_array_title = None,
        ## other options
        s_queue = 'normal',
        memory_MB = 4000,
        ## wait for previous job if run in sequence?
        bool_wait = True,
        ):

        ## initiate line
        s = 'bsub'

        ##  job name, initiate
        if job_array_title:
            s += ' -J"%s' %(job_array_title)
        else:
            ## first ten characters instead of last 10 characters
            s += ' -J"%s' %(prefix_out[:10])
        ## job array
        if job_array_intervals:
            s += '[1-%i]' %(job_array_intervals,)
        ## job name, terminate
        s += '"'

        ## queue
        s += " -q %s" %(s_queue)
        ## memory
        s += " -M%i -R'select[mem>%i] rusage[mem=%i]'" %(
            memory_MB*1000, memory_MB, memory_MB,
            )
        ## project
        s += ' -P %s' %(self.s_project)
        ## stdout
##        s += '-o out/UnifiedGenotyper_$(date +%Y%m%d)_$(date +%H%M%S).out.%J.%I '
        s += ' -o stdout/%s' %(prefix_out,)
        if job_array_intervals:
##            s += '.%J.%I' ## job ID and job array index
            s += '.%I' ## job ID and job array index
        s += '.out'
        ## stderr
##        s += '-e err/UnifiedGenotyper_$(date +%Y%m%d)_$(date +%H%M%S).err.%J.%I '
        s += ' -e stderr/%s' %(prefix_out,)
        if job_array_intervals:
##            s += '.%J.%I' ## job ID and job array index
            s += '.%I' ## job ID and job array index
        s += '.err'

        ## sequential
        if bool_wait == True:
            s += ' -w "done($LSB_JOBID)"'
        
        ## command / shell script
        s += ' ./shell/%s.sh' %(prefix_out,)
        ## EOL
        s += '\n'

        if prefix_shell:
            ## initiate shell script
            lines = ['#!/bin/sh']
            ## append bsub line to shell script
            lines += [s]

            ## append new line breaks
            lines = [line+'\n' for line in lines]

            fp = 'shell/%s.sh' %(prefix_shell,)
            fd = open(fp,'w')
            fd.writelines(lines)
            fd.close()

            os.system('chmod +x %s' %(fp)) ## only works in unix environment...
##        ## Execute by owner
##        os.chmod(fn,stat.S_IEXEC)
##        ## Read, write, and execute by group
##        os.chmod(fn,stat.S_IRWXG)

        return s


    def bsub_all_chromosomes_all_steps(self,d_array,d_shell,):

        self.write_nonarray_scripts(d_shell,)

        self.write_array_scripts_chromosome(d_shell,)

        self.write_array_scripts_subchromosome(d_array,)

        if self.bool_sequential == True:
            self.write_initial_script()

        return


    def write_initial_script(self):

        ##
        ## write bsub wrapper; i.e. initiate submission of sequential steps
        ##
        bsub_main = '#!/bin/sh\n\n'

        ##
        ## 1) initiate with step 1
        ##
        prefix = self.append_prefix(self.l_steps[0])
        bsub_main += self.generate_bsub_line(
            prefix,
            bool_wait = False,
            )

        ##
        ## 2) append a check that all output has been generated before proceeding
        ##
        bsub_main += 'python ~/github/ms23/GATK_pipeline/check_output.py\n'

        ##
        ## 3) append step 2 following step 1
        ## e.g. CombineVariants following UnifiedGenotyper
        ##
        step_next = self.l_steps[0+1]
        s_queue = 'normal'
        memory_MB = 4000
        prefix_next = self.append_prefix(step_next)
        bsub_line_next = self.generate_bsub_line(
            prefix_next,
            s_queue = s_queue,
            memory_MB = memory_MB,
            bool_wait = False,
            )
        bsub_main += bsub_line_next

        prefix = 'bsub'
        prefix = self.append_prefix(prefix)
        fp = '%s.sh' %(prefix)
        fd = open(fp,'w')
        fd.write(bsub_main)
        fd.close()
        ## make bsub wrapper executable
        os.system('chmod +x %s' %(fp)) ## only works in unix environment...

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
                    prefix_next = self.append_prefix(step_next)
                    if (
                        (self.bool_join_chromosomes == False and step_next in self.l_steps_array_chromosome)
                        or
                        (self.bool_join_chromosomes == True and step_next == 'BEAGLE')
                        ):
                        if step_next in self.d_memory.keys():
                            s_queue = self.d_queues[step_next]
                            memory_MB = self.d_memory[step_next]
                        else:
                            s_queue = 'normal'
                            memory_MB = 4000
                        bsub_line_next = self.generate_bsub_line(
                            prefix_next,
                            job_array_intervals = 24,
                            job_array_title = '%s' %(step_next[:6],),
                            bool_wait = False,
                            s_queue = s_queue,
                            memory_MB = memory_MB,
                            )
                    else:
                        if (
                            self.bool_join_chromosomes == True
                            and
                            step_next == 'VariantRecalibrator'
                            ):
                            fn_out = self.d_file_outputs[step_next][0].replace('.$CHROMOSOME','')
                            dn_out = self.d_dir_outputs[step_next]
                            fp_out = os.path.join(dn_out,fn_out,)
                            s_queue = 'normal'
                            if os.path.isfile(fp_out):
                                memory_MB = 4000 ## temporary ugly solution...
                            else:
                                memory_MB = 12000
                        elif step_next == 'ApplyRecalibration' and self.bool_join_chromosomes == True:
##                            ## 32GB not enough
##                            s_queue = 'long'
##                            memory_MB = 32000
                            s_queue = 'hugemem'
                            memory_MB = 64000
                        else:
                            s_queue = 'normal'
                            memory_MB = 4000
                        bsub_line_next = self.generate_bsub_line(
                            prefix_next,
                            s_queue = s_queue,
                            memory_MB = memory_MB,
                            )
##                    d_bsub_lines[step_next] = bsub_line_next
            else:
                bsub_line_next = None

            ## skip array steps
            if step in self.l_steps_array_subchromosome+self.l_steps_array_chromosome:
                continue

            prefix = self.append_prefix(step)
            s = self.check_input_output_shell(step)
            s += d_shell[step]
            s += '\n\nfi\n\n'
            self.write_shell_script(
                prefix,
                s,
                line_additional = bsub_line_next,
                )
##            print d_shell[step]
##            stop

        return


    def write_array_scripts_chromosome(self,d_shell,):

        instance_GATK = GATK()

        if self.bool_join_chromosomes == True:
            chromosome_prefix = ''
        else:
            chromosome_prefix = '$CHROMOSOME'
        
        ##
        ## write shell script for UnifiedGenotyper (job array and first step)
        ##
        for step in self.l_steps_array_chromosome:
               
            if step == 'ApplyRecalibration' and self.bool_join_chromosomes == True:
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
            if step == 'ApplyRecalibration' and self.bool_join_chromosomes == True:
                s_condition = s_condition.replace('.$CHROMOSOME','')
            s += s_condition
            if step == 'ApplyRecalibration':
                s += instance_GATK.ApplyRecalibration(chromosome_prefix)
            elif step == 'BEAGLE':
##                s += self.BEAGLE(chromosome_prefix)
                s += self.BEAGLE('$CHROMOSOME')
            else:
                print step
                stop
            s += '\n\nfi\n'
            d_shell[step] = s

            ## append step following UnifiedGenotyper; e.g. CombineVariants
            if self.bool_sequential == True:
                step_next = self.l_steps[self.l_steps.index(step)+1]
                prefix_next = self.append_prefix(step_next)

                bsub_line_next = ''
                ## is LSB_JOBINDEX always equal to 1, if job is not submitted as an array?
                bool_wrap = False
                if not (
                    self.bool_join_chromosomes == True
                    and
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
                    prefix_next,
                    s_queue = 'normal',
                    memory_MB = 4000,
                    )
                if bool_wrap == True:
                    bsub_line_next += 'fi\n'

            prefix = self.append_prefix(step)
            self.write_shell_script(
                prefix,
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
            for chromosome in ['X','Y',]+range(1,22+1,):
##                if chromosome != 22: continue ## tmp
                chromosome = str(chromosome)
##                ## replace if statements with dictionary...
##                if step == 'UnifiedGenotyper':
##                    fp = 'out_GATK/UnifiedGenotyper.%s.vcf.1' %(chromosome)
##                elif step == 'IMPUTE2':
##                    fp = 'out_IMPUTE2/chromosome%s.impute2.1' %(chromosome)
##                else:
##                    print step
##                    print self.d_file_outputs[step]
##                    stop

##                dn_out = instance_main.d_dir_outputs[step]
##                fn_out = instance_main.d_file_outputs[step][0].replace('$CHROMOSOME',chromosome).replace('.$LSB_JOBINDEX','')
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
                s += '%s\n' %(d_array[step][chromosome])
####                if self.l_steps.index(step) != 0:
####                    s += 'fi\n'
                s += 'fi\n'

    ##        ## append step following UnifiedGenotyper; e.g. CombineVariants
    ##        if self.bool_sequential == True:
    ##            s_UG += d_bsub_lines[
    ##                self.l_steps[self.l_steps.index('UnifiedGenotyper')+1]
    ##                ]

            prefix = self.append_prefix(step)

            self.write_shell_script(prefix,s,)

        return


    def check_input_output_shell(self,step,bool_input=True,chromosome=None,):

        s = ''
        s += 'if [ '
        s_output = ' ! -s %s/%s ' %(
            self.d_dir_outputs[step],
            self.d_file_outputs[step][0].replace('$LSB_JOBINDEX','1'),
            )
        if (
            self.bool_join_chromosomes == True
            and
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
##            if self.bool_join_chromosomes == True and step in ['ProduceBeagleInput','vcf2gen',]:
            if (
                self.bool_join_chromosomes == True
                and
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
        s += 'then\n\n\n'

        if chromosome:
            s = s.replace('$CHROMOSOME',chromosome)

        return s


    def append_prefix(self,prefix,):

        if instance_main.bool_join_chromosomes == True:
            prefix += '_join'
        else:
            prefix += '_sep'

        return prefix


    def parse_chromosome_ranges(self):

        d_chromosome_lengths = {}

        ## 1000G chromosome ranges
        fd = open('%s.fai' %(self.fp_FASTA_reference_sequence))
        lines = fd.readlines()
        fd.close()
        for line in lines:
            l = line.strip().split()
            chromosome = l[0]
            chromosome_length = int(l[1])
##            if not chromosome in ['X','Y',]:
##                chromosome = int(chromosome)
            d_chromosome_lengths[chromosome] = chromosome_length
            ## break, when we reach the Y chromosome
            if chromosome == 'Y':
                break
        
        return d_chromosome_lengths


    def init_dirs(self):


        l_dn = ['stdout','stderr','out_GATK','shell','out_BEAGLE','out_Tommy','out_IMPUTE2',]

        cwd = os.getcwd()
        if cwd[cwd.rindex('/')+1:] in l_dn:
            bool_discontinue = True
            print 'You are in a subdir, so I will not run.'
        else:
            bool_discontinue = False
            ## create subdirs
            for dn in l_dn:
                if not os.path.isdir(dn):
                    os.mkdir(dn)

        return bool_discontinue


    def __init__(self,):

        ##
        ## GATK resource bundle
        ## http://www.broadinstitute.org/gsa/wiki/index.php/GATK_resource_bundle
        ##

        ## curl -u gsapubftp-anonymous ftp.broadinstitute.org/bundle/1.5/b37/human_g1k_v37.fasta.gz -o human_g1k_v37.fasta.gz; gunzip human_g1k_v37.fasta.gz
        self.fp_FASTA_reference_sequence = '/lustre/scratch111/resources/vrpipe/ref/Homo_sapiens/1000Genomes/human_g1k_v37.fasta'
        if not '--bamdir' in sys.argv:
            print 'specify --bamdir (e.g. \
\n/lustre/scratch111/projects/uganda/release/20120610/chromosome_bams \
\nor\n/lustre/scratch111/projects/uganda/release/downsample/20120921/chromosome_bams \
\n)'
            sys.exit(0)
        self.dn_BAM_input_file = sys.argv[sys.argv.index('--bamdir')+1]
##        self.dn_BAM_input_file = '/lustre/scratch111/projects/uganda/release/20120610/chromosome_bams'

        ## dbSNP file. rsIDs from this file are used to populate the ID column of the output.
##        self.dbsnp = '/lustre/scratch107/projects/uganda/users/tc9/in_GATK/dbsnp_135.b37.vcf'
        ## curl -u gsapubftp-anonymous ftp.broadinstitute.org/bundle/1.5/b37/dbsnp_135.b37.vcf.gz -o dbsnp_135.b37.vcf.gz; gunzip dbsnp_135.b37.vcf.gz
        self.dbsnp = '/lustre/scratch107/projects/uganda/users/tc9/in_GATK/ALL.wgs.dbsnp.build135.snps.sites.vcf'

        ## VariantRecalibrator
        self.vcf_hapmap = '/lustre/scratch107/projects/uganda/users/tc9/in_GATK/hapmap_3.3.b37.sites.vcf'
        self.vcf_omni25 = '/lustre/scratch107/projects/uganda/users/tc9/in_GATK/1000G_omni2.5.b37.sites.vcf'
        self.vcf_dbsnp = '/lustre/scratch107/projects/uganda/users/tc9/in_GATK/dbsnp_135.b37.vcf'

        ## IMPUTE2 input files
        self.hap = '/nfs/t149_1kg/ALL_1000G_phase1integrated_v3_impute/ALL_1000G_phase1integrated_v3_chr$CHROMOSOME_impute.hap.gz'
        self.legend = '/nfs/t149_1kg/ALL_1000G_phase1integrated_v3_impute/ALL_1000G_phase1integrated_v3_chr$CHROMOSOME_impute.legend.gz'

        ##
        ## software paths
        ##
        self.fp_GATK = '/software/varinf/releases/GATK/GenomeAnalysisTK-1.4-15-gcd43f01/GenomeAnalysisTK.jar'

        ## A node on the cluster does approximately 4Mbp in 8 hours;
        ## i.e. UnifiedGenotyper, ugandan dataset, chromosome1=249Mbp=3.1weeks=504hours,
        ## so do intervals of 4Mbp to avoid cluster time out after 12 hours
        ## this should be set as a flag with optparse as it differs from dataset to dataset...
        ## this one varies from chromosome to chromosome...
        ## todo: collect stats...
        ##
        ## the same variable is also used for IMPUTE2, so I decided to lower it from 10Mbp to 5Mbp
        self.bps_per_interval = 1.*10**6
        self.bps_per_interval = 1000.*10**6
        self.bps_per_interval = 100.*10**6
        self.bps_per_interval = 10.*10**6
        self.bps_per_interval = 10.*10**6
        self.bps_per_interval_IMPUTE2 = 5.*10**6

        self.bool_discontinue = self.init_dirs()

        ## sequential order in which to run GATK walkers and other steps
        self.l_steps = [
            'UnifiedGenotyper',
            'CombineVariants',
            'VariantRecalibrator',
            'determine_TS_level', ## not a GATK step
            'ApplyRecalibration',
##            'PrintReads', ## split VCFs by chromosome if bool_join_chromosomes == True
            ## Imputation, BEAGLE
            'ProduceBeagleInput',
            'BEAGLE', ## not a GATK step
##            'BeagleOutputToVCF',
            ## Imputation, IMPUTE2
##            'vcf2gen', ## if IMPUTE2 not preceded by BEAGLE
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
                ], ## not a GATK step
            'BeagleOutputToVCF':[
                'BeagleOutputToVCF.$CHROMOSOME.vcf',
                'BeagleOutputToVCF.$CHROMOSOME.vcf.idx',
                ],
            'IMPUTE2':[
                'chromosome$CHROMOSOME.impute2.gen.$LSB_JOBINDEX_summary',
                'chromosome$CHROMOSOME.impute2.gen.$LSB_JOBINDEX',
                ], ## not a GATK step
            }

        self.d_dir_outputs = {
            'UnifiedGenotyper':'out_GATK',
            'CombineVariants':'out_GATK',
            'VariantRecalibrator':'out_GATK',
            'determine_TS_level':'out_Tommy',
            'ApplyRecalibration':'out_GATK',
            'ProduceBeagleInput':'out_GATK',
            'BEAGLE':'out_BEAGLE',
            'BeagleOutputToVCF':'out_GATK',
            'IMPUTE2':'out_IMPUTE2',
            }

        self.d_queues = {
            'BEAGLE':'long'
            }

        self.d_memory = {
            'BEAGLE':16000,
            'determine_TS_level':100,
            'CombineVariants':4000,
            'UnifiedGenotyper':2000,
            'VariantRecalibrator':12000,
            }

        ## at which steps should a break be inserted
        ## i.e. when should a new shell script be submitted to the cluster
        ## in order to avoid hitting the wall time
        ## i.e. which of the steps are the time consuming ones
        ## I'm not using this at the moment, but I should...
        self.l_breaks = ['UnifiedGenotyper',]

        ##
        ##
        ##
        self.s_project = 'uganda'

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

        self.bool_join_chromosomes = True
##        self.bool_join_chromosomes = False ## tmp

        ## not currently used, but will be implemented...
        self.bool_send_mail_upon_job_completion = True

        for k in self.d_dir_outputs.keys():
            self.d_dir_outputs[k] = self.d_dir_outputs[k]

        return


class GATK():


    def BeagleOutputToVCF(self,chromosome,):

        '''
this script takes less than 10 minutes to run for chromosome 22,
so it should be combined with other chromosomes and/or steps
instead of being run in an independent bsub
'''

        analysis_type = 'BeagleOutputToVCF'

        lines = self.GATK_initiate(analysis_type,)

        ## GATKwalker, required, in
        lines += ['--beaglePhased:BEAGLE out_BEAGLE/BeagleOutput.%s.bgl.ProduceBeagleInput.%s.bgl.phased \\' %(chromosome,chromosome,)]
        lines += ['--beagleProbs:BEAGLE out_BEAGLE/BeagleOutput.%s.bgl.ProduceBeagleInput.%s.bgl.gprobs \\' %(chromosome,chromosome,)]
        lines += ['--beagleR2:BEAGLE out_BEAGLE/BeagleOutput.%s.bgl.ProduceBeagleInput.%s.bgl.r2 \\' %(chromosome,chromosome,)]
        ## GATKwalker, required, out
        fp_out = 'out_GATK/%s.%s.vcf' %(analysis_type,chromosome,)
        lines += ['--out %s \\' %(fp_out,)]
        ## GATKwalker, required, in
        lines += ['--variant out_GATK/ApplyRecalibration.recalibrated.filtered.%s.vcf \\' %(chromosome)]

        s = self.GATK_terminate(lines,fp_out,)

        return s


    def ProduceBeagleInput(self,chromosome,):

        '''this walker takes approximately 5-10 minutes per chromosome to run'''

        ## http://www.broadinstitute.org/gsa/wiki/index.php/Interface_with_BEAGLE_imputation_software
        ## "After variants are called and possibly filtered,
        ## the GATK walker ProduceBeagleInputWalker will take the resulting VCF as input,
        ## and will produce a likelihood file in BEAGLE format.

        if chromosome != '':
            chromosome_suffix = '.%s' %(chromosome)
        else:
            chromosome_suffix = ''

        analysis_type = 'ProduceBeagleInput'

        lines = self.GATK_initiate(analysis_type,)

        ## GATKwalker, required, out
        fp_out = 'out_GATK/%s%s.bgl' %(analysis_type,chromosome_suffix,)
        lines += ['--out %s \\' %(fp_out,)]
        ## GATKwalker, required, in
        lines += ['--variant out_GATK/ApplyRecalibration.recalibrated.filtered%s.vcf \\' %(chromosome_suffix)]
##        if instance_main.bool_join_chromosomes == True:
##            ## GATKwalker, optional, out
##            lines += ['--intervals %s \\' %(chromosome)]

        s = self.GATK_terminate(lines,fp_out,)

        return s


    def GATK_terminate(self,lines,fp_out,):

        ## make sure there is a white space between parameters
        lines = [' '+line for line in lines]
        ## write shell script
        s = '\n'.join(lines)+'\n\n'

        ## replace this with a file check in the shell script...
        if instance_main.bool_skip_if_output_exists == True:
            fp_out = fp_out.replace('$LSB_JOBINDEX','1') ## just do first part of array...
            if (
                instance_main.bool_join_chromosomes == False
                and
                instance_main.bool_skip_if_output_exists == True
                and
                os.path.isfile(fp_out)
                and
                os.path.getsize(fp_out) > 0
                ):
                s = ''

        return s


    def ApplyRecalibration(self,chromosome,):

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

        analysis_type = 'ApplyRecalibration'

        if chromosome != '':
            chromosome_suffix = '.%s' %(chromosome)
        else:
            chromosome_suffix = ''

        fp_in = 'out_Tommy/determine_TS_level%s.txt' %(chromosome_suffix)

##        ## something went wrong in the previous step...
##        ## avoid raising an error, if something went wrong in the previous step...
##        ## NO! Need to always run this step!!!
##        if not os.path.isfile(fp_in):
##            return ''
##        if os.path.getsize(fp_in) == 0:
##            return ''

        lines = ['ts_filter_level=($(cat %s))\n' %(fp_in)]

        lines += self.GATK_initiate(analysis_type,)

##        fp_out = 'out_GATK/ApplyRecalibration.recalibrated.filtered.%s.vcf' %(chromosome)
        dn_out = instance_main.d_dir_outputs['ApplyRecalibration']
        fn_out = instance_main.d_file_outputs['ApplyRecalibration'][0].replace('.$CHROMOSOME',chromosome_suffix)
        fp_out = os.path.join(dn_out,fn_out,)

        ## GATKwalker, required, in
        lines += ['--input out_GATK/CombineVariants%s.vcf \\' %(chromosome_suffix,)]
        ## GATKwalker, required, out
        lines += ['--out %s \\' %(fp_out,)]
        ## GATKwalker, required, in
        lines += ['--recal_file out_GATK/VariantRecalibrator%s.recal \\' %(chromosome_suffix)]
        lines += ['--tranches_file out_GATK/VariantRecalibrator%s.tranches \\' %(chromosome_suffix)]
        ## GATKwalker, optional, in
        ## ts_filter_level should correspond to targetTruthSensitivity prior to novelTiTv dropping (see tranches plot)
        ## default is 99.00
##        lines += ['--ts_filter_level 99.0 \\']
        lines += ['--ts_filter_level $ts_filter_level \\']

        s = self.GATK_terminate(lines,fp_out,)

        return s


    def VariantRecalibrator(self,chromosome,):

        ## http://www.broadinstitute.org/gsa/gatkdocs/release/org_broadinstitute_sting_gatk_walkers_variantrecalibration_VariantRecalibrator.html

        if chromosome != '':
            chromosome_suffix = '.%s' %(chromosome)
        else:
            chromosome_suffix = ''

        ## http://www.broadinstitute.org/gsa/wiki/index.php/Best_Practice_Variant_Detection_with_the_GATK_v3#Whole_Genome_Shotgun_experiments

        ## http://www.broadinstitute.org/gsa/wiki/index.php/Variant_quality_score_recalibration
        ## Can I use the variant quality score recalibrator with my small sequencing experiment?
        ## This tool is expecting thousands of variant sites in order to achieve decent modeling with the Gaussian mixture model. Whole exome call sets work well, but anything smaller than that scale might run into difficulties. One piece of advice is to turn down the number of Gaussians used during training and to turn up the number of variants that are used to train the negative model. This can be accomplished by adding --maxGaussians 4 --percentBad 0.05 to your command line.
        ## http://www.broadinstitute.org/gsa/gatkdocs/release/org_broadinstitute_sting_gatk_walkers_variantrecalibration_VariantRecalibrator.html#--maxGaussians
        ## http://www.broadinstitute.org/gsa/gatkdocs/release/org_broadinstitute_sting_gatk_walkers_variantrecalibration_VariantRecalibrator.html#--percentBadVariants

        analysis_type = 'VariantRecalibrator'

        lines = self.GATK_initiate(analysis_type)

        ## GATKwalker, required, in
        lines += ['--input out_GATK/CombineVariants%s.vcf \\' %(chromosome_suffix,)]
        ## GATKwalker, required, out
        fp_out = 'out_GATK/%s%s.recal' %(analysis_type,chromosome_suffix,)
        lines += ['--recal_file out_GATK/VariantRecalibrator%s.recal \\' %(chromosome_suffix,)]
        lines += ['--tranches_file out_GATK/VariantRecalibrator%s.tranches \\' %(chromosome_suffix,)]
        ## GATKwalker, required, in (ask Deepti about this...)
        lines += ['--use_annotation QD --use_annotation HaplotypeScore\
        --use_annotation MQRankSum --use_annotation ReadPosRankSum --use_annotation MQ\
        --use_annotation FS --use_annotation DP \\']
        ## GATKwalker, optional, in
        ## Ugandan QCed to be added..!
        lines += ['-resource:hapmap,known=false,training=true,truth=true,prior=15.0 \
%s \\' %(instance_main.vcf_hapmap)] ## ask Deepti where the master copy is located (GATK resource bundle)
        lines += ['-resource:omni,known=false,training=true,truth=false,prior=12.0 \
%s \\' %(instance_main.vcf_omni25)] ## ask Deepti which vcf file to use as input
        lines += ['-resource:dbsnp,known=true,training=false,truth=false,prior=8.0 \
%s \\' %(instance_main.vcf_dbsnp)] ## ask Deepti which vcf file to use as input

        s_TStranches = ''
        l_TStranches = [100,]
        l_TStranches += [99+i/10. for i in range(9,0,-1,)]
        l_TStranches += [90+i/2. for i in range(18,-1,-1,)]
##        l_TStranches += [80+i/2. for i in range(9,-1,-1,)]
##        print l_TStranches
##        stop
        for TStranche in l_TStranches:
            s_TStranches += '--TStranche %.1f ' %(TStranche)
        lines += ['%s \\' %(s_TStranches)]
##        lines += ['--TStranche 100.0 --TStranche 99.9 --TStranche 99.8 --TStranche 99.5 --TStranche 99.0 --TStranche 95.0 --TStranche 90.0 \\']

        ## GATKwalker, optional, out
##        ## generating the VariantRecalibrator.plots.R.pdf file is time consuming, when you are eager to see your results...
##        lines += ['--rscript_file out_GATK/VariantRecalibrator.%s.plots.R \\' %(chromosome,)]

        s = self.GATK_terminate(lines,fp_out,)

        if chromosome == 'Y': s = '' ## tmp (Y chromosome fails because VariantRecalibrator is meant to be run on the whole genome and not separate chromosomes)

        return s


    def GATK_initiate(self,analysis_type,):

        ## run GATK
        s = 'java '
        ## set maximum heap size
        if analysis_type == 'VariantRecalibrator' and instance_main.bool_join_chromosomes == True:
            s += '-Xmx12g ' ## Max Memory :     11856 MB
        elif analysis_type == 'ApplyRecalibration' and instance_main.bool_join_chromosomes == True:
            s += '-Xmx62g ' ## 28GB not enough
        else:
            s += '-Xmx4g '
##        s += '-Xmx4g '
        s += '-jar %s \\' %(self.fp_GATK)
        lines = [s]

        ## CommandLineGATK, required, in
        lines += ['--analysis_type %s \\' %(analysis_type)]
        ## CommandLineGATK, optional, in
        lines += ['--reference_sequence %s \\' %(instance_main.fp_FASTA_reference_sequence)]

        return lines


    def CombineVariants(self,chromosome,d_chromosome_lengths,):

        '''write CombineVariants shell script.
It takes a while to run, so it really shouldn't be run sequentially for each chromosome...
But there is a lot of file I/O, so it's probably safer to run it one chromosome at a time...
It takes 5.5 hours to run for the entire Uganda exome (100 samples)
'''

        analysis_type = 'CombineVariants'

        ## http://www.broadinstitute.org/gsa/wiki/index.php/CombineVariants

        ## how many intervals to cover the entire range of the chromosome?
        intervals = int(math.ceil(d_chromosome_lengths[chromosome]/instance_main.bps_per_interval))

        if instance_main.bool_join_chromosomes == True:
            lines = []
        else:
            lines = self.GATK_initiate(analysis_type,)
        
        ##
        ## required
        ##

        ## File to which variants should be written
        if instance_main.bool_join_chromosomes == False:
            fp_out = 'out_GATK/%s.%s.vcf' %(analysis_type,chromosome,)
            lines += [' --out %s \\' %(fp_out)]

        ## Input VCF file
        for interval in range(1,intervals+1,):
            lines += [
                ' --variant out_GATK/UnifiedGenotyper.%s.vcf.%i \\' %(
                    chromosome, interval,
                    )
                ]

        ## File to which variants should be written
        if instance_main.bool_join_chromosomes == True:
            fp_out = 'out_GATK/%s.vcf' %(analysis_type,)

        s = self.GATK_terminate(lines,fp_out,)

        return s


    def UnifiedGenotyper(self,chromosome,d_chromosome_lengths,):

        ## http://www.broadinstitute.org/gsa/wiki/index.php/Unified_genotyper
        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_genotyper_UnifiedGenotyper.html

        analysis_type = 'UnifiedGenotyper'

        ##
        ## initiate shell script
        ##
        lines = []
        lines += ['#!/bin/sh']

        ##
        ## do not allow interval to exceed the length of the chromosome
        ## otherwise it will raise an I/O error in GATK 1.4-15
        ##
        lines += ['max=$((($LSB_JOBINDEX+0)*%i))' %(int(instance_main.bps_per_interval))]
        lines += ['if test $max -ge %i' %(d_chromosome_lengths[chromosome])]
        lines += ['then max=%i' %(d_chromosome_lengths[chromosome])]
        lines += ['fi']

        lines += self.GATK_initiate(analysis_type,)

        ##
        ## required
        ##
        ## File to which variants should be written
##        lines += [' --out UnifiedGenotyper_$(date +%Y%m%d)_$(date +%H%M%S).vcf.$LSB_JOBINDEX \\']
        fp_out = 'out_GATK/UnifiedGenotyper.%s.vcf.$LSB_JOBINDEX' %(chromosome)
        lines += [' --out %s \\' %(fp_out)]

        ##
        ## CommandLineGATK, optional
        ##
        lines += [' --input_file %s/chrom%s.bam \\' %(instance_main.dn_BAM_input_file,chromosome,)]

##        lines += [' --intervals chr%i:%i-%i \\']
##        lines += [
##            ' --intervals %s:$((($LSB_JOBINDEX-1)*%i+1))-$((($LSB_JOBINDEX+0)*%i)) \\' %(
##                chromosome,int(self.bps_per_interval),int(self.bps_per_interval),
##                )
##            ]
##        lines += [
##            ' --intervals chromosome%s.interval_list \\' %(
##                chromosome,
##                )
##            ]
        lines += [
            ' --intervals %s:$((($LSB_JOBINDEX-1)*%i+1))-$max \\' %(
                chromosome,int(instance_main.bps_per_interval),
                )
            ]

        ##
        ## UnifiedGenotyper, optional
        ##
        ## dbSNP file. rsIDs from this file are used to populate the ID column of the output.
        lines += [' --dbsnp %s \\' %(instance_main.dbsnp)]
        ## http://www.broadinstitute.org/gsa/wiki/index.php/Best_Practice_Variant_Detection_with_the_GATK_v3#Selecting_an_appropriate_quality_score_threshold
        ## http://www.broadinstitute.org/gatk/guide/topic?name=best-practices
        ## Selecting an appropriate quality score threshold
        ## A common question is the confidence score threshold to use for variant detection. We recommend:
        ## Deep (> 10x coverage per sample) data:
        ## we recommend a minimum confidence score threshold of Q30.
        ## Shallow (< 10x coverage per sample) data:
        ## because variants have by necessity lower quality with shallower coverage we recommend
        ## a minimum confidence score of Q4 in projects with 100 samples or fewer
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

        ##
        ## write shell script
        ##
        s = self.GATK_terminate(lines,fp_out,)

        ## unlike the other GATK walkers,UnifiedGenotyper.%
        ## UG requires individual shell scripts,
        ## because job arrays are used for each chromosome
        fp = 'shell/UnifiedGenotyper.%s.sh' %(chromosome)
        fd = open(fp,'w')
        fd.write(s)
        fd.close()
        os.system('chmod +x %s' %(fp)) ## only works in unix environment...

        return


    def __init__(self,):

        self.fp_GATK = '/software/varinf/releases/GATK/GenomeAnalysisTK-1.4-15-gcd43f01/GenomeAnalysisTK.jar'

        return


if __name__ == '__main__':
    instance_main = main()
    instance_main.main()
