#!/software/bin/python

## T. Carstensen (tc9), M.S. Sandhu (ms23), D. Gurdasani (dg11)
## Wellcome Trust Sanger Institute, 2012

import os, sys, time
sys.path.append(os.path.dirname(sys.argv[0]))
import QC
sys.path.append('/nfs/users/nfs_t/tc9/github/sandbox')
import gnuplot


hours = 8


def main():

    l_populations = [
        'Ethiopia', ## 'Amhara', 'Somali', 'Oromo',
        'Fula', 'Ga-Adangbe', 'Ibo', 'Jola', 'Kalenjin', 'Kikuyu', 'Mandinka', 'Muganda_octo', 'Muganda_quad', 'Munyarwanda_octo', 'Munyarwanda_quad', 'Murundi', 'Sotho', 'Wolloff', 'Zulu',
        'YRI','TSI','PUR','PEL','MXL','MKK','LWK','KHV','JPT','IBS','IBO','GIH','GBR','FIN','CLM','CHS','CHD','CHB','CEU','CDX','ASW','ACB',
        ]
    l_populations.sort()


    allelic_concordance_Baganda(
##        'pops/Baganda_quad/Baganda_quad', ## tmp
##        'pops/Baganda_octo/Baganda_octo', ## tmp
##        'pops/Baganda_quad/Baganda_quad.sampleQC', ## tmp
##        'pops/Baganda_octo/Baganda_octo.sampleQC', ## tmp
        'pops/Baganda_quad/Baganda_quad.SNPQC',
        'pops/Baganda_octo/Baganda_octo.SNPQC',
        )

    frq_discordant()    

##    MDS_Baganda()
##    sys.exit(0)


##    tables(l_populations)
##    plots(l_populations)
##    sys.exit(0)


    return


def frq_discordant():

    for chip in ['quad','octo',]:

        cmd = 'cat pops/Baganda_%s/Baganda_%s.SNPQC.fam | wc -l' %(chip,chip,)
        n_samples = int(os.popen(cmd).read())

        for flag,prefix in [
            ['','discordant',],
            ['-v','concordant',],
            ]:

            for suffix,column in [
                ['frq','$5',]
                ['lmiss','$5',],
                ]:

                cmd = 'fgrep %s -w -f disconcordant.SNPs pops/Baganda_%s/Baganda_%s.SNPQC.%s' %(
                    flag,chip,chip,suffix,)
                cmd += " | awk '{print %s}' > %s_%s.%s" %(
                    column,prefix,chip,suffix,)
                execmd(cmd)

                cmd = 'cat %s_%s.%s | wc -l' %(prefix,chip,suffix,)
                n_SNPs = int(os.popen(cmd).read())

                gnuplot.histogram2(
                    '%s_%s.%s' %(prefix,chip,suffix,),
                    x_step = 0.01,
                    xlabel='MAF after sample QC',
                    title='Baganda %s (n_{samples}=%i, n_{SNPs}=%i)' %(
                        chip,n_samples,n_SNPs,
                        ),
                    )

    return


def execmd(cmd):

    if cmd.split()[0] == 'cat':
        if not os.path.isfile(cmd.split()[1]):
            print cmd
            print 'does not exist:', cmd.split()[1]
            sys.exit(0)
    if cmd.split()[0] in ['grep','fgrep',] and cmd.split()[1] == '-f':
        if not os.path.isfile(cmd.split()[2]):
            print cmd
            print 'does not exist:', cmd.split()[2]
            sys.exit(0)
    print cmd
    os.system(cmd)

    return


def MDS_Baganda():

    import gnuplot

    prefix = 'Baganda29_quad29_octo200'
    prefix = 'Baganda29_quad29_octo200_excldisconcordant'
    prefix = 'Baganda29_quad29_octo200_SNPQCtogether'

    fn_ld_regions = 'pops/Baganda_octo/Baganda_octo.ldregions.SNPs'

    ##
    ## find common SNPs post QC
    ##
    if not os.path.isfile('%s.extract' %(prefix)):
        if prefix == 'Baganda29_quad29_octo200_SNPQCtogether':
            cmd = 'cat pops/Baganda_quad29octo200/Baganda_quad29octo200.SNPQC.bim > %s.extract' %(prefix)
            execmd(cmd)
        else:
            cmd = "cat pops/Baganda_quad/Baganda_quad.SNPQC.bim | awk '{print $2}' | sort > Baganda_quad.SNPs"
            execmd(cmd)
            cmd = "cat pops/Baganda_octo/Baganda_octo.SNPQC.bim | awk '{print $2}' | sort > Baganda_octo.SNPs"
            execmd(cmd)
            cmd = 'comm -12 Baganda_quad.SNPs Baganda_octo.SNPs > %s.extract' %(prefix)
            execmd(cmd)
            os.remove('Baganda_quad.SNPs')
            os.remove('Baganda_octo.SNPs')

    ##
    ## find common samples post QC
    ##
    cmd = 'cat pops/Baganda_quad/Baganda_quad.SNPQC.fam | sort > Baganda_quad.fam'
    execmd(cmd)
    cmd = 'cat pops/Baganda_octo/Baganda_octo.SNPQC.fam | sort > Baganda_octo.fam'
    execmd(cmd)
    
    cmd = "cat Baganda_quad.fam | awk '{print substr($1,12,10)}' | sort > Baganda_quad.samples"
    execmd(cmd)
    cmd = "cat Baganda_octo.fam | awk '{print substr($1,12,10)}' | sort > Baganda_octo.samples"
    execmd(cmd)
    cmd = 'comm -12 Baganda_quad.samples Baganda_octo.samples > Baganda29.samples'
    execmd(cmd)
    os.remove('Baganda_quad.samples')
    os.remove('Baganda_octo.samples')
    cmd = 'fgrep -f Baganda29.samples Baganda_quad.fam | sort > Baganda29_quad29_octo0.fam'
    execmd(cmd)
    cmd = 'fgrep -f Baganda29.samples Baganda_octo.fam | sort > Baganda29_quad0_octo29.fam'
    execmd(cmd)
    os.remove('Baganda29.samples')

    cmd = 'comm -23 Baganda_quad.fam Baganda29_quad29_octo0.fam > Baganda29_quad71_octo0.fam'
    execmd(cmd)
    cmd = 'comm -23 Baganda_octo.fam Baganda29_quad0_octo29.fam > Baganda29_quad0_octo200.fam'
    execmd(cmd)
    cmd = 'cat Baganda29_quad0_octo200.fam Baganda29_quad29_octo0.fam > Baganda29_quad29_octo200.fam'
    execmd(cmd)
    cmd = 'cat Baganda29_quad0_octo200.fam Baganda29_quad29_octo0.fam > Baganda29_quad29_octo200_excldisconcordant.fam'
    execmd(cmd)
    if prefix == 'Baganda29_quad29_octo200_SNPQCtogether':
        cmd = 'cat pops/Baganda_quad29octo200/Baganda_quad29octo200.fam > Baganda29_quad29_octo200_SNPQCtogether.fam'
        execmd(cmd)

    ##
    ## --bmerge
    ##
    if not os.path.isfile('%s.bed' %(prefix)):
        cmd = 'plink \\\n'
        if prefix == 'Baganda29_quad29_octo200_SNPQCtogether':
            cmd += '--bfile pops/Baganda_quad29octo200/Baganda_quad29octo200.SNPQC \\\n'
        else:
            cmd += '--bfile pops/Baganda_quad/Baganda_quad.SNPQC \\\n'
            cmd += '--bmerge \\\n'
            cmd += 'pops/Baganda_octo/Baganda_octo.SNPQC.bed \\\n'
            cmd += 'pops/Baganda_octo/Baganda_octo.SNPQC.bim \\\n'
            cmd += 'pops/Baganda_octo/Baganda_octo.SNPQC.fam \\\n'
        cmd += '--keep %s.fam \\\n' %(prefix)
        cmd += '--extract %s.extract \\\n' %(prefix)
##        cmd += '--exclude 26diff_and_monomorphic.SNPs \\\n' ## tmp
        if 'excldisconcordant' in prefix:
            cmd += '--exclude disconcordant.SNPs \\\n' ## tmp
        cmd += '--make-bed --out %s \\\n' %(prefix)
        execmd(cmd)

    ##
    ## --indep-pairwise
    ##
    if not os.path.isfile('%s.prune.in' %(prefix)):
        cmd = 'plink \\\n'
        cmd += '--bfile %s \\\n' %(prefix)
        cmd += '--out %s \\\n' %(prefix)
        ## settings
        cmd += '--indep-pairwise 50 5 0.2 \\\n'
        cmd += '--maf 0.05 \\\n'
        ## SNP exclusion
        cmd += '--exclude %s \\\n' %(fn_ld_regions)
        execmd(cmd)

    ##
    ## --genome
    ##
    if not os.path.isfile('%s.genome' %(prefix)):    
        cmd = 'plink \\\n'
        cmd += '--bfile %s \\\n' %(prefix)
        cmd += '--out %s \\\n' %(prefix)
        cmd += '--genome \\\n'
        ## SNP exclusion
        cmd += '--extract %s.prune.in \\\n' %(prefix)
        cmd += '--exclude %s \\\n' %(fn_ld_regions)
        execmd(cmd)

    ##
    ## --cluster
    ##
    if not os.path.isfile('%s.mds' %(prefix)):
        cmd = 'plink \\\n'
        cmd += '--bfile %s \\\n' %(prefix)
        cmd += '--out %s \\\n' %(prefix)
        cmd += '--cluster \\\n'
        cmd += '--mds-plot 4 \\\n'
        cmd += '--read-genome %s.genome \\\n' %(prefix)
        ## SNP exclusion
        cmd += '--extract %s.prune.in \\\n' %(prefix)
        cmd += '--exclude %s \\\n' %(fn_ld_regions)
        execmd(cmd)

##    ##
##    ## EIGENSOFT
##    ##
##    eigensoft(prefix,fn_ld_regions,)

    ##
    ## plot
    ##
    if not os.path.isfile('%s.mds' %(prefix)):
        sys.exit(0)

    cmd = "cat Baganda29_quad29_octo0.fam | awk '{print $1}' > Baganda29_quad29.samples"
    execmd(cmd)
    cmd = "cat Baganda29_quad0_octo29.fam | awk '{print $1}' > Baganda29_octo29.samples"
    execmd(cmd)
    for suffix in ['quad','octo',]:
        cmd = 'fgrep -f Baganda29_%s29.samples %s.mds' %(suffix,prefix,)
        cmd += " | awk '{print substr($1,12,10),$4,$5}'"
        cmd += ' | sort -k1,1'
        cmd += ' > %s_%s29.mds' %(prefix,suffix)
        execmd(cmd)
        if suffix == 'quad':
            continue
        cmd = "cat pops/Baganda_%s/Baganda_%s.fam | awk '{print $1}' > Baganda29_%s.samples" %(
            suffix,suffix,suffix,)
        execmd(cmd)
        cmd = 'fgrep -f Baganda29_%s.samples %s.mds' %(suffix,prefix,)
        cmd += " | awk '{print substr($1,12,10),$4,$5}'"
        cmd += ' | sort -k1,1'
        cmd += ' > %s_%s.mds' %(prefix,suffix)
        execmd(cmd)

##    lines_extra = ['set key out\n']
##    fd = open('%s_quad29.mds' %(prefix))
##    lines4 = fd.readlines()
##    fd.close()
##    fd = open('%s_octo29.mds' %(prefix))
##    lines8 = fd.readlines()
##    fd.close()
##    for i in xrange(len(lines4)):
##        l4 = lines4[i].split()
##        l8 = lines8[i].split()
##        x4 = float(l4[1])
##        y4 = float(l4[2])
##        x8 = float(l8[1])
##        y8 = float(l8[2])
##        lines_extra += ['set arrow from %f,%f to %f,%f\n' %(x4,y4,x8,y8,)]

    n_samples = int(os.popen('cat %s.mds | wc -l' %(prefix)).read())-1
    ## without pruning
    n_SNPs = int(os.popen('cat %s.bim | wc -l' %(prefix)).read())
    ## with pruning
    n_SNPs = int(os.popen('cat %s.prune.in | wc -l' %(prefix)).read())

    line_plot = 'plot '
    line_plot += '"%s_quad29.mds" u 2:3 ps 2 pt 7 lc 1 t "quad",' %(prefix)
##    line_plot += '"%s_octo29.mds" u 2:3 ps 3 pt 7 lc 3 t "octo",' %(prefix)
    line_plot += '"%s_octo.mds" u 2:3 ps 2 pt 7 lc rgb "#0000FF" t "octo",' %(prefix)
    line_plot = line_plot[:-1]

    gnuplot.scatter_plot_2d(
        '%s.mds' %(prefix),
        line_plot = line_plot,
##        column1 = 2, column2 = 3,
        xlabel = 'C1',
        ylabel = 'C2',
        title='%s (n_{samples}=%i, n_{SNPs}=%i)' %(
            'Baganda29',n_samples,n_SNPs,
            ),
        prefix_out='%s.mds' %(prefix),
##        lines_extra=lines_extra,
        )

    return


def eigensoft(prefix,fn_ld_regions,):

    bfile = prefix
    out_prefix = prefix

    cmd = 'cat %s.fam' %(out_prefix)
    cmd += ' | '
    cmd += "awk '{print $1,$2,substr($1,12,10),substr($2,12,10)}'"
    cmd += ' > '
    cmd += '%s.recoded.txt' %(out_prefix)
    execmd(cmd)

    cmd = 'plink \\\n'
    cmd += '--bfile %s \\\n' %(out_prefix,)
    cmd += '--exclude %s \\\n' %(fn_ld_regions,)
    cmd += '--extract %s.extract \\\n' %(out_prefix)
##    cmd += '--extract %s.prune.in \\\n' %(prefix)
    cmd += '--update-ids %s.recoded.txt \\\n' %(out_prefix,)
    cmd += '--make-bed \\\n'
    cmd += '--out %s.longLDexcl \\\n' %(out_prefix)
    execmd(cmd)

    path_EIG = '/nfs/team149/Software/EIG4.2/bin'
##    path_EIG = '/software/varinf/bin/eigensoft/bin'
##    path_EIG = '/nfs/team149/Software/usr/share/eigensoft/EIG4.2/bin'
    cmd = '%s/smartpca.perl \\\n' %(path_EIG)
    ## genotype file in any format (see ../CONVERTF/README)
    cmd += '-i %s.longLDexcl.bed \\\n' %(os.path.join(os.getcwd(),out_prefix))
    ## snp file in any format (see ../CONVERTF/README)
    cmd += '-a %s.longLDexcl.bim \\\n' %(os.path.join(os.getcwd(),out_prefix))
    ## indiv file in any format (see ../CONVERTF/README)
    cmd += '-b %s.longLDexcl.fam \\\n' %(os.path.join(os.getcwd(),out_prefix))
    ## (Default is 10) number of principal components to output
    cmd += '-k 10 \\\n'
    ## output file of principal components.  Individuals removed
    ## as outliers will have all values set to 0.0 in this file.
    cmd += '-o %s.pca \\\n' %(out_prefix.replace('.','_'))
    ## prefix of output plot files of top 2 principal components.
    ## (labeling individuals according to labels in indiv file)
    cmd += '-p %s.plot \\\n' %(out_prefix.replace('.','_'))
    ## output file of all eigenvalues
    cmd += '-e %s.eval \\\n' %(out_prefix.replace('.','_'))
    ## output logfile        
    cmd += '-l %s.eigensoft.log \\\n' %(out_prefix.replace('.','_'))
    ## (Default is 5) maximum number of outlier removal iterations.
    ## To turn off outlier removal, set -m 0.
    cmd += '-m 5 \\\n'
    ## (Default is 10) number of principal components along which
    ## to remove outliers during each outlier removal iteration.
    cmd += '-t 10 \\\n'
    ## (Default is 6.0) number of standard deviations which an
    ## individual must exceed, along one of topk top principal
    ## components, in order to be removed as an outlier.
    cmd += '-s 6 \\\n'
##    ## verbose
##    cmd += '-V \\\n'
    cmd += '-snpweightoutname %s.snpweight \\\n' %(out_prefix)
    cmd += '-badsnpname %s.badsnpname \\\n' %(out_prefix)
##    cmd += ' -p Muganda29_octo200_quad29.pca.par'
    print cmd
    os.system(cmd)

    return

def tmp_check_allele_concordance_bim_strand(
    l,strand_chr,strand_pos,l_set_alleles,fdstrand8,):

    while not (l[0] == strand_chr and l[3] == strand_pos):
        line_strand = fdstrand8.readline()
        l_strand = line_strand.split()
        strand_chr = l_strand[1]
        strand_pos = l_strand[2]
##        if rsID_strand == 'kgp7547300':
##            print rsID
##            stop
    allele_strand = list(l_strand[5])
    if l_strand[4] == '-':
        d = {'A':'T','G':'C','T':'A','C':'G'} ## move outside of loop if this becomes a permanent part of the code...
        for i in xrange(2):
            allele_strand[i] = d[allele_strand[i]]
    if len(set(l_set_alleles)-set(allele_strand)) > 0:
        print strand_chr, strand_pos
        print l_strand, allele_strand, l_set_alleles
        print l
        stop

    return strand_chr, strand_pos


def allelic_concordance_Baganda(bfile1,bfile2,):

    for bfile in [bfile1,bfile2,]:
        if not os.path.isfile('%s.bed' %(bfile)):
            print 'bed not found:', bfile
            sys.exit(0)

    import math

    basename1 = os.path.basename(bfile1)
    basename2 = os.path.basename(bfile2)

    ##
    ## find common samples
    ##
    fn_common_samples = '%s.%s.comm.samples' %(basename1,basename2,)
    cmd = "cat %s.fam | awk '{print substr($1,12,10)}'" %(bfile1)
    l_quad = os.popen(cmd).read().strip().split('\n')
    cmd = "cat %s.fam | awk '{print substr($1,12,10)}'" %(bfile2)
    l_octo = os.popen(cmd).read().strip().split('\n')
    l = list(set(l_quad)&set(l_octo))
    n_samples = len(l)
    s = '\n'.join(l)
    fd = open('%s' %(fn_common_samples),'w')
    fd.write(s)
    fd.close()

    ##
    ## find common SNPs
    ##
    fn_common_SNPs = '%s.%s.comm.SNPs' %(basename1,basename2,)
    if not os.path.isfile(fn_common_SNPs):
        for bfile, basename in [[bfile1,basename1,],[bfile2,basename2,],]:
            cmd = "cat %s.bim | awk '{print $2}' | sort > %s.SNPs" %(bfile,basename,)
            execmd(cmd)
        cmd = 'comm -12 %s.SNPs %s.SNPs > %s' %(
            basename1,basename2,fn_common_SNPs,)
        execmd(cmd)
        pass

    ##
    ## find common samples
    ##
    for bfile in [bfile1,bfile2,]:
        basename = os.path.basename(bfile)
        cmd = 'fgrep -f %s %s.fam' %(
            fn_common_samples,bfile,)
        cmd += " | awk '{print $1,$2}'"
        cmd += ' > %s.comm.samples' %(basename)
        execmd(cmd)
        pass

    ##
    ## recode from bed to ped and extract common samples and SNPs
    ##
    for bfile in [bfile1,bfile2,]:
        basename = os.path.basename(bfile)
        if os.path.isfile('%s.tped' %(basename)):
            deltat = time.time()-os.path.getmtime('%s.tped' %(basename))
            if deltat < hours*60*60:
                print int(deltat/float(60*60)), 'hrs'
                continue
        cmd = 'plink --bfile %s ' %(bfile)
        cmd += '--keep %s.comm.samples ' %(basename,)
        cmd += ' --recode --transpose --out %s' %(basename,)
        cmd += ' --extract %s' %(fn_common_SNPs)
        execmd(cmd)

        cmd = 'sort -k1,1 -k4,4 %s.tped -o %s.tped' %(basename,basename,)
        execmd(cmd)

    d_fam = {}
    for bfile in [bfile1,bfile2,]:
        basename = os.path.basename(bfile)
        cmd = "cat %s.tfam | awk '{print substr($1,12,10)}'" %(basename)
        l = os.popen(cmd).read().strip().split('\n')
        d_fam[basename] = l
    l_fam4 = d_fam[basename1]
    l_fam8 = d_fam[basename2]

    d_stats = {'sumx':{},'sumy':{},'sumxx':{},'sumyy':{},'sumxy':{},}
    for sample in l_fam4:
        d_stats['sumx'][sample] = 0
        d_stats['sumy'][sample] = 0
        d_stats['sumxx'][sample] = 0
        d_stats['sumyy'][sample] = 0
        d_stats['sumxy'][sample] = 0

    count_quad0octo1 = 0
    count_quad1octo0 = 0
    count_quad0octo0 = 0
    count_quad1octo1 = 0
    d_diff = {'monomorphic':{},'all':{},'nonmonomorphic':{},}
    for k in  d_diff.keys():
        for i in xrange(n_samples+1):
            d_diff[k][i] = 0
    i_line = 0
    sum_xy = 0
    sum_xx = 0
    sum_yy = 0
    sum_x = 0
    sum_y = 0
    n = 0

    l_monomorphic = []
    l_26diff_and_monomorphic = []
    l_disconcordant = []

    fd4 = open('%s.tped' %(basename1),'r')
    fd8 = open('%s.tped' %(basename2),'r')

    dn_strand = 'preQC'
    fn_strand = 'HumanOmni2.5-8v1_A-b37.strand'
    fn_strand = 'HumanOmni2.5-8v1_A-b37_exclmissmult_deduplicated_common.strand'
    fp_strand = os.path.join(dn_strand,fn_strand)
    bool_run = True
    if os.path.isfile('%s.sorted' %(fn_strand)):
        if os.path.getmtime('%s.sorted' %(fn_strand)) > 60*60:
            bool_run = False
    if bool_run == True:
        cmd = "cat %s | awk '{" %(fp_strand)
        cmd += 'sub(/X/,23,$2);sub(/Y/,24,$2);sub(/XY/,25,$2);sub(/MT/,26,$2);print $0'
        cmd += "}' | sort -k2,2 -k3,3 -o %s.sorted" %(fn_strand)
        execmd(cmd)
    fdstrand8 = open('%s.sorted' %(fn_strand),'r')

    strand_chr = None
    strand_pos = None
    for line4 in fd4:
        line8 = fd8.readline()
        l4 = line4.split()
        if int(l4[0]) not in range(1,22+1,): continue ## autosomal SNPs only
        l8 = line8.split()
        i_line += 1
        if i_line % 100000 == 0: print i_line
        l_set_alleles = list(set(l4[6:])-set(['0']))
        ## no calls (only if pre SNP QC data)
        if len(l_set_alleles) == 0:
            print l[4]
            print l[8]
            stop
            continue
        if (
            (len(set(l4[6:])-set(['0'])) == 0 and len(set(l8[6:])-set(['0'])) != 0)
            or
            (len(set(l4[6:])-set(['0'])) != 0 and len(set(l8[6:])-set(['0'])) == 0)
            ):
            print l4
            print l8
            stop

##        strand_chr, strand_pos = tmp_check_allele_concordance_bim_strand(
##            l4,strand_chr,strand_pos,l_set_alleles,fdstrand8,)

        if len(l_set_alleles) == 2:
            allele1 = l_set_alleles[0]
            allele2 = l_set_alleles[1]
        else:
            allele1 = l_set_alleles[0]
            allele2 = None
        if l4[1] != l8[1]:
            print line4
            print line8
            stop
##        sum_xy = 0
##        sum_xx = 0
##        sum_yy = 0
##        sum_x = 0
##        sum_y = 0
##        n = 0
        count_diff = 0
        bool_noncall_quad = False
        bool_noncall_octo = False
        l_x = []
        l_y = []
        l_gt = []
        if 2*n_samples*[l4[4]]==l4[4:] or 2*n_samples*[l8[4]]==l8[4:]:
            l_monomorphic += [l4[1]]
        for i4 in xrange(n_samples):
            sample = l_fam4[i4]
            i8 = l_fam8.index(sample)
            genotype4 = l4[4+2*i4:4+2*(i4+1)]
            genotype8 = l8[4+2*i8:4+2*(i8+1)]
            allele4A = genotype4[0]
            allele4B = genotype4[1]
            allele8A = genotype8[0]
            allele8B = genotype8[1]
            if (
                allele4A == '0'# or allele4B == '0'
                or
                allele8A == '0'# or allele8B == '0'
                ):
                ## octo called, quad not called
                if (
                    allele4A == '0' and allele4B == '0'
                    and
                    allele8A != '0' and allele8B != '0'
                    ):
##                    count_quad0octo1 += 1
                    bool_noncall_quad = True
                ## quad called, octo not called
                elif (
                    allele4A != '0' and allele4B != '0'
                    and
                    allele8A == '0' and allele8B == '0'
                    ):
##                    count_quad1octo0 += 1
                    bool_noncall_octo = True
                ## quad not called, octo not called
                elif (
                    allele4A == '0' and allele4B == '0'
                    and
                    allele8A == '0' and allele8B == '0'
                    ):
##                    count_quad0octo0 += 1
                    bool_noncall_quad = True
                    bool_noncall_octo = True
                else:
                    print genotype4, genotype8
                    print sample
                    print l4
                    stop
                continue
            ## genotype difference
##            elif genotype4 != genotype8:
            elif set(genotype4) != set(genotype8):
                count_diff += 1
##                count_quad1octo1 += 1
##            else:
##                count_quad1octo1 += 1
##                pass

            ##
            ## "dosages"
            ##
            d = {}
            for k,alleleA,alleleB in [
                [4,allele4A,allele4B,],
                [8,allele8A,allele8B,],
                ]:
                dosage = 0
                if alleleA == allele1:
                    dosage += 1
                if alleleB == allele1:
                    dosage += 1
                d[k] = dosage
            x = d[4]
            y = d[8]
            sum_xy += x*y
            sum_xx += x*x
            sum_yy += y*y
            sum_x += x
            sum_y += y
            n += 1
            l_x += [x]
            l_y += [y]
            l_gt += [[genotype4,genotype8,]]
            ## continue loop over samples
            continue
        ##
        ## count
        ##
        if bool_noncall_quad == True and bool_noncall_octo == True:
            count_quad0octo0 += 1
        elif bool_noncall_quad == True:
            count_quad0octo1 += 1
        elif bool_noncall_octo == True:
            count_quad1octo0 += 1
        else:
            count_quad1octo1 += 1
        ##
        ## continue if no call
        ##
        if bool_noncall_quad == True or bool_noncall_octo == True:
            ## tmp
            if count_diff > 0 or i4 == 0:
                l_disconcordant += [l4[1]]
            continue
##        if n > 1:
##            denominator = math.sqrt((sum_xx-sum_x**2/n)*(sum_yy-sum_y**2/n))
##            nominator = (sum_xy-sum_x*sum_y/n)
##            if denominator == 0:
##                r = 1.
##            else:
##                r = nominator/denominator
##            r2 = r**2
##            if r == 1 and count_diff == 29 and l4[1][:2] == 'rs':
##                for i4 in xrange(n_samples):
##                    sample = l_fam4[i4]
##                    i8 = l_fam8.index(sample)
##                    genotype4 = l4[4+2*i4:4+2*(i4+1)]
##                    genotype8 = l8[4+2*i8:4+2*(i8+1)]
##                    print genotype4, genotype8
##                print l4[:4]
##                stop
        d_diff['all'][count_diff] += 1
        if count_diff > 0:
            l_disconcordant += [l4[1]]

##        ## all different
##        if count_diff == n_samples:
##            if l4[1][:2] == 'rs':
##                print
##                print line4
##                print line8

        ## monomorphic
        if 2*n_samples*[l4[4]]==l4[4:] or 2*n_samples*[l8[4]]==l8[4:]:
            if count_diff == n_samples:
                l_26diff_and_monomorphic += [l4[1]]
##                if (
##                    l4[1][:2] == 'rs' and n_samples*l4[4:6]==l4[4:] and n_samples*l8[4:6]==l8[4:]
##                    and
##                    2*n_samples*[l4[4]]==l4[4:] and 2*n_samples*[l8[4]]==l8[4:]
##                    ):
##                    print line4
##                    print line8
##                    stop1_both_monomorphic
##                if l4[1][:2] == 'rs' and n_samples*l4[4:6]==l4[4:] and n_samples*l8[4:6]!=l8[4:]:
##                    print line4
##                    print line8
##                    stop2_quad_monomorphic
##                if l4[1][:2] == 'rs' and n_samples*l4[4:6]!=l4[4:] and n_samples*l8[4:6]==l8[4:]:
##                    print line4
##                    print line8
##                    stop3_octo_monomorphic
                pass
            d_diff['monomorphic'][count_diff] +=1
        ## not monomorphic
        else:
            d_diff['nonmonomorphic'][count_diff] +=1

        ## continue loop over lines
        continue

    fd4.close()
    fd8.close()

    fd = open('26diff_and_monomorphic.SNPs','w')
    fd.write('\n'.join(l_26diff_and_monomorphic)+'\n')
    fd.close()

    fd = open('monomorphic.SNPs','w')
    fd.write('\n'.join(l_26diff_and_monomorphic)+'\n')
    fd.close()

    fd = open('disconcordant.SNPs','w')
    fd.write('\n'.join(l_disconcordant)+'\n')
    fd.close()

    denominator = math.sqrt((sum_xx-sum_x**2/n)*(sum_yy-sum_y**2/n))
    nominator = (sum_xy-sum_x*sum_y/n)
    r = nominator/denominator
    r2 = r**2
    print 'n', n, 'r', r, 'r2', r2

    for count_diff in xrange(n_samples+1):
        print count_diff,
        for k in ['all','monomorphic','nonmonomorphic',]:
            print d_diff[k][count_diff],
        print
    for k in ['all','monomorphic','nonmonomorphic',]:
        print k, sum(d_diff[k].values())-d_diff[k][0]
    print '01,10,00,11', count_quad0octo1, count_quad1octo0, count_quad0octo0, count_quad1octo1
    print 'all', i_line
    print 'all', int(os.popen('cat %s.tped | wc -l' %(basename1)).read())
    print 'all', int(os.popen('cat %s.tped | wc -l' %(basename2)).read())
    
    return


def plots(l_populations):

    print 'plots'

    for suffix in [
        ## samples (concatenate)
        'imiss','het','sexcheck','genome',
        ## SNPs (paste/join)
##        'frq','hwe','SNPQC.lmiss',
        'fam','sampleQC.samples',
##        'mds',
        ]:
        
        bool_continue = False
        for population in l_populations:
            if not os.path.isfile('%s.%s' %(population,suffix,)):
                bool_continue = True
                break
            continue
        if bool_continue == True:
            print 'skip', suffix
            continue
        else:
            print 'concatenate', suffix

        fd = open('agv.%s' %(suffix),'w')
        fd.close()

        if not suffix in ['fam','sampleQC.samples',]:
##            cmd = 'head -1 %s.%s > agv.%s' %(l_populations[0],suffix,suffix,)
            if suffix in ['het','sexcheck',]:
                cmd = 'head -1 %s.sequential.%s > agv.%s' %(
                    l_populations[0],suffix,suffix,)
            else:
                cmd = 'head -1 %s.%s > agv.%s' %(l_populations[0],suffix,suffix,)
            execmd(cmd)
        for population in l_populations:
            ## no header
            if suffix in ['fam','sampleQC.samples',]:
                cmd = "cat %s.%s >> agv.%s" %(population,suffix,suffix,)
            ## header
            else:
                cmd = "sed '1d' %s.%s >> agv.%s" %(population,suffix,suffix,)
            execmd(cmd)
            continue

        continue

    sys.argv += ['--bfile','agv','--project','agv',]
    instanceQC = QC.main()
##    instanceQC.plink_plots('agv',i_wait=0)

    ## samples
    instanceQC.histogram_imiss('agv',)
    instanceQC.histogram_het('agv',bool_with_stddev=False,)
##    instanceQC.histogram_genome('agv',)
    instanceQC.scatter_het_call('agv',bool_with_stddev=False,)
    if os.path.isfile('agv.mds'):
        instanceQC.scatter_mds('agv')

##    ## SNPs
##    instanceQC.scatter_lmiss_frq('agv')
##    instanceQC.histogram_lmiss('agv')
##    instanceQC.histogram_frq('agv')
##    instanceQC.histogram_hwe('agv')

    return


def tables(l_populations):

    print 'tables'

    s_sample = ''
    s_SNP = ''
    for population in l_populations:
        print 'table', population
        s_sample += '%s\t' %(population)
        s_SNP += '%s\t' %(population)

        l = []
        for suffix in ['imiss','het','sexcheck','sampleQC','genome.0.10']:
            cmd = 'cat pops/%s/%s.%s.samples | wc -l' %(
                population,population,suffix,)
            l += ['%s' %(os.popen(cmd).read().strip()),]
        s_sample += '\t'.join(l)

        l = []
        cmd = 'cat %s.bim | wc -l' %(population,)
        l += ['%s' %(os.popen(cmd).read().strip()),]
        for suffix in [
            'X','autosomes','lmiss','hwe',
            ]:
            cmd = 'cat pops/%s/%s.%s.SNPs | wc -l' %(
                population,population,suffix,)
            l += ['%s' %(os.popen(cmd).read().strip()),]
        s_SNP += '\t'.join(l)

        s_sample += '\n'
        s_SNP += '\n'

    fd = open('summary_samples.table','w')
    fd.write(s_sample)
    fd.close()

    fd = open('summary_SNPs.table','w')
    fd.write(s_SNP)
    fd.close()

    return


if __name__ == '__main__':
    main()
