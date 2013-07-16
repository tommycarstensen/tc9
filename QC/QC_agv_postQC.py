#!/software/bin/python

## Tommy Carstensen (tc9)
## Wellcome Trust Sanger Institute, 2012-2013

## built-ins
import os, sys, time, inspect
## add-ons
import numpy
from xml.dom import minidom
## my own modules
sys.path.append(os.path.dirname(sys.argv[0]))
import QC
sys.path.append('/nfs/users/nfs_t/tc9/github/sandbox')
import gnuplot

hours = 24000


def main():

    l_populations, d_pops2coords = init()
    l_bfiles = [
        '../preQC/omni2.5-4_20120904_agv_gtu_sexupdated_exclconsent_commonID_build37_renamed_flipped_commonpos',
        '../preQC/omni2.5-8_agv_20120910_gtu_exclconsent_commonID_flipped_commonpos',
        ]

    if False:
        plots(l_populations)
##        sys.exit(0)
        tables(l_populations)
##        sys.exit(0)

    ## 15) analyze discordance
    analyze_discordance(
        l_populations,d_pops2coords,l_bfiles,)

    return


def parse_peakware():

    ## http://peakware.com/googleearth.html?map=area&id=432
    ## http://en.wikipedia.org/wiki/List_of_mountain_ranges#Africa
    ## http://en.wikipedia.org/wiki/List_of_mountains#Africa
    lines_mountains = []
    fn_xml = 'peakware-africa.kml'
    xmldoc = minidom.parse(fn_xml)
    placemarks = xmldoc.getElementsByTagName('Placemark')
    for placemark in placemarks:
        LookAts = placemark.getElementsByTagName("LookAt")
        NodeList_name = placemark.getElementsByTagName("name")
        for node_name in NodeList_name:
            name = node_name.childNodes[0].nodeValue
        for LookAt in LookAts:
            ## loop over NodeList
            for node_longitude in LookAt.getElementsByTagName('longitude'):
                longitude = node_longitude.childNodes[0].nodeValue
            ## loop over NodeList
            for node_latitude in LookAt.getElementsByTagName('latitude'):
                latitude = node_latitude.childNodes[0].nodeValue
##        for description in NodeList_description:
##            print description.getElements
##            for node_li in description.getElementsByTagName('li'):
##                li = node_li.childNodes[0].nodeValue
##                print li
        lines_mountains += ['%s %s %s\n' %(longitude,latitude,name,)]

    fd = open('mountains.txt','w')
    fd.writelines(lines_mountains)
    fd.close()

    return lines_mountains


def write_eigensoft_parameter_file(bfile_out):

    ##
    ## input
    par = 'genotypename: %s.bed\n' %(bfile_out)
    par += 'snpname: %s.bim\n' %(bfile_out)
    par += 'indivname: %s.fam\n' %(bfile_out)
    ##
    ## output
    par += 'evecoutname: %s.evec\n' %(bfile_out)
    par += 'evaloutname: %s.eval\n' %(bfile_out)
    ##
    ## optional parameters
    ##
    ## numoutevec:     number of eigenvectors to output.  Default is 10.
    par += 'numoutevec: 10\n'
    ##
    ## maximum number of outlier removal iterations.
    ## Default is 5.  To turn off outlier removal, set this parameter to 0.
    par += 'numoutlieriter: 0\n'
    ##
    ## snpweightoutname: output file containing SNP weightings of each
    ## principal component.  Note that this output file does not contain entries
    ## for monomorphic SNPs from the input .snp file.
    par += 'snpweightoutname: %s.snpweight\n' %(bfile_out)
    ##
    ## http://helix.nih.gov/Applications/README.eigenstrat
    ## -q YES/NO        :
    ## If set to YES, assume that there is a single population and
    ## the population field contains real-valued phenotypes.
    ## (Corresponds to qtmode parameter in smartpca program.)
    ## The default value for this parameter is NO.
    par += 'qtmode: 0\n'

    fd = open('%s.par' %(bfile_out),'w')
    fd.write(par)
    fd.close()

    return


def run_eigensoft(population,QC,):

    bfile_out = '%s_%s' %(population,QC,)

    if os.path.isfile('%s.par' %(bfile_out,)):
        return
    if os.path.isfile('%s.snpweight' %(bfile_out,)):
        return

    l_cmds = []

    ##
    ## 1) rename samples for EIGENSOFT
    ##
    l_basenames = []
    for chip in ['quad','octo',]:
        bfile_in = '../pops/%s_%s/%s_%s.%s' %(population,chip,population,chip,QC,)
        basename = os.path.basename(bfile_in)
        l_basenames += [basename]
        cmd = 'cat %s.fam' %(bfile_in)
        cmd += ' | '
##        cmd += "awk '{print $1,$2,substr($1,12,10),substr($2,12,10)}'"
        cmd += "awk '{print $1,$2,substr($1,length($2)-9,10),substr($2,length($2)-9,10)}'"
        cmd += ' > %s_%s.recoded.txt' %(basename,chip)
        l_cmds += [cmd]
    l_cmds += ['cat %s_quad.recoded.txt %s_octo.recoded.txt > %s.recoded.txt' %(
        l_basenames[0],l_basenames[1],bfile_out,)]

    ##
    ## 2) find common SNPs
    ##
    fn_common_SNPs = '%s.comm.SNPs' %(bfile_out,)
    l_basenames = []
    for chip in ['quad','octo',]:
        bfile_in = '../pops/%s_%s/%s_%s.%s' %(population,chip,population,chip,QC,)
        basename = os.path.basename(bfile_in)
        l_basenames += [basename]
        cmd = "cat %s.bim | awk '{print $2}' | sort > %s.SNPs.sorted" %(
            bfile_in,basename,)
        l_cmds += [cmd]
    cmd = 'comm -12 %s.SNPs.sorted %s.SNPs.sorted > %s' %(
        l_basenames[0],l_basenames[1],fn_common_SNPs,)
    l_cmds += [cmd]

    ##
    ## 3) merge and extract common SNPs and rename samples for EIGENSOFT
    ##
    cmd = 'plink \\\n'
    cmd += '--bfile ../pops/%s_quad/%s_quad.%s \\\n' %(population,population,QC)
    cmd += '--bmerge \\\n'
    cmd += '../pops/%s_octo/%s_octo.%s.bed \\\n' %(population,population,QC)
    cmd += '../pops/%s_octo/%s_octo.%s.bim \\\n' %(population,population,QC)
    cmd += '../pops/%s_octo/%s_octo.%s.fam \\\n' %(population,population,QC)
    if population == 'Baganda':
        execmd('cat Baganda_quad.%s.comm.samples Baganda_octo.%s.comm.samples > Baganda.%s.comm.samples' %(QC,QC,QC))
        cmd += '--remove Baganda.%s.comm.samples \\\n' %(QC)
    cmd += '--extract %s \\\n' %(fn_common_SNPs)
    cmd += '--update-ids %s.recoded.txt \\\n' %(bfile_out,)
    cmd += '--make-bed --out %s \\\n' %(bfile_out)
    l_cmds += [cmd]

##    path_EIG = '/nfs/team149/Software/EIG4.2/bin'
####    path_EIG = '/software/varinf/bin/eigensoft/bin'
####    path_EIG = '/nfs/team149/Software/usr/share/eigensoft/EIG4.2/bin'
##    cmd = '%s/smartpca.perl \\\n' %(path_EIG)
##    ## genotype file in any format (see ../CONVERTF/README)
##    cmd += '-i %s.longLDexcl.bed \\\n' %(os.path.join(os.getcwd(),out_prefix))
##    ## snp file in any format (see ../CONVERTF/README)
##    cmd += '-a %s.longLDexcl.bim \\\n' %(os.path.join(os.getcwd(),out_prefix))
##    ## indiv file in any format (see ../CONVERTF/README)
##    cmd += '-b %s.longLDexcl.fam \\\n' %(os.path.join(os.getcwd(),out_prefix))
##    ## (Default is 10) number of principal components to output
##    cmd += '-k 10 \\\n'
##    ## output file of principal components.  Individuals removed
##    ## as outliers will have all values set to 0.0 in this file.
##    cmd += '-o %s.pca \\\n' %(out_prefix.replace('.','_'))
##    ## prefix of output plot files of top 2 principal components.
##    ## (labeling individuals according to labels in indiv file)
##    cmd += '-p %s.plot \\\n' %(out_prefix.replace('.','_'))
##    ## output file of all eigenvalues
##    cmd += '-e %s.eval \\\n' %(out_prefix.replace('.','_'))
##    ## output logfile        
##    cmd += '-l %s.eigensoft.log \\\n' %(out_prefix.replace('.','_'))
##    ## (Default is 5) maximum number of outlier removal iterations.
##    ## To turn off outlier removal, set -m 0.
##    cmd += '-m 5 \\\n'
##    ## (Default is 10) number of principal components along which
##    ## to remove outliers during each outlier removal iteration.
##    cmd += '-t 10 \\\n'
##    ## (Default is 6.0) number of standard deviations which an
##    ## individual must exceed, along one of topk top principal
##    ## components, in order to be removed as an outlier.
##    cmd += '-s 6 \\\n'
####    ## verbose
####    cmd += '-V \\\n'
##    cmd += '-snpweightoutname %s.snpweight \\\n' %(out_prefix)
##    cmd += '-badsnpname %s.badsnpname \\\n' %(out_prefix)
####    cmd += ' -p Muganda29_octo200_quad29.pca.par'

    ##
    ## parameter file
    ## http://computing.bio.cam.ac.uk/local/doc/popgen.txt
    write_eigensoft_parameter_file(bfile_out)

    cmd = '/nfs/team149/Software/EIG4.2/bin/smartpca -p %s.par' %(bfile_out)
##    cmd = '/software/varinf/bin/EIG3.0/bin/smartpca -p %s.par' %(bfile_out)
##    execmd(cmd)
    l_cmds += [cmd]

    fn_sh = '%s_eigensoft.sh' %(bfile_out,)
    write_shell(fn_sh,l_cmds,)

    cmd = LSF_append('./%s' %(fn_sh),mem=4000,)
    execmd(cmd)

    return


def allelic_concordance_EIGENSOFT():

    import numpy
    from scipy import stats

    QC = 'SNPQC'

    ##
    ## merge quad/octo
    ##
    for population in ['Zulu','Banyarwanda','Baganda','Ga-Adangbe',]:
        print population
        fn_ld_regions = 'ldregions.SNPs'
        run_eigensoft(population,QC,)

    for population in ['Baganda','Zulu','Banyarwanda','Ga-Adangbe',]:
        print population
        bfile = '%s_%s' %(population,QC,)
        fn_snpweight = '%s.snpweight' %(bfile,)

        if not os.path.isfile(fn_snpweight):
            print 'missing', fn_snpweight
            sys.exit(0)

        ## sort by rsID
        cmd = "cat %s.snpweight | awk '{print $1,$4}' | sort -k1,1" %(bfile)
        cmd += ' > %s.snpweight.sorted' %(bfile)
        execmd(cmd)
        ## all SNPs
        cmd = "cat %s.snpweight.sorted | awk '{print $1}' > %s.EIGENSOFT.SNPs" %(bfile,bfile,)
        execmd(cmd)
        ## all weights
        cmd = "cat %s.snpweight.sorted | awk '{print $2}' > %s.snpweight.all" %(bfile,bfile,)
        execmd(cmd)
        for threshold in range(10):
            ## "heavy" SNPs
            cmd = "cat %s.snpweight.sorted | awk '{if($2>=%i||$2<=-%i) print $1}' > %s_EIGENSOFT_heavy_%i.SNPs" %(bfile,threshold,threshold,bfile,threshold,)
            execmd(cmd)
        ## discordant SNPs
        cmd = 'comm -12 %s.EIGENSOFT.SNPs discordant_sampleQC.SNPs.sorted > %s.discordant.SNPs' %(bfile,bfile,)
        execmd(cmd)
        ## concordant SNPs
        cmd = 'comm -23 %s.EIGENSOFT.SNPs discordant_sampleQC.SNPs.sorted > %s.concordant.SNPs' %(bfile,bfile,)
        execmd(cmd)
        ## discordant weights
        cmd = 'join %s.discordant.SNPs %s.snpweight.sorted -o 2.2 > %s.snpweight.discordant' %(bfile,bfile,bfile,)
        execmd(cmd)
        ## concordant weights
        cmd = 'join %s.concordant.SNPs %s.snpweight.sorted -o 2.2 > %s.snpweight.concordant' %(bfile,bfile,bfile,)
        execmd(cmd)

    for population in ['Zulu','Banyarwanda','Baganda','Ga-Adangbe',]:
        bfile = '%s_%s' %(population,QC,)
        for k in ['all','concordant','discordant',]:
            fn = '%s.snpweight.%s' %(bfile,k,)
##            fd = open(fn,'r')
##            s = fd.read()
##            fd.close()
####            cmd = "cat %s"  %(fn,)
####            s = os.popen(cmd).read()
##            array = numpy.fromstring(s,sep='\n')
##            print k, stats.describe(array)
            n_samples = int(os.popen('cat %s.fam | wc -l' %(bfile)).read())
            n_SNPs = int(os.popen('cat %s | wc -l' %(fn)).read())
            gnuplot.histogram2(
                fn,
                prefix_out = 'eig.%s' %(fn),
                x_step = 0.5,
##                x_min = x_min, x_max = x_max,
                xlabel='EIGENSOFT SNP weight',
                title='%s (n_{samples}=%i, n_{SNPs}=%i)' %(
                    population,n_samples,n_SNPs,
                    ),
                )

    for threshold in range(10):
        fn1='Baganda_SNPQC_EIGENSOFT_heavy_%i.SNPs' %(threshold)
        fn2='Banyarwanda_SNPQC_EIGENSOFT_heavy_%i.SNPs' %(threshold)
        fn3='Zulu_SNPQC_EIGENSOFT_heavy_%i.SNPs' %(threshold)
        fn4='Ga-Adangbe_SNPQC_EIGENSOFT_heavy_%i.SNPs' %(threshold)
        text1='Baganda'
        text2='Banyarwanda'
        text3='Zulu'
        text4='Ga-Adangbe'
        gnuplot.venn4(
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            fn1=fn1,
            fn2=fn2,
            fn3=fn3,
            fn4=fn4,
            text1=text1,
            text2=text2,
            text3=text3,
            text4=text4,
            suffix='agv_EIGENSOFT_monomorphic_%i' %(threshold),
            bool_remove = False,
            )

        execmd('comm -12 11xx xx11 > EIGENSOFT_heavy_%i_intersection.SNPs' %(threshold))

        execmd('comm -12 Baganda_SNPQC_EIGENSOFT_heavy_%i.SNPs Banyarwanda_SNPQC_EIGENSOFT_heavy_%i.SNPs > EIGENSOFT_heavy_%i_intersectBanBag.SNPs' %(threshold,threshold,threshold))

        execmd('cat %s %s %s %s | sort -u > EIGENSOFT_heavy_%i_union.SNPs' %(fn1,fn2,fn3,fn4,threshold))

        for pop1 in ['Zulu','Banyarwanda','Baganda','Ga-Adangbe',]:
            fn_out = '%s_SNPQC_EIGENSOFT_heavy_%i_union3.SNPs' %(pop1,threshold,)
            cmd = 'cat'
            for pop2 in ['Zulu','Banyarwanda','Baganda','Ga-Adangbe',]:
                fn_in = '%s_SNPQC_EIGENSOFT_heavy_%i.SNPs' %(pop2,threshold,)
                if pop1 == pop2:
                    continue
                cmd += ' %s' %(fn_in)
            cmd += ' | sort -u > %s' %(fn_out)
            execmd(cmd)

    return


def frqlmisshwe_discordant():

    fn_discordant_SNPs = 'discordant_sampleQC.SNPs'

    for chip in ['quad','octo',]:

        cmd = 'cat ../pops/Baganda_%s/Baganda_%s.postQC.autosomes.fam | wc -l' %(chip,chip,)
        n_samples = int(os.popen(cmd).read())

        for flag,prefix in [
            ['','discordant',],
            ['-v','concordant',],
            ]:

            for suffix,column,x_min,x_max,x_step,xlabel,condition in [
                ['SNPQC.frq','$5',0,0.5,0.01,'MAF','',],
                ['SNPQC.lmiss','1-$5',0.90,1.00,0.01,'SNP call rate','',],
                ['hwe','-log($9)/log(10)',0,10,1,'p_{HWE}','if(NR%3==2)',],
                ]:

                cmd = 'cat ../pops/Baganda_%s/Baganda_%s.%s' %(chip,chip,suffix,)
                cmd += " | awk '{%s print}'" %(condition,)
                cmd += ' | fgrep %s -w -f %s' %(flag,fn_discordant_SNPs)
                cmd += " | awk '{print %s}'" %(column,)
                cmd += ' > %s.%s.%s' %(suffix,chip,prefix,)
                execmd(cmd)

                cmd = 'cat %s.%s.%s | wc -l' %(suffix,chip,prefix,)
                n_SNPs = int(os.popen(cmd).read())
                if flag == '-v':
                    n_SNPs -= 1

                gnuplot.histogram2(
                    '%s.%s.%s' %(suffix,chip,prefix,),
                    x_step = x_step,
                    x_min = x_min, x_max = x_max,
                    xlabel=xlabel,
                    title='Baganda %s (n_{samples}=%i, n_{SNPs}=%i)' %(
                        chip,n_samples,n_SNPs,
                        ),
                    )

                os.remove('%s.%s.%s' %(suffix,chip,prefix,))

    return


def Baganda_frq_scatter(QC):

    if os.path.isfile('frq.%s.incldiscordant.joined.png' %(QC)):
        return

    frq1 = '../pops/Baganda_quad/Baganda_quad.%s.frq' %(QC)
    frq2 = '../pops/Baganda_octo/Baganda_octo.%s.frq' %(QC)

    basename1 = os.path.basename(frq1)
    basename2 = os.path.basename(frq2)

    for frq in [frq1,frq2,]:
        cmd = 'cat %s | sort -k2,2 > %s.sorted' %(frq,os.path.basename(frq),)
        execmd(cmd)

    cmd = 'join -1 2 -2 2 -o 0,1.3,1.4,1.5,2.3,2.4,2.5'
    cmd += ' %s.sorted %s.sorted' %(basename1,basename2,)
    cmd += " | awk '{if($3==$6) {print $4,$7,$1} else {print $4,1-$7,$1}}'"
    cmd += ' > frq.incldiscordant.joined'
    execmd(cmd)

    cmd = 'fgrep -w -v -f discordant_%s.SNPs frq.incldiscordant.joined > frq.excldiscordant.joined' %(QC)
    execmd(cmd)

    for s in ['incldiscordant','excldiscordant',]:

##        if os.path.isfile('frq.%s.%s.joined.png' %(QC,s)):
##            continue

    ##    n_samples = int(os.popen('cat %s.fam | wc -l').read())-1
        n_SNPs = int(os.popen('cat frq.%s.joined | wc -l' %(s)).read())-1

        cmd = 'cat frq.%s.joined' %(s)
        cmd += " | awk '{if($1!=0&&$2!=0&&$1!=0.5&&$2!=1&&(($1-$2)>0.2||($2-$1)>0.2)) print}'"
        cmd += ' > frq.%s.joined.labels' %(s)
        execmd(cmd)

        cmd = 'cat frq.%s.joined' %(s)
        cmd += " | awk '{if((($1-$2)>0.2||($2-$1)>0.2)) print}'"
        cmd += ' > frq.%s.joined.outliers' %(s)
        execmd(cmd)

        line_plot = 'plot [0:0.5][0:1]'
        line_plot += '"frq.%s.joined" u 1:2 pt 7 lc 0 t ""' %(s)
        if os.path.getsize('frq.%s.joined.outliers' %(s)) > 0:
            line_plot += ',"frq.%s.joined.labels" u 1:($2+0.005):3 w labels t "" font "Verdana,20"' %(s)
            line_plot += ',"frq.%s.joined.outliers" u 1:2 pt 7 lc 1 t ""' %(s)

        if s == 'incldiscordant':
            title_concordance = 'incl discordant SNPs'
        else:
            title_concordance = 'excl discordant SNPs'

        if QC == 'SNPQC':
            title_QC = 'post SNP QC'
        else:
            title_QC = 'pre SNP QC'

        gnuplot.scatter_plot_2d(
            'frq.%s.%s.joined' %(QC,s),
            line_plot = line_plot,
    ##        line_plot = line_plot,
    ##        column1 = 1, column2 = 2,
            xmin=0,xmax=0.5,ymin=0,ymax=1,
            xlabel = 'MAF quad',
            ylabel = 'MAF octo',
            title='MAF, %s, %s (n_{SNPs}=%i)' %(
                title_concordance, title_QC, n_SNPs,
                ),
    ##        prefix_out='%s.mds' %(prefix),
    ##        lines_extra=lines_extra,
            bool_remove=False,
            )

        os.remove('frq.%s.joined.outliers' %(s))
        os.remove('frq.%s.joined.labels' %(s))
        os.remove('frq.%s.joined' %(s))
    
    return


def write_shell(fn_sh,l_cmds):

    fd = open(fn_sh, 'w')
    fd.write('\n'.join(l_cmds))
    fd.close()

    execmd('chmod +x %s' %(fn_sh))

    return


def LSF_append(cmd,job_name=None,mem=3000,queue='normal',):

    cmd_lsf = 'bsub -G agv'
    cmd_lsf += ' -q %s' %(queue)
    cmd_lsf += " -M%i000 -R'select[mem>%i] rusage[mem=%i]'" %(mem,mem,mem,)
    if job_name:
        cmd_lsf += ' -J"%s"' %(job_name)
    cmd_lsf += ' %s' %(cmd)

    return cmd_lsf


def MDS(population,QC,bool_excl_discordant,bool_excl_MAF_below_5,set_name,threshold,prefix,):

    if True:
        for suffix in ['mds','genome','prune.in','bed','log',]:
            if os.path.isfile('%s.%s' %(prefix,suffix,)):
                return

    if population in ['Baganda']:
        bool_merge = True
    elif population in ['Ga-Adangbe','Zulu','Banyarwanda',]:
        bool_merge = False
        bool_merge = True
    else:
        print population
        bool_merge = False

    fn_ld_regions = 'ldregions.SNPs'

    l_cmds = []
    ##
    ## find common SNPs
    ##
    if bool_merge == True:
        l_basenames = []
        for chip in ['quad','octo',]:
            bfile = '../pops/%s_%s/%s_%s.%s' %(population,chip,population,chip,QC,)
            basename = os.path.basename(bfile)
            l_basenames += [basename]
            cmd = "cat %s.bim | awk '{print $2}' | sort > %s.SNPs" %(bfile,basename,)
            l_cmds += [cmd]
        fn_common_SNPs = '%s.comm.SNPs' %(prefix,)
        cmd = 'comm -12 %s.SNPs %s.SNPs > %s' %(
            l_basenames[0],l_basenames[1],fn_common_SNPs,)
        l_cmds += [cmd]
##    fn_common_SNPs = 'postQC.incl.SNPs.intersection'

    ## cat ac20_QC98/Zulu.snpweights.txt | awk '{if($4<-4) print $1}' | sort > Zulu_EIGENSOFT_heavy.SNPs
    
    ##
    ## --bmerge
    ##
##    if not os.path.isfile('%s.bed' %(prefix)):
    if True:
        cmd = 'plink \\\n'
        if bool_merge == True:
            cmd += '--bfile ../pops/%s_quad/%s_quad.%s \\\n' %(population,population,QC)
            cmd += '--bmerge \\\n'
            cmd += '../pops/%s_octo/%s_octo.%s.bed \\\n' %(population,population,QC)
            cmd += '../pops/%s_octo/%s_octo.%s.bim \\\n' %(population,population,QC)
            cmd += '../pops/%s_octo/%s_octo.%s.fam \\\n' %(population,population,QC)
            if population == 'Baganda':
                cmd += '--remove Baganda.%s.comm.samples \\\n' %(QC)
            cmd += '--extract %s \\\n' %(fn_common_SNPs)
        elif bool_merge == False:
            stop
            cmd += '--bfile ../pops/%s/%s.%s \\\n' %(population,population,QC,)
##            cmd += '--extract Baganda.%s.comm.SNPs \\\n' %(QC) ## ms23, tc9
##            cmd += '--extract %s \\\n' %(fn_common_SNPs) ## 8intersect
        if bool_excl_discordant == True:
            if set_name == 'intersection':
                fn_exclude = 'EIGENSOFT_heavy_%i_intersection.SNPs' %(threshold) ## 3b) EIGENSOFT
            elif set_name == 'union3pops':
                fn_exclude = '%s_%s_EIGENSOFT_heavy_%i_union3.SNPs' %(population,QC,threshold,) ## 3a) EIGENSOFT
            elif set_name == 'intersect_and_postSNPQCdisc':
                fn_exclude = 'intersect%i_and_postSNPQCdisc.SNPs' %(threshold)
                cmd_sub = 'cat EIGENSOFT_heavy_%i_intersection.SNPs discordant_SNPQC.SNPs' %(threshold)
                cmd_sub += ' | sort -u > %s' %(fn_exclude)
                execmd(cmd_sub)
            elif set_name == 'intersect_and_postSNPQCdisc_and_HWE':
                fn_exclude = 'intersect%i_and_postSNPQCdisc.SNPs' %(threshold)
                cmd_sub = 'cat EIGENSOFT_heavy_%i_intersection.SNPs discordant_SNPQC.SNPs %s_SNPQC.hwe.SNPs' %(threshold,population,)
                cmd_sub += ' | sort -u > %s' %(fn_exclude)
                execmd(cmd_sub)
            elif set_name == 'intersectBanBag':
                fn_exclude = 'EIGENSOFT_heavy_%i_intersectBanBag.SNPs' %(threshold) ## 3b) EIGENSOFT
            elif set_name == 'union':
                fn_exclude = 'EIGENSOFT_heavy_%i_union.SNPs' %(threshold) ## 3c) EIGENSOFT
            else:
                fn_exclude = '%s_%s_EIGENSOFT_heavy_%i.SNPs' %(population,QC,threshold,) ## 3a) EIGENSOFT
            cmd += '--exclude %s \\\n' %(fn_exclude)
##            cmd += '--exclude discordant_%s.SNPs \\\n' %(QC) ## 2) postQC ms23
##            cmd += '--exclude discordant_sampleQC.SNPs \\\n' ## 1) preQC tc9
        cmd += '--make-bed --out %s \\\n' %(prefix)
##        execmd(cmd)
        l_cmds += [cmd]

    ##
    ## --indep-pairwise
    ##
##    if not os.path.isfile('%s.prune.in' %(prefix)):
    if True:
        cmd = 'plink \\\n'
        cmd += '--bfile %s \\\n' %(prefix)
        cmd += '--out %s \\\n' %(prefix)
        ## settings
        cmd += '--indep-pairwise 50 5 0.2 \\\n'
        if bool_excl_MAF_below_5 == True:
            cmd += '--maf 0.05 \\\n'
        ## SNP exclusion
        cmd += '--exclude %s \\\n' %(fn_ld_regions)
##        execmd(cmd)
        l_cmds += [cmd]

    ##
    ## --genome
    ##
##    if not os.path.isfile('%s.genome' %(prefix)):    
    if True:
        cmd = 'plink \\\n'
        cmd += '--bfile %s \\\n' %(prefix)
        cmd += '--out %s \\\n' %(prefix)
        cmd += '--genome \\\n'
        ## SNP exclusion
        cmd += '--extract %s.prune.in \\\n' %(prefix)
        cmd += '--exclude %s \\\n' %(fn_ld_regions)
##        execmd(cmd)
        l_cmds += [cmd]

    ##
    ## --cluster
    ##
##    if not os.path.isfile('%s.mds' %(prefix)):
    if True:
        cmd = 'plink \\\n'
        cmd += '--bfile %s \\\n' %(prefix)
        cmd += '--out %s \\\n' %(prefix)
        cmd += '--cluster \\\n'
        cmd += '--mds-plot 4 \\\n'
        cmd += '--read-genome %s.genome \\\n' %(prefix)
        ## SNP exclusion
        cmd += '--extract %s.prune.in \\\n' %(prefix)
        cmd += '--exclude %s \\\n' %(fn_ld_regions)
##        execmd(cmd)
        l_cmds += [cmd]

    ##
    ## run PLINK
    ##
    fn_sh = '%s_mds.sh' %(prefix,)

    write_shell(fn_sh,l_cmds,)

    cmd = LSF_append('./%s' %(fn_sh))
    execmd(cmd)

    return


def analyze_discordance(
    l_populations,d_pops2coords,l_bfiles,):

##    ##
##    ## find discordant SNPs with EIGENSOFT
##    ##
##    allelic_concordance_EIGENSOFT()

##    ##
##    ## find discordant SNPs with Baganda29
##    ##
##
##    ## 15a) identify discordant SNPs pre SNP QC
##    n_samples_sampleQC, fn_common_SNPs_sampleQC = allelic_concordance_Baganda(
##        '../pops/Baganda_quad/Baganda_quad.sampleQC',
##        '../pops/Baganda_octo/Baganda_octo.sampleQC',
##        )
##
##    ## 15b) identify discordant SNPs post SNP QC
##    n_samples_SNPQC, fn_common_SNPs_SNPQC = allelic_concordance_Baganda(
##        '../pops/Baganda_quad/Baganda_quad.postQC.autosomes',
##        '../pops/Baganda_octo/Baganda_octo.postQC.autosomes',
##        )

##    ##
##    ## calculate HWE for merged datasets
##    ##
##    HWE_merged()
##    stop

    ##
    ## PCA for Africa and all
    ##
    super_MDS(d_pops2coords,l_bfiles,l_populations,)

    ##
    ## intersection between EIGENSOFT discordant and Baganda discordant
    ##
    ##
    for plink_suffix in ['hwe','lmiss',]:
        cmd = 'cat'
        for chip in ['quad','octo',]:
            cmd += ' ../pops/Baganda_%s/Baganda_%s.%s.SNPs' %(chip,chip,plink_suffix)
        cmd += ' | sort -u > Baganda.union.%s.SNPs' %(plink_suffix)
        execmd(cmd)
        cmd = 'comm -12 discordant_sampleQC.SNPs.sorted Baganda.union.%s.SNPs > Baganda.union.%s.discordant.SNPs' %(plink_suffix,plink_suffix,)
        execmd(cmd)
    ##
    d_samples = {'sampleQC':29,'SNPQC':26,} ## tmp!!!
    l_fn = []
    for pop in ['Baganda','Banyarwanda','Zulu','Ga-Adangbe',]:
        l_fn += ['%s_SNPQC_EIGENSOFT_heavy_4.SNPs' %(pop)]
    l_fn += ['EIGENSOFT_heavy_4_intersection.SNPs']
    l_fn += ['EIGENSOFT_heavy_4_union.SNPs']
    for QC in ['sampleQC','SNPQC',]:
        s = ''
        n_samples = d_samples[QC]
        ## append rows
        for i in xrange(0,n_samples+1,):
            fn1 = 'discordant_%s_%i.SNPs' %(QC,i,)
            if i == 0:
                execmd('comm -23 Baganda.%s.comm.SNPs discordant_%s.SNPs.sorted > discordant_%s_0.SNPs.sorted' %(QC,QC,QC,))
            else:
                execmd('cat %s | sort > %s.sorted' %(fn1,fn1))
            ## append columns to row
            for fn2 in l_fn:
                cmd = 'comm -12 %s.sorted %s | wc -l' %(fn1,fn2,)
                i_comm = int(os.popen(cmd).read())
                s += '\t%i' %(i_comm)
            s += '\n'
##        ## append last row
##        s += 'Baganda QCed and Baganda29 discordant \t\t\t\t\t\t\t\t'
##        ## append columns to last row
##        for fn2 in l_fn:
##            cmd = 'comm -23 %s Baganda.SNPQC.comm.SNPs | wc -l' %(fn2)
##            s += '\t%i' %(int(os.popen(cmd).read()))
##        s += '\n'
        for plink_suffix in ['lmiss','hwe',]:
            ## append last row
            s += 'Baganda QCed (%s) and Baganda29 discordant\t\t\t\t\t\t\t' %(plink_suffix)
            ## append columns to last row
            for fn2 in l_fn:
                cmd = 'comm -12 %s Baganda.union.%s.SNPs | wc -l' %(fn2,plink_suffix)
                s += '\t%i' %(int(os.popen(cmd).read()))
            s += '\n'
        ## write table
        fd = open('table%s.append.txt' %(QC),'w')
        fd.write(s)
        fd.close()
        cmd = 'paste table%s.txt table%s.append.txt' %(QC,QC,)
        cmd += ' > table%s.extended.txt' %(QC,)
        execmd(cmd)
        execmd('unix2dos table%s.extended.txt' %(QC))

####    ##
####    ## non-founders
####    ##
####    for chip in ['quad','octo',]:
######        cmd = 'cat Baganda_%s.sampleQC.tped' %(chip)
######        cmd += "| awk '{if($1>=1&&$1<=22) print}'
######        cmd = ' | cut -d " " -f 2,5- '
######        cmd += " | awk '{"
######        cmd += 'rsID=$1; gsub(/0/,""); if(NF==1) print rsID'
######        cmd += "}' | sort > Baganda_call0%s.SNPs" %(chip)
####        cmd = 'cut -d " " -f 2,5- Baganda_%s.sampleQC.tped' %(chip)
####        cmd += " | awk '{"
####        cmd += 'rsID=$1; gsub(/0/,""); if(NF==1) print rsID'
####        cmd += "}' | sort > Baganda_call0%s.SNPs" %(chip)
####        execmd(cmd)
####
####        i = int(os.popen('cat Baganda_call0%s.SNPs | wc -l' %(chip)).read())
####        print 'not called for any samples', chip, i
####    print 'overlap', int(os.popen('comm -12 Baganda_call0quad.SNPs Baganda_call0octo.SNPs | wc -l').read())

##    ##
####    ## either
####    ## 1) find intersection between quad union and octo union
####    ## or
##    ## 2) find intersection between quad and octo for 4 split pops
##    ##
##    cmd = 'cat ../pops/Baganda_quad/Baganda_quad.sampleQC.bim'
##    cmd += " | awk '{print $2}'"
##    cmd += ' | sort'
##    cmd += ' > autosomal.SNPs.sorted'
##    execmd(cmd)
######    ##
######    ## 1) find intersection between quad union and octo union
######    ##
########    for chip in ['quad','octo',]:
########        execmd('cp autosomal.SNPs.sorted %s.postQC.excl.SNPs' %(chip,))
########    for bfile,chip in d_bfiles2chips.items():
########        cmd = 'sort ../pops/%s/%s.SNPQC.SNPs > postQC.SNPs.sorted' %(bfile,bfile,)
########        execmd(cmd)
########        cmd = 'comm -12 postQC.SNPs.sorted %s.postQC.excl.SNPs > comm' %(chip,)
########        execmd(cmd)
########        execmd('mv comm %s.postQC.excl.SNPs' %(chip,))
########        print int(os.popen('cat %s.postQC.excl.SNPs | wc -l' %(chip)).read())
########    cmd = 'cat quad.postQC.excl.SNPs octo.postQC.excl.SNPs | sort -u > postQC.excl.SNPs.union'
########    execmd(cmd)
########    cmd = 'comm -23 autosomal.SNPs.sorted postQC.excl.SNPs.union > postQC.incl.SNPs.intersection'
########    execmd(cmd)
########    ##
########    ## 2) find intersection between quad/octo Baganda, Banyarwanda, Zulu, Ga-Adangbe
########    ##
######    if os.path.isfile('postQC.excl.SNPs.union'):
######        os.remove('postQC.excl.SNPs.union')
######    for pop in ['Baganda','Banyarwanda','Zulu','Ga-Adangbe',]:
######        for chip in ['quad','octo',]:
######            bfile = '%s_%s' %(pop,chip,)
######            cmd = 'cat ../pops/%s/%s.SNPQC.SNPs >> postQC.excl.SNPs.union' %(bfile,bfile,)
######            execmd(cmd)
######    cmd = 'sort -u postQC.excl.SNPs.union -o postQC.excl.SNPs.union'
######    execmd(cmd)
######    cmd = 'comm -23 autosomal.SNPs.sorted postQC.excl.SNPs.union > postQC.incl.SNPs.intersection'
######    execmd(cmd)
##    ##
##    ## 3) find intersection between quad union and octo union
##    ##
##    for chip in ['quad','octo',]:
##        execmd('cp autosomal.SNPs.sorted %s.postQC.excl.SNPs' %(chip,))
##    for bfile,chip in d_bfiles2chips.items():
##        cmd = 'sort ../pops/%s/%s.SNPQC.SNPs > postQC.SNPs.sorted' %(bfile,bfile,)
##        execmd(cmd)
##        cmd = 'comm -12 postQC.SNPs.sorted %s.postQC.excl.SNPs > comm' %(chip,)
##        execmd(cmd)
##        execmd('mv comm %s.postQC.excl.SNPs' %(chip,))
##        print int(os.popen('cat %s.postQC.excl.SNPs | wc -l' %(chip)).read())
##    cmd = 'cat quad.postQC.excl.SNPs octo.postQC.excl.SNPs | sort -u > postQC.excl.SNPs.union'
##    execmd(cmd)
##    cmd = 'comm -23 autosomal.SNPs.sorted postQC.excl.SNPs.union > postQC.incl.SNPs.intersection'
##    execmd(cmd)

##    ##
##    ## table
##    ##
##    s = '''Baganda_quad Baganda_octo Kikuyu Mandinka Zulu_quad Zulu_octo Kalenjin Fula Barundi Banyarwanda_quad Banyarwanda_octo Sotho Wolof Jola Ga-Adangbe_quad Ga-Adangbe_octo Ethiopia YRI TSI PUR PEL MXL LWK KHV JPT IBS IBO GIH GBR FIN CLM CHS CHB CEU CDX ASW ACB'''
##    l_bfiles = s.split(' ')
##    ##
##    ## sort bim files
##    ##
##    for bfile in l_bfiles:
##        if os.path.isfile('%s.SNPQC.bim.sorted' %(bfile)):
##            if os.path.getsize('%s.SNPQC.bim.sorted' %(bfile)):
##                continue
##        cmd = 'cat ../pops/%s/%s.SNPQC.bim' %(bfile,bfile,)
##        cmd += " | awk '{print $2}'"
##        cmd += ' | sort > %s.SNPQC.bim.sorted' %(bfile,)
##        execmd(cmd)
####        newpid = os.fork()
####        if newpid == 0:
####            execmd(cmd)
####            os._exit(0)
##    i_SNPs_autosomal_preQC = int(os.popen('cat ../pops/%s/%s.sampleQC.bim | wc -l' %(l_bfiles[0],l_bfiles[0],)).read())
##    execmd('sort discordant_sampleQC.SNPs > discordant_sampleQC.SNPs.sorted')
##    execmd('sort discordant_SNPQC.SNPs > discordant_SNPQC.SNPs.sorted')
##    ## check that child sorting processes have finished
##    for bfile in l_bfiles:
##        while not os.path.isfile('%s.SNPQC.bim.sorted' %(bfile)):
##            print 'isfile', bfile
##            time.sleep(10)
##        while not os.path.getsize('%s.SNPQC.bim.sorted' %(bfile)):
##            print 'getsize', bfile
##            time.sleep(10)
##    ## initiate table
##    table = ''
##    ## table headers
##    table += 'a) samples before QC\n'
##    table += 'b) autosomal SNPs before QC\n'
##    table += 'c) samples after QC\n'
##    table += 'd) autosomal SNPs after QC\n'
##    table += 'e) autosomal SNPs excluded during QC (d-b)\n'
##    table += u"f) intersection of 1) autosomal post QC SNPs (A) and 2) pre QC discordant SNPs (B) - A \u2229 B\n"
##    table += u"g) intersection of 1) autosomal post QC SNPs (A) and 2) post QC discordant SNPs (C) - A \u2229 C\n"
##    table += u"h) intersection of 1) autosomal post QC SNPs (A) and the relative complement of 2) pre QC discordant SNPs with respect to 3) post QC discordant SNP set - A \u2229 (C \ B) or (g-f)\n"
##    table += u"i) autosomal SNPs after QC \u2229 quad BBGZ post QC intersection \u2229 octo BBGZ post QC intersection - aka Manj bullet point 4\n"
##    for bfile in l_bfiles:
##        ## init row
##        table += '%-16s\t' %(bfile)
##        ## a
##        cmd = 'cat ../pops/%s/%s.fam | wc -l' %(bfile,bfile,)
##        i_samples_preQC = int(os.popen(cmd).read())
##        table += '%3i\t' %(i_samples_preQC)
##        ## b
##        table += '%7i\t' %(i_SNPs_autosomal_preQC)
##        ## c
##        cmd = 'cat ../pops/%s/%s.SNPQC.fam | wc -l' %(bfile,bfile,)
##        i_samples_postQC = int(os.popen(cmd).read())
##        table += '%3i\t' %(i_samples_postQC)
##        ## d
##        cmd = 'cat ../pops/%s/%s.SNPQC.bim | wc -l' %(bfile,bfile,)
##        i_SNPs_autosomal_postQC = int(os.popen(cmd).read())
##        table += '%7i\t' %(i_SNPs_autosomal_postQC)
##        ## e
##        table += '%6i\t' %(i_SNPs_autosomal_preQC-i_SNPs_autosomal_postQC)
##        ## f
##        cmd = 'comm -12 %s.SNPQC.bim.sorted discordant_sampleQC.SNPs.sorted | wc -l' %(bfile)
##        table += '%5i\t' %(int(os.popen(cmd).read()))
##        ## g
##        cmd = 'comm -12 %s.SNPQC.bim.sorted discordant_SNPQC.SNPs.sorted | wc -l' %(bfile)
##        table += '%5i\t' %(int(os.popen(cmd).read()))
##        ## h (just do the diff!!!)
##        execmd('comm -12 %s.SNPQC.bim.sorted discordant_sampleQC.SNPs.sorted > tmp4' %(bfile))
##        table += '%5i\t' %(int(os.popen('comm -23 tmp4 discordant_SNPQC.SNPs.sorted | wc -l').read()))
##        ## i
##        cmd = 'comm -12 %s.SNPQC.bim.sorted postQC.incl.SNPs.intersection | wc -l' %(bfile)
##        table += '%7i\t' %(int(os.popen(cmd).read()))
##        ## term row
##        table += '\n'
####        os.remove('%s.SNPQC.bim.sorted' %(bfile))
##    import codecs
####    table = table.replace(u'\u2229',u'2229')
####    table = table.replace(u'\u2229',u'2229')
##    fn_table = 'ms23_allpopsQC98_table.txt'
##    fd = codecs.open(fn_table,'w',encoding='utf-8')
##    fd.write(table)
##    fd.close()
##    execmd('unix2dos %s' %(fn_table))
####    cmd = 'echo "Here is the ascii encoded table.\n\n%s" | mail -s "Table" tc9@sanger.ac.uk' %(str(table))
####    execmd(cmd)

    ##
    ## MDS plots
    ##
    execmd("cat omni2.5-4_20120904_agv_gtu.fam | awk '{print $2}' > quad.samples")
    execmd("cat omni2.5-8_agv_20120910_gtu.fam | awk '{print $2}' > octo.samples")
    for QC in [
##        'sampleQC',
        'SNPQC',
        ]:
        cmd = 'cat Baganda_quad.%s.comm.samples Baganda_octo.%s.comm.samples' %(QC,QC,)
        cmd += ' > Baganda.%s.comm.samples' %(QC)
        execmd(cmd)
        for population in ['Baganda','Ga-Adangbe','Banyarwanda','Zulu',]:
##        for population in l_populations:
##            population = d_populations[population]
            for bool_excl_discordant in [
##                False,
                True,
                ]:
                for bool_excl_MAF_below_5 in [
                    False,
##                    True,
                    ]:
##                    for threshold in range(0,10):
                    for threshold in [4,]:
                        for set_name in [
                            'single','union','intersection',
##                            'intersectBanBag',
##                            'intersect_and_postSNPQCdisc',
                            'union3pops', ## label should have been union3pops
                            'intersect_and_postSNPQCdisc_and_HWE',
                            ]:
                            prefix = '%s_%s_%s_%s_%s_%i' %(population,QC,bool_excl_discordant,bool_excl_MAF_below_5,set_name,threshold,)
                            MDS(population,QC,
                                bool_excl_discordant,bool_excl_MAF_below_5,
                                set_name,threshold,prefix,)
                            prefix_out = '%s_%s_%i_%s_%s_%s' %(
                                bool_excl_discordant,set_name,threshold,
                                population,
                                bool_excl_MAF_below_5,QC,)
                            MDS_plot(population,QC,
                                     bool_excl_discordant,prefix,prefix_out,)

##    ##
##    ## frq.png
##    ##
##    for QC in ['sampleQC','SNPQC',]:
##        Baganda_frq_scatter(QC)
##
##    ##
##    ## png
##    ##
##    frqlmisshwe_discordant()

    return


def allelic_concordance_Baganda(bfile1,bfile2,):

    for bfile in [bfile1,bfile2,]:
        if not os.path.isfile('%s.bed' %(bfile)):
            print 'bed not found:', bfile
            sys.exit(0)

    basename1 = os.path.basename(bfile1)
    basename2 = os.path.basename(bfile2)

    suffix = basename1[basename1.rindex('.')+1:]

    ##
    ## find common post QC samples (chip/lane independent)
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
    ## find common post QC samples (chip/lane specific)
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
    ## find common post QC SNPs
    ##
    fn_common_SNPs = 'Baganda.%s.comm.SNPs' %(suffix,)
##    if not os.path.isfile(fn_common_SNPs):
    for bfile, basename in [[bfile1,basename1,],[bfile2,basename2,],]:
        cmd = "cat %s.bim | awk '{print $2}' | sort > %s.SNPs" %(bfile,basename,)
        execmd(cmd)
    cmd = 'comm -12 %s.SNPs %s.SNPs > %s' %(
        basename1,basename2,fn_common_SNPs,)
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

    ##
    ## parse sample IIDs
    ##
    d_fam = {}
    for bfile in [bfile1,bfile2,]:
        basename = os.path.basename(bfile)
        cmd = "cat %s.tfam | awk '{print substr($2,12,10)}'" %(basename)
        l = os.popen(cmd).read().strip().split('\n')
        d_fam[basename] = l
    l_fam4 = d_fam[basename1]
    l_fam8 = d_fam[basename2]

    ##
    ## initiate count
    ##
    d_sum = {}
    for k in ['either','confirmed','possible',]:
        d_sum[k] = {}
        for count_disc in range(n_samples+1):
            d_sum[k][count_disc] = {'count':0,'sum_conc':{},}
##                                    'count_cum':0,'sum_conc_cum':{},}
            for i_sample in xrange(n_samples):
                d_sum[k][count_disc]['sum_conc'][i_sample] = 0
##                d_sum[k][count_disc]['sum_conc_cum'][i_sample] = 0

    for i_sample in xrange(n_samples):
        fn = 'discordant_%s_%i.SNPs' %(suffix,i_sample,)
        fd = open(fn,'w')
        fd.close()

    l_discordant_SNPs = []

    i_line = 0

    print 'loop over lines/SNPs'

    fd4 = open('%s.tped' %(basename1),'r')
    fd8 = open('%s.tped' %(basename2),'r')
    for line4 in fd4:
        ## read corresponding line/SNP
        line8 = fd8.readline()

        ## parse line
        l4 = line4.split()

        ## skip if not autosomal or X chromosome SNP
        chromosome = int(l4[0])
        if chromosome < 1 or chromosome > 22:
            continue

        ## parse line
        l8 = line8.split()

        ## count line/SNP
        i_line += 1
        if i_line % 100000 == 0: print i_line

        ## initiate counts per line/SNP
        d_disc_counts = {
            'either':0,
            'confirmed':0,
            'possible':0,
            }
        bool_disc_confirmed = False
        for i4 in xrange(n_samples):
            sample = l_fam4[i4]
            ## samples in different order, so...
            ## identify index of octo sample equivalent to quad sample
            i8 = l_fam8.index(sample)
            genotype4 = l4[4+2*i4:4+2*(i4+1)]
            genotype8 = l8[4+2*i8:4+2*(i8+1)]
            ## concordant
            if (
                (genotype4[0] == genotype8[0] and genotype4[1] == genotype8[1])
                or
                (genotype4[0] == genotype8[1] and genotype4[1] == genotype8[0])
                ):
                pass
            ## discordant
            else:
                d_disc_counts['either'] += 1
                if genotype4[0] == '0' or genotype8[0] == '0':
                    d_disc_counts['possible'] += 1
                else:
                    d_disc_counts['confirmed'] += 1

            ## continue loop over samples
            continue

##        l = ['either']
##        if count_diff_confirmed == 0:
##            l += ['possible']
##        if count_diff_possible == 0:
##            l += ['confirmed']
        l = ['either','confirmed','possible',]
        for k in l:
            d_sum[k][d_disc_counts[k]]['count'] += 1
##            for count in range(0,d_disc_counts[k]+1-1):
            d_sum[k][d_disc_counts[k]]['sum_conc'][i_sample] += n_samples-d_disc_counts[k]

        if d_disc_counts['confirmed'] > 0:
            l_discordant_SNPs += [l4[1]]
            ## append to file within loop - slow but less than 50k SNPs
            fn = 'discordant_%s_%i.SNPs' %(suffix,d_disc_counts['confirmed'],)
            fd = open(fn,'a')
            fd.write('%s\n' %(l4[1]))
            fd.close()

        ## continue loop over lines/SNPs
        continue

    fd4.close()
    fd8.close()

    fd = open('discordant_%s.SNPs' %(suffix),'w')
    fd.write('\n'.join(l_discordant_SNPs)+'\n')
    fd.close()

##    s_table = 'confirmed discordances (e.g. TC!=CC or TT!=CC, TC and CT considered to be equivalent)\n'
##    s_table += 'possible discordances (e.g. TC!=00, 00 and 00 considered to be equivalent)\n'
    print 'all', i_line
    print 'all', int(os.popen('cat %s.tped | wc -l' %(basename1)).read())
    print 'all', int(os.popen('cat %s.tped | wc -l' %(basename2)).read())

    ## initiate table rows
    lines_table = ['%2i' %(count) for count in xrange(n_samples+1)]
    ## loop over columns
    for k in ['either','confirmed','possible',]:
        sum_conc_cum = 0
        n_SNPs_cum = 0
        ## loop over rows
        for count in xrange(n_samples+1):
            ## calculate concordance
            sum_conc = 0
            n_SNPs = d_sum[k][count]['count']
            for i in xrange(n_samples):
                sum_conc += d_sum[k][count]['sum_conc'][i]
            ## calculate cumulated concordance
            sum_conc_cum += sum_conc
            avg_sum_conc_cum = float(sum_conc_cum)/n_samples
            n_SNPs_cum += n_SNPs
            if n_SNPs_cum == 0:
                avg_conc_cum = 'N/A'
            else:
                avg_conc_cum = '%.4f' %(float(avg_sum_conc_cum)/n_SNPs_cum)
            ## append to row
            lines_table[count] += '\t%7i\t%7i\t%s' %(
                n_SNPs,n_SNPs_cum,avg_conc_cum,)
            ## continue loop over rows
            continue
        ## continue loop over columns
        continue
    ## terminate table rows
    for count in xrange(n_samples+1):
        lines_table[count] += '\n'
##    ## append table header
##    lines_table = [s_table]+lines_table
    ## write table to file
    fd = open('table%s.txt' %(suffix),'w')
    fd.writelines(lines_table)
    fd.close()
    ####

##    cmd = 'cat'
##    for i_sample in xrange(n_samples):
##        fn = 'discordant_%s_%i.SNPs' %(suffix,i_sample,)
##        cmd += ' %s' %(fn)
##    cmd += ' > discordant_%s.SNPs' %(suffix)

    cmd = 'sort discordant_%s.SNPs > discordant_%s.SNPs.sorted' %(suffix,suffix,)
    execmd(cmd)
    cmd = 'sort %s > %s.sorted' %(fn_common_SNPs,fn_common_SNPs,)
    execmd(cmd)
    execmd('comm -23 %s.sorted discordant_%s.SNPs.sorted > discordant_%s_0.SNPs' %(
        fn_common_SNPs,suffix,suffix,))
    cmd = "awk '{print $1,substr(FILENAME,18,length(FILENAME)-22)}'"
    cmd += ' discordant_%s_*.SNPs > discordant_%s.dg11' %(suffix,suffix,)
    execmd(cmd)

    return n_samples, fn_common_SNPs


def execmd(cmd,bool_exe=True,):

    l_indexes = []
    l_cmd = cmd.split()

    if '%' in cmd:
        print cmd
        stop

    if cmd.split()[0] in ['cat','mv','cp','rm',]:
        l_indexes = [1]
    elif cmd.split()[0] == 'paste':
        l_indexes = [1,2,]
    elif cmd.split()[0] == 'comm':
        l_indexes = [2,3,]
    elif cmd.split()[0] in ['grep','fgrep',] and '-f' in cmd.split():
        l_indexes = [l_cmd.index('-f')+1]
    for index in l_indexes:
        if not os.path.isfile(cmd.split()[index]):
            print cmd
            print inspect.stack()[1][3]
            print 'does not exist:', cmd.split()[index]
            sys.exit(0)

    if bool_exe == True:
        print inspect.stack()[1][3]
        print cmd
        os.system(cmd)

    return


def HWE_merged():

    QC = 'SNPQC'
    for population in ['Zulu','Banyarwanda','Baganda','Ga-Adangbe',]:

        bfile_out = '%s_%s' %(population,QC,)

        l_cmds = []
    
        ##
        ## 1) find common SNPs
        ##
        fn_common_SNPs = '%s.comm.SNPs' %(bfile_out,)
        l_basenames = []
        for chip in ['quad','octo',]:
            bfile_in = '../pops/%s_%s/%s_%s.%s' %(population,chip,population,chip,QC,)
            basename = os.path.basename(bfile_in)
            l_basenames += [basename]
            cmd = "cat %s.bim | awk '{print $2}' | sort > %s.SNPs.sorted" %(
                bfile_in,basename,)
            l_cmds += [cmd]
        cmd = 'comm -12 %s.SNPs.sorted %s.SNPs.sorted > %s' %(
            l_basenames[0],l_basenames[1],fn_common_SNPs,)
        l_cmds += [cmd]

        ##
        ## 2) merge and extract common SNPs
        ##
        cmd = 'plink \\\n'
        cmd += '--bfile ../pops/%s_quad/%s_quad.%s \\\n' %(population,population,QC)
        cmd += '--bmerge \\\n'
        cmd += '../pops/%s_octo/%s_octo.%s.bed \\\n' %(population,population,QC)
        cmd += '../pops/%s_octo/%s_octo.%s.bim \\\n' %(population,population,QC)
        cmd += '../pops/%s_octo/%s_octo.%s.fam \\\n' %(population,population,QC)
        if population == 'Baganda':
            cmd += '--remove Baganda.%s.comm.samples \\\n' %(QC)
        cmd += '--extract %s \\\n' %(fn_common_SNPs)
        cmd += '--out %s \\\n' %(bfile_out)
        cmd += '--hardy \\\n'
        l_cmds += [cmd]

        ##
        ## 3) find HWE outliers
        ##
        cmd = 'cat %s.hwe' %(bfile_out)
        cmd += " | awk 'NR>1&&NR%3==2&&$9<0.0001{print $2}'"
        cmd += ' | sort'
        cmd += ' > %s.hwe.SNPs' %(bfile_out)
        l_cmds += [cmd]

        fn_sh = '%s_hwe.sh' %(bfile_out,)
        write_shell(fn_sh,l_cmds,)

        cmd = LSF_append('./%s' %(fn_sh),queue='yesterday',)
        execmd(cmd)

    return


def super_MDS(d_pops2coords,l_bfiles,l_populations,):

    d_l_popset2bfiles = init_MDS(l_populations)

    fn_ld_regions = 'ldregions.SNPs'

    l_population_sets = list(d_l_popset2bfiles.keys())
    l_population_sets = [
        'Africa',
##        'Africa1000G',
##        'BBGZ',
##        '1000G',
        ]

    ##
    ## 1) perform MDS
    ##
    for population_set in l_population_sets:

        l_populations = d_l_popset2bfiles[population_set]

        for bool_exclude in [
            False,
            True,
            ]:

            bfile = '%s_%s' %(population_set,str(bool_exclude))

            super_MDS_PLINK(bfile,l_populations,bool_exclude,)

            ## continue loop over False/True exclusion of discordant SNPs
            continue

        ## continue loop over populations
        continue

##    ##
##    ## 3) find component with greatest "quad/octo correlation"
##    ##
##    for bfile in l_bfiles:
##        cmd = 'cat %s.fam | sort -k1,1 > %s.fam.sorted' %(bfile,bfile,)
##        execmd(cmd)
##    d_max_correlation = {}
##    for population_set in l_population_sets:
##
##        for bool_exclude in [
##            False,
##            True,
##            ]:
##
##            bfile = '%s_%s' %(population_set,str(bool_exclude))
##
##            if not os.path.isfile('%s.mds' %(bfile)):
##                continue
##
##            cmd = 'cat %s.mds | sort -k1,1 > %s.mds.sorted' %(bfile,bfile,)
##            execmd(cmd)
##            l_keys = ['quad','octo',]
##            d_arrays = {'quad':None,'octo':None,}
##            for i in xrange(2):
##                bfile_main = l_bfiles[i]
##                fn_array = '%s.%s.mds' %(bfile,l_keys[i],)
##                cmd = 'join -1 1 -2 1 %s.fam.sorted %s.mds.sorted' %(
##                    bfile_main,bfile,)
##                cmd += ' > %s' %(fn_array)
##                execmd(cmd)
##                cmd = 'cat %s | wc -l' %(fn_array)
##                i_lines = int(os.popen(cmd).read())
##                if i_lines == 0:
##                    continue
##                ## read array
##                cmd = 'cat %s' %(fn_array)
##                cmd += " | sed 's/ /\\t/g' | cut -f%i-" %(1+(6-1)+(3-1)+1)
##                s = os.popen(cmd).read()
######                print s[:12]
######                print s[-12:]
##    ##            array = numpy.fromstring(s,sep='\n').reshape((1+(6-1)+(3-1)+1+4,i_lines))
######                print len(numpy.fromstring(s,sep='\n')), 500, i_lines
##                array = numpy.fromstring(s,sep='\n').reshape((i_lines,500))
##                d_arrays[l_keys[i]] = array
##            l = []
##            for component in xrange(500):
##                l_quad = d_arrays['quad'][:,component]
##                l_octo = d_arrays['octo'][:,component]
##                l1 = list(l_quad)+list(l_octo) ## numpy.column_stack
##                l2 = len(l_quad)*[0]+len(l_octo)*[1] ## numpy.tile
##    ##            instance_tests = statistics.tests()
##    ##            a,b,r,p_correlation = instance_tests.do_regression(l1,l2)
##    ##            print a,b,r,p_correlation
##    ##            r = instance_tests.correlation(l1,l2)
##    ##            print r
##                ndarray = numpy.corrcoef(l1,l2)
##                r = ndarray[0][1]
##                l += [[abs(r),component]]
##                l.sort()
##            d_max_correlation[bfile] = [
##                l[-1][1],l[-2][1],
##                l[-1][0],l[-2][0],]
##            fd = open('tmp.txt','a')
##            fd.write('%s %s\n' %(bfile,str(d_max_correlation[bfile])))
##            fd.close()

    ##
    ## 2) plot MDS
    ##
    l_colors, d_bfiles2chips, d_bfiles2pts = init_MDS_plot()

    for population_set in l_population_sets:

        for bool_exclude in [
            False,
            True,
            ]:

            bfile = '%s_%s' %(population_set,str(bool_exclude))

            if not os.path.isfile('%s.mds' %(bfile)):
                continue
             
            ##
            ## 6) plot
            ##
            MDS_super_plot(
                population_set,bool_exclude,
                d_l_popset2bfiles,
                l_colors, d_bfiles2chips, d_bfiles2pts,
                d_pops2coords,
                )

    stop_end
    return


def genome(bfile):

    if not os.path.isdir('genome'):
        os.mkdir('genome')

    prefix = bfile

    l_cmds = []

    ## parse FID and IID columns
    cmd = 'cat %s.fam' %(bfile,)
    cmd += ' | '
    cmd += " awk '{print $1,$2}' "
    ## and pipe them to individual fam files
    cmd += ' | split -d -a 2 -l %i - genome/%s.fam.' %(
        400,bfile,
        )
    l_cmds += [cmd]

    ##
    ## count number of fam files generated
    ##
    cmd = 'n=$(cat %s.fam | wc -l)' %(bfile,)
    l_cmds += [cmd]
    cmd = 'cnt=$(((${n}+%i-1)/%i))' %(
        400, 400,
        )
    l_cmds += [cmd]

    ##
    ## loop 1 - execute plink
    ##
    cmd = '\n##\n## loop 1 - run plink\n##\n'
    cmd += 'for i in $(seq -f "%02g" 0 $(($cnt-1)))\ndo\n'
    cmd += 'for j in $(seq -f "%02g" 0 $(($cnt-1)))\ndo\n'
    cmd += 'if [ $j -gt $i ]\nthen continue\nfi\n\n'
    cmd += '## continue if output exists\n'
    cmd += 'if [ -s genome/%s.$i.$j.genome ]\nthen continue\nfi\n\n' %(
        bfile,
        )
    l_cmds += [cmd]

##    if not os.path.isfile('%s.genome' %(prefix)):    
    cmd = 'plink \\\n'
    cmd += '--bfile %s \\\n' %(prefix)
    cmd += '--out genome/%s.$i.$j \\\n' %(prefix)
    cmd += '--genome \\\n'
    ## SNP exclusion
    cmd += '--extract %s.prune.in \\\n' %(prefix)
    cmd += '--genome-lists genome/%s.fam.$i genome/%s.fam.$j \\\n' %(bfile,bfile,)

    cmd = LSF_append(cmd)
    l_cmds += [cmd]

    l_cmds += ['\ndone\ndone\n']

    ##
    ## loop2 - count lines of plink output
    ##
    cmd = '\n##\n## loop 2 - count lines\n##\n'
    cmd += 'nlines=0\n'
    cmd += 'break=0\n\n'
    cmd += 'for i in $(seq -f "%02g" 0 $(($cnt-1)))\ndo\n'
    cmd += 'if [ $break -eq 1 ]\nthen\nbreak\nfi\n'
    cmd += 'for j in $(seq -f "%02g" 0 $(($cnt-1)))\ndo\n'
    cmd += 'if [ $j -gt $i ]\nthen continue\nfi\n\n'
    cmd += 'if [ ! -s genome/%s.$i.$j.genome ]\nthen\n' %(bfile)
    cmd += 'break=1\nbreak\nfi\n\n'
    cmd += 'let nlines=$nlines+$('
    cmd += 'cat genome/%s.$i.$j.genome | wc -l' %(bfile)
    cmd += ')-1\n'
    cmd += '\ndone\ndone\n'
    l_cmds += [cmd]

    cmd += 'echo actual $nlines expected $(($n*($n-1)/2))\n'
    cmd += 'if [ $nlines -ne $(($n*($n-1)/2)) ];then\nexit\nfi\n'
    l_cmds += [cmd]

    ##
    ## loop 3 - combine plink output lines
    ##
    s = '                   FID1                   IID1                   FID2                   IID2 RT    EZ      Z0      Z1      Z2  PI_HAT PHE       DST     PPC   RATIO'
    cmd = '\n##\n## loop 3 - concatenate genome files\n##\n'
    cmd += 'echo "%s" > %s.genome\n' %(s,bfile,)
    ## loop .genome files
    cmd += 'for i in $(seq -f "%02g" 0 $(($cnt-1)))\ndo\n'
    cmd += 'for j in $(seq -f "%02g" 0 $(($cnt-1)))\ndo\n'
    cmd += 'if [ $j -gt $i ]\nthen continue\nfi\n\n'
    ## append to .genome file
    cmd += "cat genome/%s.$i.$j.genome | awk 'NR>1' >> %s.genome\n" %(
        bfile,bfile,
        )
    cmd += '\ndone\ndone\n\n'
    l_cmds += [cmd]

    return l_cmds


def init_MDS_plot():

    l_colors = [
        [255,0,0,],
        [255,85,0,],
        [255,170,0,],
        [255,255,0,],
        [170,255,0,],
        [85,255,0,],
        [0,255,0,],
        [0,255,85,],
        [0,255,170,],
        [0,255,255,],
        [0,170,255,],
        [0,85,255,],
        [0,0,255,],
        [85,0,255,],
        [170,0,255,],
        [255,0,255,],
        [255,0,170,],
        [255,0,85,],
        [0,0,0,],
        [85,85,85,],
        [170,170,170,],
        ]

    d_bfiles2chips = {
        ## Uganda
        'Baganda_quad':'quad', 'Baganda_octo':'octo',
        'Barundi':'octo',
        'Banyarwanda_quad':'quad', 'Banyarwanda_octo':'octo',
        ## Kenya
        'Kalenjin':'quad', 'Kikuyu':'quad',
        ## South Africa
        'Zulu_quad':'quad', 'Zulu_octo':'octo', 'Sotho':'octo',
        ## Nigeria
        'Igbo':'quad',
        ## Ghana
        'Ga-Adangbe_quad':'quad', 'Ga-Adangbe_octo':'octo',
        ##
        'Mandinka':'octo', 'Wolof':'octo', 'Fula':'octo', 'Jola':'octo',
        ## Ethiopia
        'Ethiopia':'octo', ## SOMALI,AMHARA,OROMO
        ## 1000g
        'YRI':'quad', 'CHB':'quad', 'ASW':'quad', 'TSI':'quad', 'CDX':'quad', 'CLM':'quad', 'CEU':'quad', 'KHV':'quad', 'PEL':'quad', 'LWK':'quad', 'MXL':'quad', 'CHS':'quad', 'GBR':'quad', 'ACB':'quad', 'IBS':'quad', 'FIN':'quad', 'JPT':'quad', 'PUR':'quad', 'GIH':'quad',
        }

    d_bfiles2pts = {
        ## Uganda
        'Baganda_quad':'East', 'Baganda_octo':'East',
        'Barundi':'East',
        'Banyarwanda_quad':'East', 'Banyarwanda_octo':'East',
        ## Kenya
        'Kalenjin':'East', 'Kikuyu':'East',
        ## South Africa
        'Zulu_quad':'South', 'Zulu_octo':'South', 'Sotho':'South',
        ## Nigeria
        'Igbo':'Central',
        ## Ghana
        'Ga-Adangbe_quad':'Central', 'Ga-Adangbe_octo':'Central',
        ##
        'Mandinka':'West', 'Wolof':'West', 'Fula':'West', 'Jola':'West',
        ## Ethiopia
        'Ethiopia':'East', ## SOMALI,AMHARA,OROMO
        ## 1000g
        'CHB':'quad', 'ASW':'quad', 'TSI':'quad', 'CDX':'quad', 'CLM':'quad', 'CEU':'quad', 'KHV':'quad', 'PEL':'quad', 'MXL':'quad', 'CHS':'quad', 'GBR':'quad', 'ACB':'quad', 'IBS':'quad', 'FIN':'quad', 'JPT':'quad', 'PUR':'quad', 'GIH':'quad',
        'YRI':'Central', 'LWK':'East',
        }

    return l_colors, d_bfiles2chips, d_bfiles2pts


def MDS_super_plot(
    population_set,bool_exclude,
    d_l_popset2bfiles,
    l_colors, d_bfiles2chips, d_bfiles2pts,
    d_pops2coords,
    ):

    bfile1 = '%s_%s' %(population_set,str(bool_exclude))
            
    d_pops2style = {
        }

    ##
    ## sort mds file before join with fam file
    ##
    cmd = "cat %s.mds | sort -k1,1 -k2,2 > %s.mds.sorted" %(
        bfile1,bfile1,)
    execmd(cmd)

    ##
    ## sort populations
    ##
##            l_populations = d_l_popset2bfiles[population_set]
    l_populations = []
##            for chip in ['quad','octo',]:
##                for population in d_l_popset2bfiles[population_set]:
##                    if d_bfiles2chips[population] == chip:
##                        l_populations += [population]
    l = []
    if 'Sotho' in d_l_popset2bfiles[population_set]:
        l_bfiles = ['Sotho','Zulu_octo','Zulu_quad',]
    else:
        l_bfiles = []
    for bfile2 in d_l_popset2bfiles[population_set]:
        if bfile2 in ['Sotho','Zulu_quad','Zulu_octo',]:
            continue
        pop = bfile2.replace('_quad','').replace('_octo','')
        if not pop in d_pops2coords.keys():
            coord = 'N/A'
        else:
            coord = d_pops2coords[pop][1]
        l += [[coord,bfile2,]]
    l.sort()
    l_bfiles += [l2[1] for l2 in l]

    if not os.path.isfile('map_%s.png' %(population_set)):
        s_map = Africa_map_init()

    for c1 in xrange(4):
        for c2 in xrange(4):
            if c2 <= c1:
                continue
            if os.path.isfile('mds.3D.%s.%i.%i.%i.png' %(bfile1,c1,c2,c2+1,)):
                continue
##    for c1 in [0]+[d_max_correlation[bfile1][0]]:
##        for c2 in [1]+[d_max_correlation[bfile1][1]]:
##            if not (
####                (c1 == 0 and c2 == 1)
######                or
######                (c1 == 1 and c2 == 2)
####                or
##                (
##                    c1 == d_max_correlation[bfile1][0]
##                    and
##                    c2 == d_max_correlation[bfile1][1]
##                    )
##                ):
##                continue
            line_plot = 'plot '
            line_splot = 'splot '
            for i in xrange(len(l_bfiles)):
                bfile2 = l_bfiles[i]
##                        if 'quad' not in population and 'octo' not in population: continue
                if d_bfiles2chips[bfile2] == 'quad':
                    ps = 2
                    l_pt = [5,7,]
                else:
                    ps = 3
                    l_pt = [9,11,]
                if d_bfiles2pts[bfile2] == 'East':
                    l_pt= [5]
                elif d_bfiles2pts[bfile2] == 'West':
                    l_pt= [7]
                elif d_bfiles2pts[bfile2] == 'South':
                    l_pt= [9]
                elif d_bfiles2pts[bfile2] == 'Central':
                    l_pt= [11]
                else:
                    l_pt = [7]
                pt = l_pt[i%len(l_pt)]
                if 'quad' in bfile2 or 'octo' in bfile2:
                    ps = 3
                color = l_colors[i%len(l_colors)]
##                if bool_exclude == False:
##                    if not bfile2.replace('_quad','').replace('_octo','') in d_pops2coords.keys():
##                        color = [0,0,0,]
##                    elif d_bfiles2chips[bfile2] == 'quad':
##                        color = [255,0,0,]
##                    else:
##                        color = [0,0,255,]
                color = "".join(map(chr, color)).encode('hex')
                d_pops2style[bfile2.replace('_quad','').replace('_octo','')] = {
                    'color':color,'pt':pt,'ps':ps,
                    }
                
                cmd = 'cat ../pops/%s/%s.postQC.autosomes.fam | sort -k1,1 -k2,2 > %s.fam.sorted' %(
                    bfile2,bfile2,bfile2,)
                execmd(cmd)
                cmd = 'join -1 2 -2 2'
                cmd += ' %s.fam.sorted %s.mds.sorted ' %(bfile2,bfile1,)
                cmd += ' > %s.%s.mds' %(bfile1,bfile2,)
                execmd(cmd)
##                os.remove('%s.fam.sorted' %(bfile2))
                line_plot += '"%s.%s.mds" u %i:%i lc rgb "#%s" pt %i ps %i t "%s (%s)",' %(
                    bfile1,bfile2,
                    1+(6-1)+(3-1)+1+c1,
                    1+(6-1)+(3-1)+1+c2,
                    color,
                    pt,
                    ps,
                    bfile2.replace('_quad','').replace('_octo',''),
                    d_bfiles2chips[bfile2],
                    )
##                if bfile2 in ['Baganda_quad','Banyarwanda_quad','Kikuyu','Kalenjin','LWK','Zulu_quad',]: color = 'ff0000' ## quad/red
##                elif bfile2 in ['Baganda_octo','Banyarwanda_octo','Barundi','Zulu_octo','Sotho',]: color = 'ff8080' ## octo/pink
##                elif bfile2 in ['Ethiopia',]: color = '00ff00'
##                elif bfile2 in ['Igbo','YRI','Ga-Adangbe_quad',]: color = '0000ff' ## blue/quad
##                elif bfile2 in ['Ga-Adangbe_octo',]: color = '8080ff' ## bright blue/octo
##                elif bfile2 in ['Mandinka','Wolof','Jola','Fula',]: color = 'ffff00' ## orange/octo
                line_splot += '"%s.%s.mds" u %i:%i:%i lc rgb "#%s" pt %i ps %i t "%s (%s)",' %(
                    bfile1,bfile2,
                    1+(6-1)+(3-1)+1+c1,
                    1+(6-1)+(3-1)+1+c2,
                    1+(6-1)+(3-1)+1+c2+1,
                    color,
                    pt,
                    ps,
                    bfile2.replace('_quad','').replace('_octo',''),
                    d_bfiles2chips[bfile2],
                    )
            line_plot = line_plot[:-1]+'\n'
##                    for z in range(5):
##                        line_splot += '0.0%i lc rgb "#bbbbbb" t "",' %(z)
##                        line_splot += '-0.0%i lc rgb "#bbbbbb" t "",' %(z)
            for z in [0.01,0.00,-0.01,-0.02,-0.03]:
                line_splot += '%f lc rgb "#bbbbbb" t "",' %(z)
            line_splot = line_splot[:-1]+'\n'

            n_samples = int(os.popen('cat %s.fam | wc -l' %(bfile1)).read())-1
            n_SNPs1 = int(os.popen('cat %s.bim | wc -l' %(bfile1)).read())-1
            n_SNPs2 = int(os.popen('cat %s.prune.in | wc -l' %(bfile1)).read())-1

            title = '%s (n_{samples}=%i, n_{SNPs_{before pruning}}=%i, n_{SNPs_{after pruning}}=%i)' %(
                population_set,
                n_samples,n_SNPs1,n_SNPs2,
                )

            gnuplot.scatter_plot_2d(
                '%s.mds' %(bfile1),
                line_plot = line_plot,
        ##        column1 = 2, column2 = 3,
                xlabel = 'C%i' %(c1+1),
                ylabel = 'C%i' %(c2+1),
                title=title,
                prefix_out='mds.2D.%s.%i.%i' %(bfile1,c1,c2),
        ##        lines_extra=lines_extra,
        ##        bool_remove=False,
        ##        bool_execute=False,
                lines_extra=['set key out\n'],
                )

            if (
                (c1 == 0 and c2 == 1)
                or
                (c1 == 1 and c2 == 2)
                ):
                gnuplot.contour_plot(
                    'mds.3D.%s.%i.%i.%i' %(bfile1,c1,c2,c2+1,),
                    line_splot = line_splot,
                    xlabel = 'C%i' %(c1+1),
                    ylabel = 'C%i' %(c2+1),
                    zlabel = 'C%i' %(c2+1+1),
    ##                        bool_remove = False,
                    title=title,
                    )

            for i in xrange(len(l_bfiles)):
                bfile2 = l_bfiles[i]
                os.remove('%s.%s.mds' %(
                    bfile1,bfile2,))
                os.remove('%s.fam.sorted' %(bfile2))

    if not os.path.isfile('map_%s.png' %(population_set)):
        Africa_map_term(
            population_set,s_map,d_pops2style,d_pops2coords,)

    return


def Africa_map_term(population_set,s,d_pops2style,d_pops2coords,):

    d_pops2style['Baganda']['color'] = '01559D'
    d_pops2style['Zulu']['color'] = '01559D'
    d_pops2style['Fula']['color'] = '01559D'
    d_pops2style['Ethiopia']['color'] = '01559D'
    d_pops2style['Baganda']['pt'] = 7
    d_pops2style['Zulu']['pt'] = 7
    d_pops2style['Fula']['pt'] = 7
    d_pops2style['Ethiopia']['pt'] = 7

    for pop in d_pops2style.keys():
        if not pop in d_pops2coords.keys():
            continue
        if not pop in ['Baganda','Zulu','Fula','Ethiopia',]: continue
##        s += '"-" lc rgb "#%s" pt %s ps %s\n%f %f,' %(
        s += '''"<echo '%f,%f'" lc rgb "#%s" pt %i ps %i,''' %(
            d_pops2coords[pop][1],d_pops2coords[pop][0],
            d_pops2style[pop]['color'],
            d_pops2style[pop]['pt'],
##            5,#d_pops2style[pop]['ps'],
            5,#d_pops2style[pop]['ps'],
            )
        x = d_pops2coords[pop][1]
        y = d_pops2coords[pop][0]
        if True:
            y += 1.25
        elif pop in ['Kalenjin','Ga-Adangbe',]:
##            y -= 1
            y -= 2.25
        elif pop in ['Jola',]:
##            y -= 1
            y -= 1.75
            x -= 1.50
        elif pop in ['LWK']:
            x += 0.25
            y += 1.25
        elif pop in ['YRI']:
            x -= 0.5
            y += 0.75
        elif pop in ['Igbo']:
            x += 0.50
            y += 1.20
        else:
            y += 1
        if True:
            pass
        elif pop in ['Mandinka',]:
##            x += 2.5
            x += 4.00
            y -= 0.25
        elif pop in ['Kikuyu',]:
            x += 2.75
            y -= 0.00
        elif pop in ['Wolof',]:
            y += 0.00
            x -= 1.50
        elif pop in ['Baganda',]:
            x -= 4.50
            y -= 0.25
        elif pop in ['Barundi',]:
            y -= 2.75
        elif pop in ['Banyarwanda',]:
##            x -= 1.75
            x -= 6.25
            y -= 0.75
        elif pop in ['Sotho',]:
            x -= 1.00
            y += 0.25
        elif pop in ['Zulu',]:
            x += 0.50
            y += 0.25
        elif pop in ['Ethiopia',]:
            y += 0.25
        s += '''"<echo '%f,%f,%s'" lc rgb "#000000" w labels font "Helvetica,36",''' %(
            x,y,
##            pop,
##            pop.replace('LWK','Luhya').replace('YRI','Yoruba'),
            pop.replace(
                'Baganda','Uganda\\\\n\\\\nOmni 4778\\\\nHiSeq 100\\\\nIntersection 94').replace(
                'Fula','Fula\\\\n\\\\nOmni 74\\\\nHiSeq 67\\\\nIntersection 0').replace(
                'Zulu','Zulu\\\\n\\\\nOmni 95\\\\nHiSeq 100\\\\nIntersection 95').replace(
                'Ethiopia','Ethiopia\\\\n\\\\nOmni 108\\\\nHiSeq 120\\\\nIntersection 63'),
            )
    s = s[:-1]+'\n'

    fd = open('map.plt','w')
    fd.writelines([
        ## set terminal and output
##        'set terminal png size 1280,1280\n',
        'set terminal postscript\n',
##        'set output "map_%s.png"\n' %(population_set),
        'set output "map_%s.ps"\n' %(population_set),
        ## hide legend keys
        'unset key\n',
        'set size square\n',
        'set size 3,3\n',
        'set encoding iso_8859_1\n', ## postscript encoding *necessary* for special characters (e.g. Angstrom)
        ## hide border
        'set noxtics\n',
        'set noytics\n',
        'set noborder\n',
    ####    ## background color
    ##    'set object 20 rect from graph 0, 0, 0 to graph 1, 1, 0 behind lw 1.0 fc  rgb "blue"  fillstyle  solid 0.15 border -1\n',
        'set datafile separator ","\n',
        s,
        ])
    fd.close()
    execmd('gnuplot map.plt')
    if os.path.isfile('map_%s.png' %(population_set)):
        os.remove('map_%s.png' %(population_set))
    execmd('convert map_%s.ps map_%s.png' %(population_set,population_set))
    execmd('rm map_%s.ps' %(population_set))

    return    


def plots(l_populations):

    print 'plots'

    ##
    ## heterozygosities
    ##
    plot_heterozygosities(l_populations)

    ##
    ##
    ##
    execmd('cp ../pops/%s/%s.autosomes.SNPs agv.autosomes.SNPs' %(
        l_populations[0],l_populations[0],))

    for suffix in [
        ## samples (concatenate)
        'imiss','het','sexcheck',
        ## 'prehardy.genome',
        ## SNPs (paste/join)
##        'frq','hwe','SNPQC.lmiss',
        'fam',
##        'mds',
        'imiss.samples','het.samples','sexcheck.samples','sampleQC.samples',
        ]:
        
        bool_continue = False
        for pop in l_populations:
            fn = '../pops/%s/%s.%s' %(pop,pop,suffix,)
            if not os.path.isfile(fn):
                print 'not found', fn
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

        ##
        ## add header if one exists
        ##
        if not suffix in ['fam',] and not '.samples' in suffix:
##            cmd = 'head -1 %s.%s > agv.%s' %(l_populations[0],suffix,suffix,)
            if suffix in ['imiss','het','sexcheck',]:
                cmd = 'head -1 ../pops/%s/%s.%s > agv.%s' %(
                    l_populations[0],l_populations[0],suffix,suffix,)
            else:
                stop1
                cmd = 'head -1 ../pops/%s/%s.%s > agv.%s' %(
                    l_populations[0],l_populations[0],suffix,suffix,)
            execmd(cmd)
        for pop in l_populations:
            fn = '../pops/%s/%s.%s' %(pop,pop,suffix,)
            ## no header
            if suffix in ['fam',] or '.samples' in suffix:
                cmd = "cat %s >> agv.%s" %(fn,suffix,)
            ## header
            else:
                cmd = "sed '1d' %s >> agv.%s" %(fn,suffix,)
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

    ##
    ## Venn diagram of excluded samples
    ##
    l_het = os.popen('cat agv.het.samples').read().split()
    l_imiss = os.popen('cat agv.imiss.samples').read().split()
    l_sexcheck = os.popen('cat agv.sexcheck.samples').read().split()
    gnuplot.venn3(
        l1=l_imiss,l2=l_het,l3=l_sexcheck,
        suffix='agv',
        text1='call rate',text2='heterozygosity',text3='sex',
##        bool_labels=False,
        )

    return


def plot_heterozygosities(l_populations):

    
##    saturation = 240
##    luminance = 120
##    l = [
##        '%s %i %i %i %i\n' %(
##            l_populations[i],i,
##            int( (240.*i)/(len(l_populations)-1) ), saturation, luminance,
##            )
##        for i in xrange(len(l_populations))
##        ]

    ##
    ## assign number (z) to each population (key) and sort populations alphabetically
    ##
    fd = open('populations.sorted','w')
    fd.writelines(
        ['%s %i\n' %(l_populations[i],i) for i in xrange(len(l_populations))])
    fd.close()
    cmd = 'sort -k1,1 populations.sorted -o populations.sorted'
    execmd(cmd)

    ##
    ## assign population (key) and heterozygosity (y) to each sample
    ##
    cmd = "awk '"
    cmd += ' FNR>1{'
    cmd += ' split(FILENAME,a,"/"); pop=a[3]; het=$5; IID=$2; print pop,het,IID'
    cmd += " }'"
    cmd += ' ../pops/*/*.het'
    ## sort by sample ID
    cmd += ' | sort -k3,3'
    ## sort 
    cmd += ' > agv.het.sorted'
    execmd(cmd)

    ##
    ## exclude low sample call rates and heterozygosity outliers
    ## before assigning numbers
    ## outer join
    ##
    cmd = 'cat agv.imiss.samples agv.het.samples | sort -u > agv.imiss.het.samples'
    execmd(cmd)
    cmd = 'join -1 3 -2 2 -v1 agv.het.sorted agv.imiss.het.samples'
    ## sort by population for inner join
    cmd += ' | sort -k2,2'
    cmd += ' > agv.het.joined1'
    execmd(cmd)

    ##
    ## assign population number (z) and heterozygosity (y) to each sample
    ## and add number (x) to each sample
    ## inner join
    ##
    cmd = 'join -1 2 -2 1 agv.het.joined1 populations.sorted'
    ## sort by pop ID and sample ID for NR=x
    cmd += ' | sort -k4n,4 -k2,2'
    cmd += " | awk '"
    cmd += ' {print NR,$0}'
    cmd += " '"
    cmd += " | sed 's/ /\\t/g' "
    cmd += ' > agv.het.joined2'
    execmd(cmd)

    cmd = 'cat agv.het.joined2'
    cmd += " | awk '"
    cmd += ' BEGIN{i1=1}'
    cmd += ' {if(NR==1) pop=$2}'
    cmd += ' {if($2!=pop) {i2=NR-1; print "\\""pop"\\"",0.5* ( i1+i2 ) ; i1=NR ; pop=$2}}'
    cmd += ' END{i2=NR; print "\\""pop"\\"",0.5* ( i1+i2 )}'
    cmd += "'"
    s_xtics = ', '.join(os.popen(cmd).read().split('\n'))[:-2].replace('_',' ')

    lines = []
    ## http://www.gnuplotting.org/using-a-palette-as-line-color/
    lines += ['h1 = 0/360.\n']
    lines += ['h2 = 320/360.\n']
    lines += ['set palette model HSV functions (1-gray)*(h2-h1)+h1,1,1.0\n']
    lines += ['set boxwidth 1.0 relative\n']
    lines += ['set style fill solid 1.0\n']
##    lines += ['unset xtics\n']
    ## hide cbrange
    lines += ['unset colorbox\n']
    lines += ['set xtics (%s)\n' %(s_xtics)] ## or xtic(col) http://stackoverflow.com/questions/4805930/making-x-axis-tics-from-column-in-data-file-in-gnuplot
    lines += ['set xtics rotate by -45 left autojustify\n']
    lines += ['set xtics font "Helvetica,30"\n']
    lines += ['set bmargin 12\n']
    lines1 = list(lines)
    lines2 = list(lines)
    ## http://gnuplot.sourceforge.net/demo_cvs/varcolor.html
    ## "2D plots cannot color by Z value" unless lc pal z instead of lc pal
    lines1 += ['plot [:][0:0.25]"agv.het.joined2" u 1:4:5 w boxes lc pal z notitle\n']
    lines2 += ['''plot [:][0:0.25]"< awk '{if($2==\"LWK\"||$2==\"YRI\"||length($2)>3) print}' agv.het.joined2" u 1:4:5 w boxes lc pal z notitle\n''']
    for lines,prefix in [
        [lines1,'all',],
        [lines2,'africa',],
        ]:
        gnuplot.plot_and_convert(
            lines,'het_impulses_%s' %(prefix),bool_remove=False,
            ylabel='heterozygosity',
            )

    ##
    ## clean up
    ##
    for fn in [
        'populations.sorted',
        'agv.het.sorted',
        'agv.het.joined1',
        'agv.het.joined2',
        ]:
        os.remove(fn)

    return


def tables(l_populations,):

    print 'tables'

    s_sample = ''
    s_SNP = ''
    d_tables = {
        'imiss':'','het':'','sexcheck':'',
        'lmiss.SNPQC':'',
        'genome':'',
        'lmiss.X.males':'',
        'lmiss.X.females':'',
        'hwe':'',
        'hwe.X.females':'',
        }
    pop0 = l_populations[0]
    cmd = 'cat ../pops/%s/%s.hwe.SNPs ../pops/%s/%s.lmiss.SNPs' %(
        pop0,pop0,pop0,pop0,)
    cmd += ' | sort > SNPQC.SNPs.comm'
    execmd(cmd)
    cmd = 'cp SNPQC.SNPs.comm SNPQC.SNPs.union'
    execmd(cmd)
    for pop in l_populations:
        print 'table', pop
        s_sample += '%20s\t' %(pop)
        s_SNP += '%20s\t' %(pop)

        ##
        ## samples
        ##
        l = []
        l_suffixes_sample = ['imiss','het','sexcheck','sampleQC','genome.0.05']
        for i_suffix in xrange(len(l_suffixes_sample)):
            suffix = l_suffixes_sample[i_suffix]
            if i_suffix == 0 or suffix == 'sampleQC':
                cmd = 'echo "" > %s.%s.samples.concatenated.sorted' %(pop,suffix,)
            else:
                cmd = 'cat'
                for i_suffix2 in xrange(i_suffix):
                    suffix2 = l_suffixes_sample[i_suffix2]
                    cmd += ' ../pops/%s/%s.%s.samples' %(
                        pop,pop,suffix2,)
                cmd += ' | sort > %s.%s.samples.concatenated.sorted' %(pop,suffix,)
            execmd(cmd)
            cmd = 'cat ../pops/%s/%s.%s.samples' %(
                pop,pop,suffix,)
            cmd += ' | sort > %s.%s.samples.sorted' %(pop,suffix,)
            execmd(cmd)
##            cmd = 'cat pops/%s/%s.%s.samples | wc -l' %(
##                population,population,suffix,)
            cmd = 'comm -23 %s.%s.samples.sorted %s.%s.samples.concatenated.sorted | wc -l' %(
                pop,suffix,pop,suffix,)
            l += ['%s' %(os.popen(cmd).read().strip()),]
            os.remove('%s.%s.samples.sorted' %(pop,suffix,))
            os.remove('%s.%s.samples.concatenated.sorted' %(pop,suffix,))
        s_sample += '\t'.join(l)
        s_sample += '\n'

        ##
        ## SNPs
        ##
        l = []
        cmd = 'cat ../pops/%s/%s.bim | wc -l' %(pop,pop,)
        execmd(cmd,bool_exe=False)
        l += ['%s' %(os.popen(cmd).read().strip()),]
        for suffix in [
##            'X.SNPs','autosomes.SNPs',
            'lmiss.SNPs','hwe.SNPs','SNPQC.bim',
            ]:
            cmd = 'cat ../pops/%s/%s.%s | wc -l' %(
                pop,pop,suffix,)
            execmd(cmd,bool_exe=False,)
            l += ['%s' %(os.popen(cmd).read().strip()),]
        s_SNP += '\t'.join(l)
        s_SNP += '\n'

        ##
        ## imiss/het
        ##
        for k in d_tables.keys():
            fn = '../pops/%s/%s.table' %(pop,k,)
            if not os.path.isfile(fn):
                print fn
                stop
            if k in ['het',]: ## no header
                cmd = 'cat %s' %(fn)
            else:
                cmd = 'cat %s | awk "NR>1"' %(fn)
            d_tables[k] += os.popen(cmd).read()

        ##
        ## SNP intersection/union
        ##
        cmd = 'cat ../pops/%s/%s.hwe.SNPs ../pops/%s/%s.lmiss.SNPs' %(
            pop,pop,pop,pop,)
        cmd += ' | sort > SNPQC.SNPs'
        execmd(cmd)
        ## SNP intersection
        cmd = 'comm -12 SNPQC.SNPs.comm SNPQC.SNPs | sort -o SNPQC.SNPs.comm'
        execmd(cmd)
        ## SNP union
        cmd = 'cat SNPQC.SNPs SNPQC.SNPs.union | sort -u -o SNPQC.SNPs.union'
        execmd(cmd)

        ## continue loop over populations
        continue

    for k in d_tables.keys():
        fd = open('summary_%s.table' %(k),'w')
        fd.write(d_tables[k])
        fd.close()

    fd = open('summary_samples.table','w')
    fd.write(s_sample)
    fd.close()

    fd = open('summary_SNPs.table','w')
    fd.write(s_SNP)
    fd.close()

    print 'fgrep -f agv.sexcheck.samples ../pops/*/*.imiss.samples'
    print 'fgrep -f agv.sexcheck.samples ../pops/*/*.het.samples'

    cmd = 'cat SNPQC.SNPs.comm | wc -l'
    print cmd
    print os.popen(cmd).read()
    cmd = 'cat SNPQC.SNPs.union | wc -l'
    print cmd
    print os.popen(cmd).read()

    print 'wc -l ../pops/*/*.sampleQC.fam'
    print 'wc -l ../pops/*/*.postQC.autosomes.fam'

    print '''awk 'BEGIN{filename="None"} {if(FILENAME!=filename) {if(filename!="None") {print m,f,filename}; filename=FILENAME;m=0;f=0}; if($5==1) m++; if($5==2) f++}' pops/*/*.X.fam'''

    return


def super_MDS_PLINK(bfile,l_populations,bool_exclude,):

    fn_ld_regions = 'ldregions.SNPs'

    if os.path.isfile('%s.log' %(bfile)):
        if time.time()-os.path.getmtime('%s.log' %(bfile))<10*60:
            return
    if os.path.isfile('%s.mds' %(bfile,)):
        return

##            bool_continue = False
##            for suffix in ['mds','prune.in','genome','bed',]:
##                if os.path.isfile('%s.%s' %(bfile,suffix,)):
##                    bool_continue = True
##                    break
##            if bool_continue == True:
##                continue

    l_cmds = []

    ##
    ## 1) find common SNPs
    ##
##    if not os.path.isfile('%s.SNPs.comm' %(bfile)):
    if True:
        l_cmds += find_common_SNPs(bfile,l_populations,)

    ##
    ## 2) merge all populations
    ##
    if not os.path.isfile('%s.bed' %(bfile)):
        cmd = merge_populations(bfile,l_populations,bool_exclude,)
        l_cmds += [cmd]

    ##
    ## 3) --indep-pairwise
    ##
    if not os.path.isfile('%s.prune.in' %(bfile)):
        cmd = 'plink \\\n'
        cmd += '--bfile %s \\\n' %(bfile)
        cmd += '--out %s \\\n' %(bfile)
        ## settings
        cmd += '--indep-pairwise 50 5 0.2 \\\n'
        cmd += '--maf 0.05 \\\n'
        ## SNP exclusion
        cmd += '--exclude %s \\\n' %(fn_ld_regions)
        l_cmds += [cmd]

    ##
    ## 4) --genome
    ##
    if not os.path.isfile('%s.genome' %(bfile)):
        l_cmds += genome(bfile)

    ##
    ## 5) --cluster
    ##
    if not os.path.isfile('%s.mds' %(bfile)):
        cmd = 'plink \\\n'
        cmd += '--bfile %s \\\n' %(bfile)
        cmd += '--out %s \\\n' %(bfile)
        cmd += '--cluster \\\n'
        cmd += '--mds-plot 500 \\\n'
        cmd += '--read-genome %s.genome \\\n' %(bfile)
        cmd += '--extract %s.prune.in \\\n' %(bfile)
        l_cmds += [cmd]

    fn_sh = '%s.sh' %(bfile)
    write_shell(fn_sh,l_cmds)

    cmd = LSF_append('./%s' %(fn_sh),mem=6000,)
    execmd(cmd)

    return


def merge_populations(prefix,l_populations,bool_exclude,):

    s = ''
    for population in l_populations[1:]:
        bfile = '../pops/%s/%s.postQC.autosomes' %(population,population,)
        fn_bed = '%s.bed' %(bfile)
        if not os.path.isfile(fn_bed):
            print fn_bed, 'not found'
            sys.exit(0)
        s += '%s.bed %s.bim %s.fam\n' %(bfile,bfile,bfile,)
    fd = open('allfiles_%s.txt' %(prefix),'w')
    fd.write(s)
    fd.close()

    bfile0 = '../pops/%s/%s.postQC.autosomes' %(l_populations[0],l_populations[0],)
    cmd = 'plink \\\n'
    cmd += '--bfile %s \\\n' %(bfile0)
    cmd += '--merge-list allfiles_%s.txt \\\n' %(prefix)
    cmd += '--extract %s.SNPs.comm \\\n' %(prefix)
    cmd += '--make-bed --out %s \\\n' %(prefix)
    cmd += '--noweb --nonfounders --allow-no-sex \\\n'
    if bool_exclude == True:
        execmd('cat EIGENSOFT_heavy_4_union.SNPs discordant_SNPQC.SNPs | sort -u > exclude.SNPs')
##        cmd += '--exclude EIGENSOFT_heavy_4_union.SNPs \\\n'
        cmd += '--exclude exclude.SNPs \\\n'
##    execmd(cmd)

    return cmd


def find_common_SNPs(prefix,l_populations,):

    l_cmds = []

    population = l_populations[0]
    bfile = '../pops/%s/%s.postQC.autosomes' %(population,population,)
    cmd = "cat %s.bim | awk '{print $2}' | sort > %s.SNPs.comm" %(
        bfile,prefix,)
##    execmd(cmd)
    l_cmds += [cmd]

##    cmd = 'cat %s.SNPs.comm | wc -l' %(prefix)
##    print population, int(os.popen(cmd).read())

    for population in l_populations[1:]:
        bfile = '../pops/%s/%s.postQC.autosomes' %(population,population,)
        ## sort
        cmd = "cat %s.bim | awk '{print $2}' | sort > %s.%s.SNPs.sorted" %(
            bfile,population,prefix,)
        l_cmds += [cmd]
##        execmd(cmd)
        ## comm
        cmd = 'comm -12 %s.SNPs.comm %s.%s.SNPs.sorted > %s.SNPs.comm.out' %(
            prefix, population,prefix, prefix,)
        l_cmds += [cmd]
##        execmd(cmd)
        ##
##        os.remove('%s.SNPs.sorted' %(population))
        cmd = 'rm %s.%s.SNPs.sorted' %(population,prefix,)
        l_cmds += [cmd]
        cmd = 'mv %s.SNPs.comm.out %s.SNPs.comm' %(prefix,prefix,)
##        execmd(cmd)
        l_cmds += [cmd]
##        cmd = 'cat %s.SNPs.comm | wc -l' %(prefix)
##        print population, int(os.popen(cmd).read())
        l_cmds += ['']

    l_cmds += ['\n\n']

    return l_cmds


def MDS_plot(population,QC,bool_excl_discordant,prefix,prefix_out,):

    ## in
    if not os.path.isfile('%s.mds' %(prefix)):
        return
    ## out
    if os.path.isfile('mds.%s.png' %(prefix_out)):
        return

    ##
    ## plot
    ##
    if not os.path.isfile('%s.mds' %(prefix)):
        return

    execmd('fgrep -w -f quad.samples %s.mds > %s_quad.mds' %(prefix,prefix,))
    execmd('fgrep -w -f octo.samples %s.mds > %s_octo.mds' %(prefix,prefix,))
    line_plot = 'plot '
    if population == 'Baganda':
        ps = 3
    else:
        ps = 2
    line_plot += '"%s_quad.mds" u 4:5 ps %i pt 7 lc 1 t "quad",' %(prefix,ps,)
    line_plot += '"%s_octo.mds" u 4:5 ps 2 pt 7 lc rgb "#0000FF" t "octo",' %(prefix)
    line_plot = line_plot[:-1]

    n_samples = int(os.popen('cat %s.mds | wc -l' %(prefix)).read())-1
    n_SNPs = int(os.popen('cat %s.prune.in | wc -l' %(prefix)).read())

    if QC == 'SNPQC':
        QC_title = 'post SNP QC'
    else:
        QC_title = 'pre SNP QC'
    if bool_excl_discordant == True:
        concordance_title = 'excl discordant SNPs'
    else:
        concordance_title = 'incl discordant SNPs'
##    if bool_excl_MAF_below_5 == True:
##        maf_title = 'excl MAF<5%'
##    else:
##        maf_title = 'incl MAF<5%'

    gnuplot.scatter_plot_2d(
        '%s.mds' %(prefix),
        line_plot = line_plot,
##        column1 = 2, column2 = 3,
        xlabel = 'C1',
        ylabel = 'C2',
        title='%s, %s, %s (n_{samples}=%i, n_{SNPs}=%i)' %(
            population,QC_title,concordance_title,
            n_samples,n_SNPs,
            ),
        prefix_out='mds.%s' %(prefix_out),
##        lines_extra=lines_extra,
##        bool_remove=False,
##        bool_execute=False,
        )

    os.remove('%s_quad.mds' %(prefix))
    os.remove('%s_octo.mds' %(prefix))

    return


def init_MDS(l_populations_all):

    ##
    ## 0) define population sets
    ##
    l_Africa_South = ['Sotho','Zulu_quad','Zulu_octo',]
    l_Africa_West = [
        'Wolof','Mandinka','Jola','Fula',
        'Ga-Adangbe_quad','Ga-Adangbe_octo',
        ## Central...
        'YRI','Igbo',
        ]
    l_Africa_East = [
        'Barundi','LWK','Kalenjin','Kikuyu',
        'Ethiopia', ## SOMALI,AMHARA,OROMO
        'Baganda_quad','Baganda_octo',
        'Banyarwanda_quad','Banyarwanda_octo',
        ]

    l_BBGZ = [
        ## South
        'Zulu_quad','Zulu_octo',
        ## East
        'Baganda_quad','Baganda_octo',
        'Banyarwanda_quad','Banyarwanda_octo',
        ## West
        'Ga-Adangbe_quad','Ga-Adangbe_octo',
        ]
        
    l_Africa = l_Africa_South+l_Africa_East+l_Africa_West

    l_Africa_exclEthiopia = list(l_Africa)
    l_Africa_exclEthiopia.remove('Ethiopia')

    l_Africa_East_exclEthiopia = list(l_Africa_East)
    l_Africa_East_exclEthiopia.remove('Ethiopia')

    l_Africa_exclEthiopiaZuluSotho = list(l_Africa_East+l_Africa_West)
    l_Africa_exclEthiopiaZuluSotho.remove('Ethiopia')

    l_1000g = list(set(l_populations_all)-set(l_Africa))+['YRI','LWK',]

    d_l_popset2bfiles = {
        'Africa1000G':list(l_populations_all),
        'Africa':l_Africa,
        'BBGZ':l_BBGZ,
        'Africa_exclEthiopia':l_Africa_exclEthiopia,
        'Africa_exclEthiopiaZuluSotho':l_Africa_exclEthiopiaZuluSotho,
        'Africa_South':l_Africa_South,
        'Africa_East_exclEthiopia':l_Africa_East_exclEthiopia,
        'Africa_East':l_Africa_East,
        'Africa_West':l_Africa_West,
        '1000G':l_1000g,
        }

    return d_l_popset2bfiles


def Africa_map_init():

    if not os.path.isfile('WDB/africa-cil.txt'):
        execmd('wget http://www.evl.uic.edu/pape/data/WDB/africa.tar.gz')
        execmd('tar -xvf africa.tar.gz')

    lines_coast = []
    lines_islands = []
    lines_lakes = []
    lines_boundary = []
    lines_rivers = []

    x_min = -20
    x_max = 54
    y_min = -36
    y_max = 38

    for fn in [
        'africa-cil.txt', ## coastlines, islands, lakes
        'africa-bdy.txt', ## national boundaries
        'africa-riv.txt', ## rivers
        ]:
        fp = 'WDB/%s' %(fn)
        fd = open(fp,'r')
        lines = fd.readlines()
        fd.close()
        k = fn[-7:-4]

        for i1 in range(len(lines)):
            line = lines[i1]
            if lines[i1][:7] == 'segment':
                segment = int(lines[i1].split()[1])
                rank = int(lines[i1].split()[3])
                points = int(lines[i1].split()[5])
                ##
                ## append line breaks between segments
                ##
                if k == 'cil' and segment in [1212,1213,1214,]: continue ## Victoria Lake
                if len(lines_coast) > 0 and not lines_coast[-1] == ['\n']:
                    lines_coast += ['\n']
                if len(lines_islands) > 0 and not lines_islands[-1] == ['\n']:
                    lines_islands += ['\n']
                if len(lines_lakes) > 0 and not lines_lakes[-1] == ['\n']:
                    lines_lakes += ['\n']
                if len(lines_rivers) > 0 and not lines_rivers[-1] == ['\n']:
                    lines_rivers += ['\n']
                if k == 'bdy':
                    if len(lines_boundary) > 0 and not lines_boundary[-1] == ['\n']:
                        lines_boundary += ['\n']
                continue
##            if k == 'bdy' and rank > 2:
##                continue
            if k == 'riv' and rank > 2:
                continue
            if k == 'cil' and rank > 1:
                continue
            if k not in ['cil','riv','bdy',]:
                continue
            y = float(lines[i1].split()[0])
            x = float(lines[i1].split()[1])
            if y > y_max:
                continue
            if y < y_min:
                continue
            if x > x_max:
                continue
            if x < x_min:
                continue
    ##        if k == 'cil' and (segment > 70 and segment < 1650):
    ##            continue
    ##        if k == 'cil' and rank == 1 and (segment < 70 or segment > 1650):
    ##            continue
    ##        if k == 'riv' and segment > 1792:
    ##            continue
            if k == 'riv' and segment == 1384:
                continue
    ##        if k == 'cil' and (segment < 800 or segment > 1300):
    ##            continue

    ##        if lines[i1-1][:7] == 'segment':
    ##            print k, s, lines[i1-1]

            ##
            ##
            ##
            
            if k == 'cil' and segment in [
    ##            1211, ## Victoria Lake (West)
    ##            1212, ## Victoria Lake (North)
    ##            1213, ## Victoria Lake (East)
    ##            1214, ## Victoria Lake (South)
                163, ## Mellita (Tunisia)
                165, ## Djerba (Tunisia)
                293,295,297,298, ## Malabo, Santo Antonio,Sao Tome and Principe, Annobon Natural Reserve
                657, ## Dahlac Marine National Park
                ]:
                continue
            ## Zanzibar, Pemba Island, Island
            if k == 'cil' and segment in [
                518,519,522,584,586,
                ]:
                continue
            ## Canarian islands
            if k == 'cil' and segment in [
                190,192,193,194,195,196,197,
                205, ## Funchal
                ]:
                continue
            if k == 'cil' and segment in [
                1,2,3, ## Pantelleria,Victoria,Malta
                ]:
                continue
            ## Egypt/Suez
            if k == 'cil' and segment in [
                156,157,1366,1368,1369,1370,1371,1372,1374,
                ]:
                continue
            ## Mediterranean islands
            if k == 'cil' and segment in [1358,1406,1407,1417,1428,1433,1436,1437,]:
                continue
            ## Madagascar
            if k == 'riv' and segment in [1831,1832,1844,1845,1846,1889,1905,]:
                continue
            ## Madagascar and surrounding islands
            if k == 'cil' and segment in [457,458,459,461,462,471,703,705,707,727,729,731,733,744,759,]:
                continue
            ## Saudia Arabia and Yemen
            if k == 'riv' and segment in [
                1990,1991,2003,2006,2007,2008,2090,2091,2092,2098,
                2115,2117,2118,2119,
                2120,2121,2122,2139,2324,2331,2332,
                ]:
                continue
            ## Saudia Arabia and Yemen (and Israel, Cyprus)
            if k == 'cil' and segment in [
                1356,1357,1362,1363,1364,1365,1367,
                1373,1375,1376,1377,1378,1379,1380,1384,1489,1490,1524,1527,
                ]:
                continue

    ##        if x > 31 and x < 37 and y > -4 and y < 4 and k == 'cil' and segment > 0 and segment < 5000:
    ##            print k, segment, x ,y
    ##            pass

            ##
            ##
            ##
                
            ## national boundaries
            if k == 'bdy':
                lines_boundary += ['%s,%s\n' %(lines[i1].split()[1],lines[i1].split()[0],)]
            ## rivers
            elif k == 'riv':
                lines_rivers += ['%s,%s\n' %(lines[i1].split()[1],lines[i1].split()[0],)]
            ## coastline
            elif k == 'cil' and rank == 1:
                if segment in [
                    502, ## Bugala Island (Victoria Lake)
                    534, ## Ukerewe Island (Victoria Lake)
                    ]:
                    lines_islands += ['%s,%s\n' %(lines[i1].split()[1],lines[i1].split()[0],)]
                if segment > 500 and segment < 1600:
                    lines_lakes += ['%s,%s\n' %(lines[i1].split()[1],lines[i1].split()[0],)]
                else:
                    lines_coast += ['%s,%s\n' %(lines[i1].split()[1],lines[i1].split()[0],)]
            else:
    ##            print k, rank, segment
                continue

    fd = open('coast.txt','w')
    fd.writelines(lines_coast)
    fd.close()
    fd = open('islands.txt','w')
    fd.writelines(lines_islands)
    fd.close()
    fd = open('lakes.txt','w')
    fd.writelines(lines_lakes)
    fd.close()
    fd = open('boundary.txt','w')
    fd.writelines(lines_boundary)
    fd.close()
    fd = open('rivers.txt','w')
    fd.writelines(lines_rivers)
    fd.close()

    parse_peakware()

    ## http://stackoverflow.com/questions/1412004/reading-xml-using-python-minidom-and-iterating-over-each-node
    lines_volcanoes = []
    for country in ['Uganda','Ethiopia','Kenya','Rwanda','Tanzania',]:
        ## kml files manually downloaded from Wikipedia
        fn_xml = 'volcanoes_%s.kml' %(country)
        xmldoc = minidom.parse(fn_xml)
        placemarks = xmldoc.getElementsByTagName('Placemark')
        for placemark in placemarks:
            points = placemark.getElementsByTagName("Point")
            names = placemark.getElementsByTagName("name")
            for node_name in names:
                name = node_name.childNodes[0].nodeValue
            if name in []:
                continue
            for point in points:
                coordinates = point.getElementsByTagName('coordinates')
                for node_coordinate in coordinates:
                    coordinate = node_coordinate.childNodes[0].nodeValue
                    l = coordinate.split(',')
                    latitude = l[1]
                    longitude = l[0]
            lines_volcanoes += ['%s %s %s\n' %(longitude,latitude,name,)]

    fd = open('volcanoes.txt','w')
    fd.writelines(lines_volcanoes)
    fd.close()

    s = 'plot [%f:%f][%f:%f]' %(x_min,x_max,y_min,y_max,)
    s += '"coast.txt" w l lw 5 lc rgb "#000000",'
    s += '"lakes.txt" w filledcurve lc rgb "#808080",'
    s += '"islands.txt" w filledcurve lc rgb "#ffffff",'
    s += '"rivers.txt" w l lc rgb "#202020",'
##    s += '"volcanoes.txt" lc rgb "#202020" pt 66 ps 3,'
####    s += '"mountains.txt" lc rgb "#202020" pt 66 ps var,'
##    s += '"mountains.txt" lc rgb "#202020" pt 66 ps 3,'
####    s += '"volcanoes.txt" w labels lc rgb "#202020" pt 66 ps 3,'
    s += '"boundary.txt" w l lc 0 lw 2 lt 1,'

    return s


def init():

    fn = 'ldregions.SNPs'
    if not os.path.isfile(fn):
        print 'not found', fn
        sys.exit(0)

    l_populations = [
        ## Africa West
        'Wolof','Mandinka','Jola','Fula',
        ## Africa Central/West
        'Ga-Adangbe_quad','Ga-Adangbe_octo','YRI','Igbo',
        ## Africa East
        'Barundi', 'Banyarwanda_quad', 'Banyarwanda_octo',
        'Baganda_quad', 'Baganda_octo', 'LWK', 'Kalenjin', 'Kikuyu',
        'Ethiopia', ## 'Amhara', 'Somali', 'Oromo',
        ## Africa South
        'Sotho','Zulu_quad','Zulu_octo',
        ## Out of Africa
        'ASW','ACB',
        ## Europe
        'IBS','TSI','CEU','GBR','FIN',
        ## Asia
        'CDX','CHB','CHS','KHV','JPT',
        ## North America
        'GIH',
        ## South America
        'PEL','MXL','CLM','PUR',
##        'CHD', ## 1 sample
##        'MKK', ## 31 samples
        ]

    d_pops2coords = {
        ## http://en.wikipedia.org/wiki/Jola_people
        ## http://en.wikipedia.org/wiki/Casamance
        ## http://en.wikipedia.org/wiki/Ziguinchor
        'Jola':[12.561944,-16.283889,],
        ## http://en.wikipedia.org/wiki/Wolof_people
        ## http://en.wikipedia.org/wiki/Senegal
        ## http://en.wikipedia.org/wiki/Dakar
        'Wolof':[14.692778,-17.446667,],
        ## http://en.wikipedia.org/wiki/Baganda
        ## http://en.wikipedia.org/wiki/Kampala
        'Baganda':[0.313611,32.581111,],
        ## http://en.wikipedia.org/wiki/Barundi
        ## http://en.wikipedia.org/wiki/Burundi
        'Barundi':[-3.5,30,],
        ## http://en.wikipedia.org/wiki/Kalenjin_people
        ## http://en.wikipedia.org/wiki/Rift_Valley_Province
        ## http://en.wikipedia.org/wiki/Nakuru
        'Kalenjin':[-0.3,36.066667,],
        ## http://en.wikipedia.org/wiki/Kikuyu_people
        ## http://en.wikipedia.org/wiki/Mount_Kenya
        'Kikuyu':[-0.150833,37.3075,],
        ## http://en.wikipedia.org/wiki/Zulu_people
        ## http://en.wikipedia.org/wiki/KwaZulu-Natal
        ## http://en.wikipedia.org/wiki/Pietermaritzburg
        'Zulu':[-29.616667,30.383333,],
        ## http://en.wikipedia.org/wiki/Sotho_people
        ## http://en.wikipedia.org/wiki/Lesotho
        ## http://en.wikipedia.org/wiki/Maseru
        'Sotho':[-29.31,27.48,],
        ## http://en.wikipedia.org/wiki/Igbo_people#Nigeria
        ## http://en.wikipedia.org/wiki/Enugu
        'Igbo':[6.452667,7.510333,],
        ## http://en.wikipedia.org/wiki/Ga-Adangbe_people
        ## http://en.wikipedia.org/wiki/Accra
        'Ga-Adangbe':[5.55,-0.2,],
        ## http://en.wikipedia.org/wiki/Mandinka_people
        ## http://en.wikipedia.org/wiki/The_Gambia
        'Mandinka':[13.466667,-16.6,],
        ## http://en.wikipedia.org/wiki/Fula_people
        ## http://en.wikipedia.org/wiki/Guinea
        ## http://en.wikipedia.org/wiki/Fouta_Djallon
        ## http://en.wikipedia.org/wiki/Lab%C3%A9
        'Fula':[11.316667,-12.283333,],
        ## Ethiopia
        ## http://en.wikipedia.org/wiki/Somali_people
        ## http://en.wikipedia.org/wiki/Amhara_people
        ## http://en.wikipedia.org/wiki/Oromo_people
        ## http://en.wikipedia.org/wiki/Oromia_Region
        ## http://en.wikipedia.org/wiki/Addis_Ababa
        'Ethiopia':[9.03,38.74,],
        ## http://en.wikipedia.org/wiki/Luhya_people
        ## http://en.wikipedia.org/wiki/Kitale (too far North according to Louise)
        ## http://en.wikipedia.org/wiki/Western_Province_(Kenya)
##        'LWK':[1.016667,35,],
        'LWK':[0.5,34.583333,],
        ## http://en.wikipedia.org/wiki/Yoruba_people#Demographics
        ## http://en.wikipedia.org/wiki/Lagos_State,_Nigeria
        ## http://en.wikipedia.org/wiki/Lagos
        'YRI':[6.453056,3.395833,],
        ## http://en.wikipedia.org/wiki/Banyarwanda
        ## http://en.wikipedia.org/wiki/Umutara_Province
        ## http://en.wikipedia.org/wiki/Nyagatare
        'Banyarwanda':[-1.3,30.325,],
        }

    return l_populations, d_pops2coords


if __name__ == '__main__':
    main()
