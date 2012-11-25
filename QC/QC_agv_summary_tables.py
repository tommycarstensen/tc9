#!/software/bin/python

## T. Carstensen (tc9), M.S. Sandhu (ms23), D. Gurdasani (dg11)
## Wellcome Trust Sanger Institute, 2012

import os, sys
sys.path.append(os.path.dirname(sys.argv[0]))
import QC

def main():

    l_populations = [
        u'Muganda', u'KIKUYU', u'Mandinka', u'ZULU', u'IBO', u'KALENJIN', u'Fula', u'Murundi', u'Munyarwanda', u'AMHARA', u'Sotho', u'Wolloff', u'SOMALI', u'Jola', u'OROMO', u'GA-ADANGBE'
        ]
    l_populations.sort()

    tables(l_populations)

    plots(l_populations)

    return


def plots(l_populations):

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
            cmd = 'head -1 %s.%s > agv.%s' %(l_populations[0],suffix,suffix,)
            os.system(cmd)
        for population in l_populations:
            ## no header
            if suffix in ['fam','sampleQC.samples',]:
                cmd = "cat %s.%s >> agv.%s" %(population,suffix,suffix,)
            ## header
            else:
                cmd = "sed '1d' %s.%s >> agv.%s" %(population,suffix,suffix,)
            os.system(cmd)
            continue

        continue

    instanceQC = QC.main()
##    instanceQC.plink_plots('agv',i_wait=0)

    ## samples
    instanceQC.histogram_imiss('agv',)
    instanceQC.histogram_het('agv',bool_with_stddev=False,)
    instanceQC.histogram_genome('agv',)
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

    s_sample = ''
    s_SNP = ''
    for population in l_populations:
        s_sample += '%s\t' %(population)
        s_SNP += '%s\t' %(population)

        l = []
        for suffix in ['imiss','het','sexcheck','sampleQC','genome0.20']:
            cmd = 'cat %s.%s.samples | wc -l' %(population,suffix,)
            l += ['%s' %(os.popen(cmd).read().strip()),]
        s_sample += '\t'.join(l)

        l = []
        cmd = 'cat %s.bim | wc -l' %(population,)
        l += ['%s' %(os.popen(cmd).read().strip()),]
        for suffix in [
            'position','miss','duplicates','X','autosomes','lmiss','hwe',
            ]:
            cmd = 'cat %s.%s.SNPs | wc -l' %(population,suffix,)
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
