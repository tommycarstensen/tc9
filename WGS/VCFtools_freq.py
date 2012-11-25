#!/software/bin/python

## T. Carstensen (tc9), M.S. Sandhu (ms23), D. Gurdasani (dg11)
## Wellcome Trust Sanger Institute, 2012

import os
import sys
sys.path.append('/nfs/users/nfs_t/tc9/github/tc9/misc')
import gnuplot

def main():

    for chromosome in range(1,22+1,)+['X','Y',]:
        fp_in = 'out_GATK/sep/ApplyRecalibration.recalibrated.filtered.%s.vcf' %(chromosome)
        fp_out = 'out_VCFtools/freq%s' %(chromosome)
        if os.path.isfile('%s.frq' %(fp_out)):
            continue
        s = 'bsub \
        -M4000000 -R\'select[mem>4000] rusage[mem=4000]\' \
        vcftools \
        --vcf %s \
        --freq \
        --out %s \
        ' %(fp_in, fp_out,)
        os.system(s)

    for chromosome in range(1,22+1,)+['X','Y',]:
        fp_out = 'out_VCFtools/freq%s.frq' %(chromosome)
##        if chromosome <= 11:
##            continue
        print fp_out
        fd = open(fp_out,'r')
        for line in fd: break
        l_pos = []
        l_MAF = []
        for line in fd:
            l = line.split()
            pos = float(l[1])/10**6
            MAF = l[-1][2:]
            l_pos += [pos]
            l_MAF += [MAF]
##            lines += ['%s %s\n' %(pos,MAF,)]
##        fd = open('gnuplot%s.data','w')
##        fd.writelines(lines)
##        fd.close()
        gnuplot.scatter_plot_2d(
            'MAF%s' %(chromosome),
            l_pos, l_MAF,
            xlabel = 'pos (Mbp)',
            ylabel = 'MAF',
            )

    return

if __name__ == '__main__':
    main()
