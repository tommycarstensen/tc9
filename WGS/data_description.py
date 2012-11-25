#!/software/bin/python

## T. Carstensen (tc9), M.S. Sandhu (ms23), D. Gurdasani (dg11)
## Wellcome Trust Sanger Institute, 2012

import os

def main():

    path = '/lustre/scratch107/projects/uganda/users/tc9'

    l_steps = [
        'UnifiedGenotyper',
        'CombineVariants',
        'ApplyRecalibration',

        'ProduceBeagleInput',
        'BEAGLE', ## not a GATK step
        'BeagleOutputToVCF',
        ## Imputation, IMPUTE2
        'IMPUTE2', ## not a GATK step
        ]

##    check_CombineVariants(path)

    check_ApplyRecalibration(path)

    return


def check_ApplyRecalibration(path):

    fp_next = os.path.join(path,'out_GATK/sep','ApplyRecalibration.recalibrated.filtered.18.vcf')
    fd_next = open(fp_next,'r')
    fp_prev = os.path.join(path,'out_GATK/sep','CombineVariants.18.vcf')
    fd_prev = open(fp_prev,'r')
    for line_prev in fd_prev:
        if line_prev[0] == '#':
            line_header_prev = line_prev
            continue
        n_columns = len(line_header_prev.split('\t'))
        for line_next in fd_next:
            if line_next[0] != '#':
                break
            else:
                line_header_next = line_next
        if line_prev != line_next:
            l_prev = line_prev.split('\t')
            l_next = line_next.split('\t')
            if len(l_prev) != len(l_next):
                print len(l_prev), len(l_next)
                stop
            for i in xrange(n_columns):
                if i == 6: ## INFO
                    if (
                        l_prev[i] != 'PASS'
                        or
                        (
                            l_next[i][:23] != 'TruthSensitivityTranche'
                            and
                            l_next[i] != 'PASS'
                            )
                        ):
                        print l_prev[i]
                        print l_next[i]
                        print line_header_prev
                        print line_header_next
                        stop6
                elif i == 7: ## FORMAT
                    suffix_prev = l_prev[i][l_prev[i].rindex(';set=variant'):]
                    suffix_next = l_next[i][l_next[i].rindex(';VQSLOD='):]
                    if not (
                        'culprit=MQ' in suffix_next
                        or
                        'culprit=DP' in suffix_next
                        or
                        'culprit=QD' in suffix_next
                        or
                        'culprit=FS' in suffix_next
                        or
                        'culprit=HaplotypeScore' in suffix_next
                        or
                        'culprit=ReadPosRankSum' in suffix_next
##                        or
##                        'PASS' == l_next[i]
                        ):
                        print 'a'+l_next[i]+'a'
                        print suffix_next
                        stop
                    if l_prev[i][:-len(suffix_prev)] != l_next[i][:-len(suffix_next)]:
                        print l_prev[i][170:]+'a'
                        print l_next[i][170:]+'b'
                        print 
                        print l_prev[i][:-len(suffix_prev)]
                        print l_next[i][:-len(suffix_next)]
                        for j in range(len(l_prev[i])):
                            if l_prev[i][:j] != l_next[i][:j]:
                                print l_prev[i][:j]
                                print l_next[i][:j]
                                stop
                        stop7
                else:
                    if l_prev[i] != l_next[i]:
                        print line_header
                        print line_next
                        print l_prev[:i+1]
                        print l_next[:i+1]
                        print l_prev[i]
                        print l_next[i]
                        for x in xrange(1000):
                            if l_prev[i][:x] == l_next[i][:x]:
                                continue
                            print l_prev[i][:x]
                            print l_next[i][:x]
                            break
                        print i, x
                        stop3
        else:
            print line_prev
            print line_next
            stop
##    print line_next
##    print line_prev
    fd_prev.close()
    fd_next.close()
        
    return


def check_CombineVariants(path):

    fp = os.path.join(path,'out_GATK/sep','CombineVariants.18.vcf')
    fd_CV = open(fp,'r')
    for x in xrange(1,100):
        fp = os.path.join(path,'out_GATK','UnifiedGenotyper.18.vcf.%s' %(x))
        if not os.path.isfile(fp):
            break
        fd_UG = open(fp,'r')
        for line_UG in fd_UG:
            if line_UG[0] == '#':
                line_UG_prev = line_UG
                n_columns = len(line_UG.split('\t'))
                continue
            ## next line in next VCF
            for line_CV in fd_CV:
                if line_CV[0] != '#':
                    break
                else:
                    line_CV_prev = line_CV
            if line_UG != line_CV:
                l_UG = line_UG.split('\t')
                l_CV = line_CV.split('\t')
                if len(l_UG) != len(l_CV):
                    print len(l_UG), len(l_CV)
                    stop
                for i in xrange(n_columns):
                    if i == 6:
                        if l_UG[i] != '.' or l_CV[i] != 'PASS':
                            stop6
                    elif i == 7:
                        if x == 1:
                            suffix = ';set=variant'
                        else:
                            suffix = ';set=variant%i' %(x)
                        if l_UG[i]+suffix != l_CV[i]:
                            print l_UG[i][170:]
                            print l_CV[i][170:]
                            stop7
                    else:
                        if l_UG[i] != l_CV[i]:
                            print line_UG
                            print line_CV
                            print l_UG[:i+1]
                            print l_CV[:i+1]
                            print l_UG[i]
                            print l_CV[i]
                            for x in xrange(1000):
                                if l_UG[i][:x] == l_CV[i][:x]:
                                    continue
                                print l_UG[i][:x]
                                print l_CV[i][:x]
                                break
                            print i, x
                            stop3
            else:
                print line_UG
                print line_CV
                stop
        fd_UG.close()
        
    return


if __name__ == '__main__':
    main()
