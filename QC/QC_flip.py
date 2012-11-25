#!/software/bin/python

## T. Carstensen (tc9), M.S. Sandhu (ms23), D. Gurdasani (dg11)
## Wellcome Trust Sanger Institute, 2012

import os, sys

def main():

    bfile = 'omni2.5-8_20120809_gwa_uganda_gtu'
    bfile = '/nfs/african_diversity/data/1KG/genotypes/release20110527/working_data/omni25_b37_bed_autosomes/1kg_r20110527_omni25_b37_autosomes'
    bfile = sys.argv[sys.argv.index('--bfile')+1]

    ##
    ## assign strand file
    ##
    n_SNPs = int(os.popen('cat %s.bim | wc -l' %(bfile)).read())
    if n_SNPs == 2450000:
        fn          = 'HumanOmni2.5M-b37-strand-v2'
        strand = 'HumanOmni2.5M-b37-v2.strand'
    elif n_SNPs == 2379855:
        fn          = 'HumanOmni2.5-8v1_A-b37-strand'
        strand = 'HumanOmni2.5-8v1_A-b37.strand'
    elif bfile == '/nfs/african_diversity/data/1KG/genotypes/release20110527/working_data/omni25_b37_bed_autosomes/1kg_r20110527_omni25_b37_autosomes':
        fn = 'HumanOmni2.5M-b37-strand-v2'
        strand = 'HumanOmni2.5M-b37-v2.strand'
    else:
        print bfile
        print n_SNPs
        stop
    if not os.path.isfile(strand):
        os.system('wget http://www.well.ox.ac.uk/~wrayner/strand/%s.zip' %(fn))
        os.system('unzip %s.zip' %(fn))
        os.remove('%s.zip' %(fn))

    ##
    ## make sure output is to cwd
    ##
    if '/' in bfile:
        out = bfile[bfile.rindex('/')+1:]
    else:
        out = bfile
        for suffix in ['bed','fam','bim',]:
                os.system('mv %s.%s %s.preflip' %(bfile,suffix,out,))

    cmd = 'echo "job initiated" | mail -s "job $LSB_JOBID initiated flip %s" tc9@sanger.ac.uk\n\n' %(bfile,)

    ##
    ## write list of SNPs to be flipped
    ##
    cmd += '''awk '{if($5=="-") print $1}' %s.strand > %s.flip\n\n''' %(strand,bfile,)

    cmd += 'plink \\\n'
    cmd += '--bfile %s.preflip \\\n' %(bfile)
    cmd += '--flip %s.flip --make-bed \\\n' %(bfile)
    cmd += '--out %s.flipped \\\n' %(out)
    cmd += '--noweb \\\n'
    cmd += '--allow-no-sex \\\n'
    cmd += '--nonfounders \\\n'

    cmd += '\necho "job finished" | mail -s "job $LSB_JOBID finished flip %s" tc9@sanger.ac.uk' %(bfile)

    fn_out = '%s.flip.sh' %(out)
    fd = open(fn_out,'w')
    fd.write(cmd)
    fd.close()

    os.system('chmod +x %s' %(fn_out))

    bsub = "bsub -M8000000 -R'select[mem>8000] rusage[mem=8000]' "
    bsub += '-P uganda_gwas '
    bsub += '-q normal '
    bsub += "-J'%s' " %(bfile)
    bsub += '%s' %(fn_out)

    os.system(bsub)

##    self.fp_flip = '%s.flip' %(bfile)
##
##        ##
##        ## flip additional SNPs in 1000g to match those of current data
##        ## should not be part of __init__ function!!!
##        ##
##        self.fp_flip1000g = 'flip1000g.txt'
##        if not os.path.isfile(self.fp_flip1000g):
##            fd1 = open('%s.bim' %(self.bed_1000g),'r')
##            fd2 = open('%s.flipped.bim' %(bfile),'r')
##            bool_read1 = True
##            bool_read2 = True
##            d_flip = {'A':'T','G':'C','T':'A','C':'G',}
##        ##            os.system('cp %s %s' %(self.fp_flip,self.fp_flip1000g))
##        ##            fd = open(self.fp_flip1000g,'w')
##            os.system('cp exclude.missnp %s' %(self.fp_flip1000g))
##            fd = open(self.fp_flip1000g,'a')
##            while True:
##                if bool_read1 == True:
##                    line1 = fd1.readline()
##                    if line1 == '':
##                        break
##                    l1 = line1.split()
##                    chrom1 = int(l1[0])
##                    pos1 = int(l1[3])
##                    bool_read1 = False
##                if bool_read2 == True:
##                    line2 = fd2.readline()
##                    if line2 == '':
##                        break
##                    l2 = line2.split()
##                    chrom2 = int(l2[0])
##                    pos2 = int(l2[3])
##                    bool_read2 = False
##                if chrom1 < chrom2:
##                    bool_read1 = True
##                elif chrom1 > chrom2:
##                    bool_read2 = True
##                elif pos1 == pos2 or l1[1] == l2[1]:
##                    bool_read1 = True
##                    bool_read2 = True
##                    if l1[1] == 'rs10241053':
##                        stop2
##                    if l1[4] == l2[4] and l1[5] == l2[5]:
##                        continue
##                    elif l1[4] == l2[5] and l1[5] == l2[4]:
##                        continue
##        ##                    elif l1[4] == 'T' and l1[5] == 'C' and l2[4] == 'A' and l2[5] == 'G':
##        ##                        fd.write('%s\n' %(l1[1]))
##        ##                    elif l1[4] == 'T' and l1[5] == 'C' and l2[4] == 'G' and l2[5] == 'A':
##        ##                        fd.write('%s\n' %(l1[1]))
##        ##                    elif l1[4] == 'T' and l1[5] == 'G' and l2[4] == 'C' and l2[5] == 'A':
##        ##                        fd.write('%s\n' %(l1[1]))
##        ##                    elif l1[4] == 'T' and l1[5] == 'G' and l2[4] == 'A' and l2[5] == 'C':
##        ##                        fd.write('%s\n' %(l1[1]))
##        ##                    elif l1[4] == 'G' and l1[5] == 'T' and l2[4] == 'A' and l2[5] == 'C':
##        ##                        fd.write('%s\n' %(l1[1]))
##        ##                    elif l1[4] == 'G' and l1[5] == 'T' and l2[4] == 'C' and l2[5] == 'A':
##        ##                        fd.write('%s\n' %(l1[1]))
##        ##                    elif l1[4] == 'C' and l1[5] == 'T' and l2[4] == 'A' and l2[5] == 'G':
##        ##                        fd.write('%s\n' %(l1[1]))
##        ##                    elif l1[4] == 'C' and l1[5] == 'T' and l2[4] == 'G' and l2[5] == 'A':
##        ##                        fd.write('%s\n' %(l1[1]))
##                    elif l1[4] == '0' and l1[5] in l2[4:6]:
##                        continue
##                    elif l2[4] == '0' and l2[5] in l1[4:6]:
##                        continue
##                    elif l2[4] == '0' and l2[5] not in l1[4:6]:
##                        fd.write('%s\n' %(l1[1]))
##                    elif l1[4] == '0' and l1[5] not in l2[4:6]:
##                        fd.write('%s\n' %(l1[1]))
##                    elif d_flip[l1[4]] == l2[4] and d_flip[l1[5]] == l2[5]:
##                        fd.write('%s\n' %(l1[1]))
##                    elif d_flip[l1[4]] == l2[5] and d_flip[l1[5]] == l2[4]:
##                        fd.write('%s\n' %(l1[1]))
##                    ## SNP IDs different
##                    elif l1[1] != l2[1]:
##                        print
##                        print l1
##                        print l2
##                        continue
##                    else:
##                        print
##                        print l1,
##                        print l2
##                        stop
##        ##                    if l1[1] == l2[1] and l1[4] != '0' and l2[4] != '0':
##        ##                        print l1
##        ##                        print l2
##        ##                        stop
##                elif pos1 < pos2:
##                    bool_read1 = True
##                    if l1[1] == l2[1]:
##                        print '%2i %9i %9i %10s %10s' %(chrom1, pos1, pos2, l1[1], l2[1],)
##                else:
##                    bool_read2 = True
##                    if l1[1] == l2[1]:
##                        print '%2i %9i %9i %10s %10s' %(chrom2, pos1, pos2, l1[1], l2[1],)
##            fd.close()
##            os.system('sort -u flip1000g.txt -o flip1000g.txt')

## http://pngu.mgh.harvard.edu/~purcell/plink/dataman.shtml#flip






##
##    sh = ''
##
##    ## http://pngu.mgh.harvard.edu/~purcell/plink/dataman.shtml#updatemap
##
##
##    ## bsub -M8000000 -R'select[mem>8000] rusage[mem=8000]'
##    ## -o update_build4.out -e update_build4.err
##    ## ./update_build.sh
##    ## /lustre/scratch107/projects/agv/data/omni2.5-4_20120904/omni2.5-4_20120904_agv_gtu_aaa
##    ## HumanOmni2.5M-b37-v2.strand
##    ## omni2.5-4_20120904_agv_gtu_aaa
##
##    fp_prefix_in = '/lustre/scratch107/projects/agv/data/omni2.5-4_20120904/omni2.5-4_20120904_agv_gtu_aaa'
##    fn_strand = 'HumanOmni2.5M-b37-v2.strand'
##
##    fp_prefix_in = '/lustre/scratch107/projects/agv/data/omni2.5-8_20120910/omni2.5-8_agv_20120910_gtu'
##    fn_strand = 'HumanOmni2.5-8v1_A-b37.strand'
##
##    if '/' in fp_prefix_in:
##        fn_prefix_in = fp_prefix_in[fp_prefix_in.rindex('/')+1:]
##        fn_prefix_out = fn_prefix_in
##    else:
##        fn_prefix_out = fp_prefix_in
##
####    if not '--prefix' in sys.argv:
####        print '--prefix missing'
####        sys.exit()
####    if not '--strand' in sys.argv:
####        print '--strand missing'
####        sys.exit()
####    fp_prefix_in = sys.argv[sys.argv.index('--prefix-in')+1]
####    fp_strand = sys.argv[sys.argv.index('--strand')+1]
##
##    ## format of chr file is two columns: rsid, chr
##    ## format of SNP dic is two columns: old rsid, new rsid
##    ## format of pos file is two columns: rsid, pos
##
##    d_cols = {
##        'update-chr':2,
##        'update-pos':3,
##        }
##
##    for plink_cmd,col in d_cols.items():
##        cmd = 'cut -f 1,%i %s > %s_%s.txt' %(
##            d_cols[plink_cmd], fn_strand, fn_strand, plink_cmd,
##            )
####        print cmd
##        sh += cmd+'\n\n'
####        os.system(cmd)
##    cmd = '''awk '{if ($5=="-") print $1}' %s > %s_%s.txt''' %(fn_strand,fn_strand,'flip',)
####    print cmd
##    sh += cmd+'\n\n'
####    os.system(cmd)
##
##    l_plink_cmds = [
##        ## update the chromosome of all variants to build 37
##        'update-chr',
####        'update-name',
##        ## update the position of all variants to build 37
##        'update-pos',
##        ## flip the orientation of all variants not aligned to the forward DNA strand
##        'flip',
##        ]
##
##    for i_plink_cmd in xrange(len(l_plink_cmds)):
##        plink_cmd = l_plink_cmds[i_plink_cmd]
##        cmd = 'plink \\\n'
##        if i_plink_cmd == 0:
##            cmd += '--bfile %s \\\n' %(fp_prefix_in)
##        else:
##            cmd += '--bfile %s_%s \\\n' %(fp_prefix_in,l_plink_cmds[i_plink_cmd-1])
##        cmd += '--make-bed \\\n'
##        cmd += '--out %s_%s \\\n' %(fn_prefix_in,plink_cmd,)
##        if plink_cmd != 'flip':
##            cmd += '--update-map %s.txt \\\n' %(plink_cmd,)
##        if plink_cmd != 'update-pos':
##            cmd += '--%s \\\n' %(plink_cmd)
##        cmd += '--noweb \\\n'
##        cmd += '--allow-no-sex \\\n'
##        cmd += '--non-founders \\\n'
##
##        sh += cmd+'\n'
##
##    ##
##    ## clean up at the end
##    ##
##    for plink_cmd in l_plink_cmds:
##        sh += '\nrm %s_%s.txt' %(fn_strand,plink_cmd,)
####        s += '\nrm %s.hh' %(fn_prefix_out)
##        if plink_cmd != 'flip':
##            for suffix in ['bed','bim','fam',]:
##                sh += '\nrm %s_%s.%s' %(fn_prefix_in,plink_cmd,suffix,)
##
##    fn_out = 'flip_%s.sh' %(fn_prefix_out)
##    fd = open(fn_out,'w')
##    fd.write(sh)
##    fd.close()
##    os.system('chmod +x %s' %(fn_out))
##
##    print '''bsub -M8000000 -R'select[mem>8000] rusage[mem=8000]' ./flip_%s.sh''' %(fn_prefix_out)

    return


if __name__ == '__main__':
    main()
