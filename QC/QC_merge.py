#!/software/bin/python

## T. Carstensen (tc9), M.S. Sandhu (ms23), D. Gurdasani (dg11)
## Wellcome Trust Sanger Institute, 2012

import os, sys

def main():

    ## merge quad and octo data for each population

    l_bfiles = [
        'omni2.5-4_20120904_agv_gtu_aaa',
        'omni2.5-8_agv_20120910_gtu',
        ]

    if not os.path.isfile('%s.bim' %(l_bfiles[0])):
        print 'wrong dir'
        return

    l_populations = [
        u'Muganda', u'KIKUYU', u'Mandinka', u'ZULU', u'IBO', u'KALENJIN', u'Fula', u'Murundi', u'Munyarwanda', u'AMHARA', u'Sotho', u'Wolloff', u'SOMALI', u'Jola', u'OROMO', u'GA-ADANGBE'
        ]

    extract_and_translate(l_bfiles,) 

    extract_and_translate_cleanup(l_bfiles,)

    bmerge(l_bfiles,l_populations,)

####    os.system('rm quad2octo.dic extract.SNPs')

    return


def extract_and_translate_cleanup(l_bfiles,):

    ##
    ## 8) clean up
    ##
    for bfile in l_bfiles:
        os.remove('%s.bim.chrpos.sorted' %(bfile))
        os.remove('%s.bim.chrpos.sorted.concatenated' %(bfile))
        os.remove('%s.bim.sorted' %(bfile))
        os.remove('%s.bim.sorted.filtered' %(bfile))
        os.remove('%s.bim.resorted.filtered.concatenated' %(bfile))
        ## all SNPs
        os.remove('%s.sorted.SNPs' %(bfile))
        ## exclusion
        os.remove('%s.duplicates.SNPs' %(bfile))
        os.remove('%s.mismatch.SNPs' %(bfile))
        ## inclusion
        os.remove('%s.inclusion.SNPs' %(bfile))
    ## exclusion
    os.remove('common_name_different_chrpos.SNPs')
    os.remove('mismatch.SNPs')
    os.remove('exclusion.SNPs')

    return


def extract_and_translate(l_bfiles,):

    '''
this function creates the files
quad2octo.dic and extract.SNPs
used by the function bmerge
'''

    col_bim_chr = 1
    col_bim_name = 2
    col_bim_pos = 4

    ##
    ## 1) sort
    ##
    sort(
        l_bfiles, '',
        'sort -k%i,%i' %(col_bim_name, col_bim_name,),
        )
    sort(
        l_bfiles, '.chrpos',
        'sort -k%i,%i -k%i,%i' %(
            col_bim_chr,col_bim_chr, col_bim_pos,col_bim_pos,
            ),
        )

    ##
    ## 2) identify duplicates (chromosome + position)
    ##
    identify_duplicates(l_bfiles,col_bim_chr,col_bim_pos,col_bim_name,)

    ##
    ## 3) identify allele mismatches not due to strand flipping
    ##
    ##
    ## 3a) concatenate chr+pos into one column
    ##
    identify_allele_mismatches(l_bfiles,col_bim_chr,col_bim_pos,col_bim_name,)

    ##
    ## 4) identify common rsID, different position
    ##
    ## 4a) join by name (i.e. identify common rsIDs)
    cmd = 'join -1 %i -2 %i -o 0 1.%i 2.%i 1.%i 2.%i ' %(
        col_bim_name, col_bim_name,
        col_bim_chr, col_bim_chr, col_bim_pos, col_bim_pos,
        )
    cmd += '%s.bim.sorted %s.bim.sorted' %(
        l_bfiles[0], l_bfiles[1],
        )
    ## 4b) identify different positions
    cmd += " | awk '{if ($2!=$3 || $4!=$5) print $1}'"
    ## sort before comm
    cmd += ' | sort'
    cmd += ' > common_name_different_chrpos.SNPs'
    print cmd
    os.system(cmd)

    ##
    ## 5) concatenate exclusion lists from steps 2,3,4
    ##
    cmd = 'cat'
    ## step 2
    cmd += ' %s.duplicates.SNPs %s.duplicates.SNPs ' %(
        l_bfiles[0], l_bfiles[1],
        )
    ## step 3
    cmd += ' %s.mismatch.SNPs %s.mismatch.SNPs' %(l_bfiles[0], l_bfiles[1],)
    ## step 4
    cmd += ' common_name_different_chrpos.SNPs'
    ## sort and redirect stdout
    cmd += ' | sort -u > exclusion.SNPs'
    print cmd
    os.system(cmd)

    ##
    ## 6) exclude I) duplicates and II) common rsID but different position
    ##
    ## N.B. the duplicates need to be excluded before quad to octo translation,
    ## otherwise there is a risk of a random and incorrect translation
    ##
    ##
    ## 6a) write name column to file
    ##
    for bfile in l_bfiles:
        cmd = "cat %s.bim.sorted | awk '{print $%i}' > %s.sorted.SNPs" %(
            bfile, col_bim_name, bfile,
            )
        print cmd
        os.system(cmd)
        continue

    ##
    ## 6b) generate inclusion lists
    ##
    for bfile in l_bfiles:
        cmd = 'comm -23 %s.sorted.SNPs exclusion.SNPs > %s.inclusion.SNPs' %(
            bfile, bfile,
            )
        print cmd
        os.system(cmd)

    ##
    ## 6c) join inclusion lists and sorted bim files;
    ## i.e. exclude SNPs not in inclusion list
    ##
    for bfile in l_bfiles:
        cmd = 'join -1 1 -2 %i' %(col_bim_name,)
        cmd += ' %s.inclusion.SNPs %s.bim.sorted -o 2.1,2.2,2.3,2.4,2.5,2.6' %(
            bfile, bfile,
            )
        cmd += ' > %s.bim.sorted.filtered' %(bfile,)
        print cmd
        os.system(cmd)

    ##
    ## 7) generate dictionary / translation table for remaining SNP names
    ##
    ## N.B. this needs to be done for unique rsIDs *after* exclusion of duplicates
    ## to avoid incorrect mappings by join
    ##
    ##
    ## 7a) concatenate chr+pos into one column
    ##
    for bfile in l_bfiles:
        cmd = 'cat %s.bim.sorted.filtered' %(bfile)
        cmd += ' | sort -k%i,%i -k%i,%i' %(
            col_bim_chr,col_bim_chr, col_bim_pos,col_bim_pos,
            )
        cmd += ''' | awk '{print $%i":"$%i,$%i}' ''' %(
            col_bim_chr,col_bim_pos,col_bim_name,
            )
        cmd += ' > %s.bim.resorted.filtered.concatenated' %(bfile,)
        print cmd
        os.system(cmd)
    ##
    ## 7b) join by concatenated chr+pos
    ##
    cmd = 'join -1 %i -2 %i ' %(1,1,)
    cmd += '-o 1.%i 2.%i ' %(2,2,)
    cmd += ' %s.bim.resorted.filtered.concatenated' %(l_bfiles[0],)
    cmd += ' %s.bim.resorted.filtered.concatenated' %(l_bfiles[1],)
    cmd += ' > quad2octo.dic'
    print cmd
    os.system(cmd)

    cmd = "cat quad2octo.dic | awk '{print $2}' > extract.SNPs"
    print cmd
    os.system(cmd)

    return


def identify_allele_mismatches(l_bfiles,col_bim_chr,col_bim_pos,col_bim_name,):

    ## concatenate chromosome and position
    for bfile in l_bfiles:
        cmd = 'cat %s.bim.chrpos.sorted' %(bfile)
        cmd += ' | sort -k%i,%i -k%i,%i' %(
            col_bim_chr,col_bim_chr, col_bim_pos,col_bim_pos,
            )
        cmd += ''' | awk '{print $0,$%i":"$%i}' ''' %(col_bim_chr, col_bim_pos,)
        cmd += ' > %s.bim.chrpos.sorted.concatenated''' %(bfile,)
        print cmd
        os.system(cmd)
    ## join on concatenated field/column
    cmd = 'join -1 %i -2 %i' %(7,7,)
    cmd += ' -o 1.%i,2.%i,1.%i,1.%i,2.%i,2.%i' %(
        col_bim_name,col_bim_name,5,6,5,6,
        )
    cmd += ' %s.bim.chrpos.sorted.concatenated' %(l_bfiles[0])
    cmd += ' %s.bim.chrpos.sorted.concatenated' %(l_bfiles[1])
    cmd += ''' | awk '{if (!('''
    cmd += '''(($3==0 || $5==0 || $3==$5) && ($4==0 || $6==0 || $4==$6))'''
    cmd += ' || '
    cmd += '''(($3==0 || $6==0 || $3==$6) && ($4==0 || $5==0 || $4==$5))'''
    cmd += ''')) print $1,$2}' '''
    cmd += ' > mismatch.SNPs'
    print cmd
    os.system(cmd)

    ## print columns of mismatch.SNPs to individual mismatch.SNPs files
    for i in range(len(l_bfiles)):
        bfile = l_bfiles[i]
        cmd = "cat mismatch.SNPs | awk '{print $%i}' > %s.mismatch.SNPs" %(
            i+1, l_bfiles[i],
            )
        print cmd
        os.system(cmd)

    return


def identify_duplicates(l_bfiles,col_bim_chr,col_bim_pos,col_bim_name,):

    for bfile in l_bfiles:
        print 'identify duplicates %s' %(bfile)
        fd_in = open('%s.bim.chrpos.sorted' %(bfile,),'r')
        fd_out = open('%s.duplicates.SNPs' %(bfile,),'w')
        for line in fd_in:
            line_prev = line
            break
        for line in fd_in:
            l = line.split()
            l_prev = line_prev.split()
            ## duplicate
            if (
                l_prev[col_bim_chr-1] == l[col_bim_chr-1]
                and
                l_prev[col_bim_pos-1] == l[col_bim_pos-1]
                ):
                if (
                    l[4] not in ['0','A','C','G','T',]
                    or
                    l[5] not in ['0','A','C','G','T',]
                    ):
                    name_exclude = l[col_bim_name-1]
                ## both are rs (exclude current)
                elif l[col_bim_name-1][:2] == 'rs' and l_prev[col_bim_name-1][:2] == 'rs':
                    name_exclude = l[col_bim_name-1]
                ## current line is rs (exclude previous)
                elif l[col_bim_name-1][:2] == 'rs':
                    name_exclude = l_prev[col_bim_name-1]
                ## prev line is rs (exclude current)
                elif l_prev[col_bim_name-1][:2] == 'rs':
                    name_exclude = l[col_bim_name-1]
                ## none are rs (exclude current)
                else:
                    name_exclude = l_prev[col_bim_name-1]
                fd_out.write('%s\n' %(name_exclude))
            ## not a duplicate
            else:
                pass
            line_prev = line
        fd_in.close()
        fd_out.close()

##    cmd = 'cat %s.duplicates.SNPs %s.duplicates.SNPs | sort -u > duplicates.SNPs'

    return


def sort(l_bfiles,suffix,sort,):

    for bfile in l_bfiles:
        l_cmds = []
        l_cmds += ['cat %s.bim' %(bfile)]
        l_cmds += [sort]
##        l_cmds += ['''awk '{print $%i":"$%i,$%i}' ''' %(
##            col_bim_chr, col_bim_pos, col_bim_name,
##            )]
        cmd = ' | '.join(l_cmds)
        cmd += ' > %s.bim%s.sorted' %(bfile,suffix,)
        print cmd
        os.system(cmd)

    return


def bmerge(l_bfiles,l_populations,):

    '''
this function uses the files 1) quad2octo.dic and 2) extract.SNPs to
1) translate quad SNP names to octo SNP names and
2) extract SNPs common (by chromosome and position) to the quad and octo chip
'''

    ##
    ## make clean up dirs
    ##
    for extension in ['nof','hh','log',]:
        if not os.path.isdir(extension):
            os.mkdir(extension)
            pass

    for population in l_populations:
        bool_both_exist = True
        for bfile in l_bfiles:
            out = '%s_%s' %(bfile,population,)
            if not os.path.isfile('%s.bed' %(out)):
                bool_both_exist = False
                break
        if bool_both_exist == False:
            for bfile in l_bfiles:
                out = '%s_%s' %(bfile,population,)
                if os.path.isfile('%s.bed' %(out)):
                    for suffix in ['bed','fam','bim',]:
                        os.system('mv %s.%s %s.%s' %(
                            out,suffix,population,suffix,
                            ))
                    break
        else:

            l_out = ['%s_%s' %(bfile,population,) for bfile in l_bfiles]

            l_cmds = []

            ##
            ## rename quad rsIDs (--update-map --update-name)
            ##
            l_cmd = ['plink']
            l_cmd += ['--bfile %s' %(l_out[0],)]
            l_cmd += ['--make-bed --out %s_renamed' %(l_out[0],)]
            ## the file quad2octo.dic is generated by the function extract_and_translate
            l_cmd += ['--update-map quad2octo.dic --update-name']
            l_cmd += ['--noweb --nonfounders']
            cmd = ' \\\n'.join(l_cmd)
            l_cmds += [cmd]

            ##
            ## merge quad and octo (--bmerge)
            ##
            l_cmd = ['plink']
            l_cmd += ['--bmerge']
            l_cmd += ['%s_renamed.bed' %(l_out[0])]
            l_cmd += ['%s_renamed.bim' %(l_out[0])]
            l_cmd += ['%s_renamed.fam' %(l_out[0])]
            l_cmd += ['--bfile %s' %(l_out[1])]
            l_cmd += ['--make-bed --out %s' %(population)]
            ## the file extract.SNPs is generated by the function extract_and_translate
            l_cmd += ['--extract extract.SNPs']
            l_cmd += ['--noweb']
            cmd = ' \\\n'.join(l_cmd)
            l_cmds += [cmd]

            ##
            ## clean up
            ##
            for extension in ['nof','hh','log',]:
                l_cmds += ['mv %s.%s %s/.' %(population,extension,extension,)]
                l_cmds += ['mv %s_%s.%s %s/.' %(l_bfiles[0],population,extension,extension,)]
                l_cmds += ['mv %s_%s.%s %s/.' %(l_bfiles[1],population,extension,extension,)]
            for extension in ['bed','bim','fam','log','nof','hh',]:
                l_cmds += ['rm %s_renamed.%s' %(
                    l_out[0],extension,
                    )]
                for bfile in l_bfiles:
                    l_cmds += ['rm %s_%s.%s' %(bfile,population,extension,)]
            l_cmds += ['rm merge_%s.sh' %(population)]
            
            bsub = "bsub -M4000000 -R'select[mem>4000] rusage[mem=4000]' "
            bsub += '-P agv '
##            bsub += '-o out.txt '
##            bsub += '-e err.txt '
            bsub += '-q normal '
            bsub += "-J'%s' " %(population)

            fd = open('merge_%s.sh' %(population),'w')
            fd.write('\n\n'.join(l_cmds))
            fd.close()

            os.system('chmod +x merge_%s.sh' %(population))

            cmd_lsf = '%s ./merge_%s.sh' %(bsub, population)
            print cmd_lsf
            os.system(cmd_lsf)

            pass

        continue

    return


if __name__ == '__main__':
    main()
