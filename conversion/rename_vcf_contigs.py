## Tommy Carstensen, Wellcome Trust Sanger Institute, 2012

import os

def main():

    l_fn = [
        'in_GATK/hapmap_3.3.b37.sites.vcf',
        'in_GATK/1000G_omni2.5.b37.sites.vcf',
        'in_GATK/dbsnp_132.b37.vcf'
        ]

    for fn in l_fn:

        print 'renaming in', fn

##        ##
##        ## 1) read all lines into memory (and run out of memory)
##        ##
##        fd = open(fn,'r')
##        lines = fd.readlines()
##        fd.close()
##
##        for i_line in xrange(len(lines)):
##            line = lines[i_line]
##            if line[0] == '#':
##                continue
##            if line[:3] == 'Chr':
##                lines[i_line] = line[3:]
##                continue
##
##        fd = open(fn,'w')
##        fd.writelines(lines)
##        fd.close()

        ##
        ## 2) use sed
        ## much slower but requires less memory, as it operates line-by-line
        ##
        ## replace Chr with empty string if line does not start with #
        ## -i = --in-place
        ## also check http://codept.blogspot.co.uk/2007/12/sed-remove-first-4-letters-in-each-line.html
        os.system("sed -i '/^#/! s/Chr//g' %s" %(fn))

    return


if __name__ == '__main__':
    main()
