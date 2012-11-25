## Tommy Carstensen, Wellcome Trust Sanger Institute, 2012

import os

def main():

    fd = open('omni_keep_flip.map','r')
    lines = fd.readlines()
    fd.close()
    for line in lines:
        print line
        stop

    fd = open('omni_keep_flip.ped','r')
    for line in fd:
        print line
        stop

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
