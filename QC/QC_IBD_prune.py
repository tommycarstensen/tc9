#!/software/bin/python

## Tommy Carstensen, Wellcome Trust Sanger Institute, 2012 and November 2013

import os
import sys
import gzip

def main():

    ## if two samples are related, then only one of them is removed
    ## but which one is removed is random
    ## i.e. the first occurence among sorted samples is removed

    d_ID2IDs = parse_genome()

    d_count2IDs = count_relations(d_ID2IDs)

    if '--imiss' in sys.argv:
        d_imiss = parse_imiss()
    else:
        d_imiss = {ID: 0 for ID in d_ID2IDs.keys()}

    l_twins = parse_twins()

    l_exclude = generate_exclusion(
        d_ID2IDs,d_count2IDs,d_imiss=d_imiss,l_twins=l_twins)

    write_exclusion_list(l_exclude)

    return


def parse_twins():

    l_twins = []

    if not '--twins' in sys.argv:
        return l_twins

    twins = sys.argv[sys.argv.index('--twins')+1]
    f = open(twins,'r')
    l_twins = [line.split()[0] for line in f.readlines()]
    f.close()

    return l_twins


def parse_imiss():

    imiss = sys.argv[sys.argv.index('--imiss')+1]

    d_imiss = {}

    fd = open('%s.imiss' %(imiss),'r')
    ## skip header
    for line in fd: break
    ## loop lines
    for line in fd:
        l = line.split()
        ID = l[0]+l[1]
        F_MISS = float(l[5])
        d_imiss[ID] = F_MISS
    fd.close()

    return d_imiss


def write_exclusion_list(l_exclude,):

    pi_hat_max = float(sys.argv[sys.argv.index('--pi_hat_max')+1])
    out = sys.argv[sys.argv.index('--out')+1]

##    s = '\n'.join(['%s %s' %(IDa,IDa,) for IDa in l_exclude])+'\n'
    l = ['%s %s\n' %(IDa,IDa,) for IDa in l_exclude]

    fd = open('%s.genome.%.2f.samples' %(out,pi_hat_max,),'w')
    fd.writelines(l)
    fd.close()

    return


def generate_exclusion(
    d_ID2IDs,d_count2IDs,d_imiss=None,verbose=False,l_twins=[],
    ):

    l_exclude = []
    while d_count2IDs != {}:
        count_max = max(d_count2IDs.keys())
        if count_max == 0:
            stop
            break
        l_IDas = list(d_count2IDs[count_max])
        ## sort so ID with largest F_MISS is deleted first
        if d_imiss:
            l = [(d_imiss[IDa],IDa,) for IDa in l_IDas]
            l.sort()
            l.reverse()
            l_IDas = [t[1] for t in l]
        else:
            ## sort so first occurence(s) are deleted first
            l_IDas.sort()
        for IDa in l_IDas:
            bool_exclude = False
            l_IDbs = list(d_ID2IDs[IDa])
            if verbose == True:
                print(count_max, IDa, l_IDbs)
                pass
            for IDb in l_IDbs:
                if verbose == True:
                    print(count_max, IDa,IDb)
                    pass
                ## only exclude if not twin
                if not (IDa in l_twins and IDb in l_twins):
                    bool_exclude = True
                count = len(d_ID2IDs[IDb])
                d_ID2IDs[IDa].remove(IDb)
                d_ID2IDs[IDb].remove(IDa)
                d_count2IDs[count].remove(IDb)
                if count > 1:
                    if not count-1 in d_count2IDs.keys():
                        d_count2IDs[count-1] = []
                        pass
                    d_count2IDs[count-1] += [IDb]
                    pass
                if d_count2IDs[count] == []:
                    del d_count2IDs[count]
                    pass
                continue
            d_count2IDs[count_max].remove(IDa)
            del d_ID2IDs[IDa]
            ## only exclude if not twin
            if bool_exclude:
                l_exclude += [IDa]
            break
        if d_count2IDs[count_max] == []:
            del d_count2IDs[count_max]
            pass
        continue

    l_exclude.sort()
    print('%i samples to be removed' %(len(l_exclude)))

    return l_exclude


def count_relations(d_ID2IDs):

    d_count2IDs = {}
    for ID,l_IDs in d_ID2IDs.items():
        count = len(l_IDs)
        if not count in d_count2IDs.keys():
            d_count2IDs[count] = []
            pass
        d_count2IDs[count] += [ID]

    return d_count2IDs


def parse_genome():

    genome = sys.argv[sys.argv.index('--genome')+1]
    pi_hat_max = float(sys.argv[sys.argv.index('--pi_hat_max')+1])

    if os.path.isfile('{}.genome'.format(genome)):
        fd = open('{}.genome'.format(genome),'r')
    elif os.path.splitext(genome)[1] == '.gz':
        fd = gzip.open(genome, 'rt')
    else:
        fd = open(genome)
    d_ID2IDs = {}
    ## skip header
    for line in fd: break
    ## loop lines
    for line in fd:
        l = line.split()
        ID1 = l[0]+l[1]
        ID2 = l[2]+l[3]
        PI_HAT = float(l[9])
##        if PI_HAT < pi_hat_max:
        if PI_HAT <= pi_hat_max:
            continue
        for IDa,IDb in [
            [ID1,ID2,],
            [ID2,ID1,],
            ]:
            if not IDa in d_ID2IDs.keys():
                d_ID2IDs[IDa] = []
                pass
            d_ID2IDs[IDa] += [IDb]
            continue
        continue
    fd.close()

    return d_ID2IDs


def unit_test():

    d_ID2IDs = {
        'a':['b','c','d','e','f',],
        'b':['a','d','e','g','h',],
        'c':['a','d','i',],
        'd':['a','b','c','h',],
        'e':['a','b','j',],
        'f':['a',],
        'g':['b',],
        'h':['b','d',],
        'i':['c',],
        'j':['e',],
        'k':['l',],
        'l':['k',],
        'm':['n','o',],
        'n':['m','o',],
        'o':['m','n',],
        }
    d_count2IDs = count_relations(d_ID2IDs)
    l_exclude = generate_exclusion(d_ID2IDs,d_count2IDs,verbose=True,)
    print(l_exclude == ['a', 'b', 'c', 'd', 'e', 'k', 'm', 'n'])

    d_ID2IDs = {
        'a':['b',],
        'b':['a',],
        }
    d_count2IDs = count_relations(d_ID2IDs)
    l_exclude = generate_exclusion(d_ID2IDs,d_count2IDs,verbose=True,)
    print(l_exclude == ['a',])

    d_ID2IDs = {
        'a':['b','c',],
        'b':['a','c',],
        'c':['a','b',],
        }
    d_count2IDs = count_relations(d_ID2IDs)
    l_exclude = generate_exclusion(d_ID2IDs,d_count2IDs,verbose=True,)
    print(l_exclude == ['a','b',])

    return


if __name__ == '__main__':
    main()
##    unit_test()
