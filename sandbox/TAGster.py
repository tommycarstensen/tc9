#!python3

## Tommy Carstensen, Wellcome Trust Sanger Institute, April 2014

import argparse
import os
import sys
import collections
import gzip


def main():

    ## Parse arguments.
    d_args = argparser()

    if d_args['pretagged']:
        set_tmp_ignore = set()
        with open(d_args['pretagged']) as f:
            set_tmp_ignore = set(line.rstrip() for line in f)
    else:
        set_tmp_ignore = set()

    if d_args['preselected']:
        with open(d_args['preselected']) as f_preselected:
            set_preselected = set(line.rstrip() for line in f_preselected)
    else:
        set_preselected = set()

    ## tag SNPs should not be ignored.
    set_tmp_ignore -= set_preselected

    ## Declare pseudo matrix file name.
    ## Declare dictionary to be populated with
    ## IDs as keys and byte positions as values
    d_tell = {}
    d_cnt = {}
    ## Open pseudo matrix file.
    with open(d_args['out']+'.LDmatrix', 'w') as f_matrix:
        ## Loop over input files.
        for file in d_args['in']:
            print(file)
            sys.stdout.flush()
            ## Declare temporary LD dictionary.
            d_LD = {}
            ## Open input file.
            with gzip.open(file,'rt') as f:
                ## Loop over input file lines.
                ## TODO: MAKE A GENERATOR TO SHORTEN!!! AFTER DEADLINE!!!
                for line in f:
                    l = line.rstrip().split()
                    ## Parse positions.
                    pos1 = int(l[2])
                    pos2 = int(l[3])
                    ## Skip if positions are distant from each other.
                    if abs(pos2-pos1) > d_args['max_window']:
                        continue
                    ## Parse LD.
                    LD = float(l[4])
                    ## Skip if LD below threshold.
                    if LD < d_args['min_LD']:
                        continue
                    ## Parse MAF.
                    MAF1 = float(l[5])
                    MAF2 = float(l[6])
                    ## Skip if either SNP below MAF threshold.
                    if MAF1 < d_args['min_MAF'] or MAF2 < d_args['min_MAF']:
                        continue
                    ## Parse SNP IDs.
                    ID1 = l[0]
                    ID2 = l[1]
                    ## Skip if tagged by IMPUTE2.
                    ## CLEAN UP THIS UGLY CODE AFTER THE DEADLINE!!!
                    if ID1 in set_tmp_ignore or ID2 in set_tmp_ignore:
                        if ID1 in set_tmp_ignore and ID2 in set_tmp_ignore:
                            continue
                        elif ID1 in set_tmp_ignore:
                            ID = ID2
                        else:
                            ID = ID1
                        try:
                            d_cnt[ID][0] += 0
                            d_cnt[ID][1] += 0
                        except KeyError:
                            d_cnt[ID] = [0,0]
                            d_LD[ID] = []
                        continue
                    ## Append count across populations.
                    for ID in (ID1,ID2):
                        try:
                            d_cnt[ID][0] += 1
                            d_cnt[ID][1] += 1
                        except KeyError:
                            d_cnt[ID] = [1,1]
                    ## Append SNPs in LD to temporary LD dictionary.
                    for IDa,IDb in ((ID1,ID2),(ID2,ID1)):
                        try:
                            d_LD[IDa].append(IDb)
                        except KeyError:
                            d_LD[IDa] = [IDb]
                    ## Continue loop over input file lines.
                    continue
                ## Close input file.
                pass
            ## Append IDs in temporary LD dictionary to pseudo matrix file.
            for ID,l_IDs in d_LD.items():
                ## Get pseudo matrix file byte position.
                tell = f_matrix.tell()
                ## Append IDs to pseudo matrix file.
                f_matrix.write('%s\n' %(' '.join(l_IDs)))
                ## Append byte position to dictionary.
                try:
                    d_tell[ID].append(tell)
                except KeyError:
                    d_tell[ID] = [tell]

            ## Continue loop over input files.
            continue
##            break ## tmp!!!
        ## Delete temporary LD dictionary from memory.
        del d_LD
        ## Close pseudo matrix file.
        pass

    ## Declare dictionary to be populated with
    ## counts as keys and sets of IDs as values
    d_cnt2IDs = {}
    for ID,t_count in d_cnt.items():
        count = t_count[0]
        try:
            d_cnt2IDs[count].add(ID)
        except KeyError:
            d_cnt2IDs[count] = set([ID])
##    del cnt

    ##
    ## Update counts based on preselected tag SNPs.
    ##
    print('len(setS)', len(d_tell.keys()), len(d_cnt.keys()))
    setT = set()
    setQ = set()
    with open(d_args['out']+'.LDmatrix') as f_matrix:
        for ID in set_preselected:
            setT.add(ID)
            assert ID not in set_tmp_ignore
##            print(
##                ID, '20:18463897' in d_cnt.keys(),
##                ID in set_tmp_ignore, ID in set_preselected,
##                ID in d_tell.keys(), ID in d_cnt.keys())
            setQ = update_counts_and_set_of_tagged_SNPs(
                ID, d_tell, f_matrix, setT, d_cnt, d_cnt2IDs, setQ)
            d_cnt2IDs[d_cnt[ID][0]].remove(ID)
            del d_cnt[ID]

    del set_tmp_ignore

    d_args['max_tagSNP'] = min(d_args['max_tagSNP'],len(d_tell.keys()))

    print('selecting tag SNPs')
    sys.stdout.flush()
    select_tag_SNPs(d_args, d_cnt2IDs, setT, setQ, d_cnt, d_tell)

    with open(d_args['out']+'.tagSNPs','w') as f:
        for ID in setT:
            f.write('%s\n' %(ID))
    with open(d_args['out']+'.tagged','w') as f:
        for ID in setQ:
            f.write('%s\n' %(ID))

    ## Clean up large temporary files.
    os.remove(d_args['out']+'.LDmatrix')

    print('len(setT)',len(setT))
    print('len(setQ)',len(setQ))
    print('len(setQ-setT)',len(setQ-setT))
    print('len(setQ+setT)',len(setQ.union(setT)))
    print('len(setS)', len(d_tell.keys()), len(d_cnt.keys()))

    return


def select_tag_SNPs(d_args, d_cnt2IDs, setT, setQ, d_cnt, d_tell):

    len_setS = len(set(d_cnt.keys()))
    
    most_common_key = max(d_cnt2IDs.keys())
    with open(d_args['out']+'.LDmatrix') as f:
        while len(setT) < d_args['max_tagSNP']:
            if len(setQ) >= d_args['max_coverage_global']*len_setS:
                break
            ## Find the most common element.
            ## collections.Counter.most_common(1) is too slow.
            while True:
                try:
                    ID = d_cnt2IDs[most_common_key].pop()
                    break
                except KeyError:
                    del d_cnt2IDs[most_common_key]
                    try:
                        most_common_key = max(d_cnt2IDs.keys())
                    except ValueError:
                        ## Dictionary exhausted.
                        if not d_cnt2IDs:
                            return setT, setQ
                        stopshouldnothappen
                        break
            if (
                (d_cnt[ID][0]-d_cnt[ID][1])/d_cnt[ID][1]
                > d_args['max_coverage_regional']):
                del d_cnt[ID]
                del d_tell[ID]
                continue
            else:
                setT.add(ID)
                del d_cnt[ID]
                pass
            setQ = update_counts_and_set_of_tagged_SNPs(
                ID, d_tell, f, setT, d_cnt, d_cnt2IDs, setQ)

    return setT, setQ


def update_counts_and_set_of_tagged_SNPs(
    ID, d_tell, f, setT, d_cnt, d_cnt2IDs, setQ):

    ## File seek is the second slowest step in this loop.
    for seek in d_tell[ID]:
        f.seek(seek)
        line = f.readline()
        for ID2 in line.rstrip().split():
            ## tag SNP
            if ID2 in setT:
                continue
            ## tagged
            try:
                count = d_cnt[ID2][0]
            ## max_coverage_regional
            except KeyError:
                continue
            d_cnt2IDs[count].remove(ID2)
            try:
                d_cnt2IDs[count-1].add(ID2)
            except KeyError:
                d_cnt2IDs[count-1] = set([ID2])
            d_cnt[ID2][0] -= 1
            if d_cnt[ID2][0] < 0: print(ID,ID2,d_cnt[ID2]); stop ## tmp!!!
            setQ.add(ID2) ## tmp!!!
    del d_tell[ID]

    return setQ


def argparser():

    parser = argparse.ArgumentParser()

    parser.add_argument('--in', required=True, nargs='+')
    parser.add_argument('--out', required=True)
    parser.add_argument('--min_LD', type=float, default=.8)
    parser.add_argument('--min_MAF', type=float, default=.01)
    parser.add_argument('--max_window', type=int, default=250000)
    parser.add_argument('--preselected')
    parser.add_argument('--pretagged','--ignore')
    parser.add_argument('--max_tagSNP', type=int, default=1000000)
    parser.add_argument('--max_coverage_regional', type=float, default=1)
    parser.add_argument('--max_coverage_global', type=float, default=1)

    d_args = vars(parser.parse_args())

    if d_args['preselected'] == 'None':
        d_args['preselected'] = None
    if d_args['pretagged'] == 'None':
        d_args['pretagged'] = None

    return d_args


if __name__ == '__main__':
    main()
