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

    ## Declare pseudo matrix file name.
    ## Declare dictionary to be populated with
    ## IDs as keys and byte positions as values
    d_tell = {}
    d_cnt1 = {}
    d_cnt2 = {}
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
                for line in f:
                    l = line.rstrip().split()
                    ## Parse positions.
                    pos1 = int(l[2])
                    pos2 = int(l[3])
                    ## Skip if positions are distant from each other.
                    if abs(pos2-pos1) > d_args['window']:
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
                    ## Append count across populations.
                    for ID in (ID1,ID2):
                        try:
                            d_cnt1[ID] += 1
                            d_cnt2[ID] += 1
                        except KeyError:
                            d_cnt1[ID] = 1
                            d_cnt2[ID] = 1
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

    print('len(setS)', len(d_tell.keys()), len(d_cnt1.keys()), len(d_cnt2.keys()))
    setT = set()
    setQ = set()
    if d_args['preselected']:
        with open(d_args['preselected']) as f:
            for line in f:
                ID = line.rstrip()
                setT.add(ID)
                setQ.append(allinLD)
                setP.remove(ID)
                del d_cnt1[ID]
                for ID2 in d_LD[ID]:
                    d_LD[ID2].remove(ID)
                    d_cnt1[ID2] -= 1
                del d_LD[ID]
##                setT = set(f.read().rstrip().split('\n'))

    ## Declare dictionary to be populated with
    ## counts as keys and sets of IDs as values
    d_cnt2IDs = {}
    for ID,count in d_cnt1.items():
        try:
            d_cnt2IDs[count].add(ID)
        except KeyError:
            d_cnt2IDs[count] = set([ID])
##    del cnt

    d_args['max_tagSNP'] = min(d_args['max_tagSNP'],len(d_tell.keys()))

    select_tag_SNPs(d_args, d_cnt2IDs, setT, setQ, d_cnt1, d_cnt2, d_tell)

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
    print(
        'len(setS)', len(d_tell.keys()),
        len(d_cnt1.keys()), len(d_cnt2.keys()))

    return


def select_tag_SNPs(d_args, d_cnt2IDs, setT, setQ, d_cnt1, d_cnt2, d_tell):
    
    most_common_key = max(d_cnt2IDs.keys())
    with open(d_args['out']+'.LDmatrix') as f:
        while len(setT) < d_args['max_tagSNP']:
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
                (d_cnt2[ID]-d_cnt1[ID])/d_cnt2[ID]
                > d_args['max_coverage_regional']):
                del d_cnt1[ID]
                del d_tell[ID]
                continue
            else:
                setT |= set([ID])
                del d_cnt1[ID]
                pass
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
                        count = d_cnt1[ID2]
                    ## max_coverage_regional
                    except KeyError:
                        continue
                    d_cnt2IDs[count].remove(ID2)
                    try:
                        d_cnt2IDs[count-1].add(ID2)
                    except KeyError:
                        d_cnt2IDs[count-1] = set([ID2])
                    d_cnt1[ID2] -= 1
                    if d_cnt1[ID2] < 0: print(ID,ID2,d_cnt1[ID2]); stop ## tmp!!!
                    setQ.add(ID2) ## tmp!!!
            del d_tell[ID]
            sys.stdout.flush()

    return setT, setQ


def argparser():

    parser = argparse.ArgumentParser()

    parser.add_argument('--in', required=True, nargs='+')
    parser.add_argument('--out', required=True)
    parser.add_argument('--min_LD', type=float, default=.8)
    parser.add_argument('--min_MAF', type=float, default=.01)
    parser.add_argument('--window', type=int, default=250000)
    parser.add_argument('--preselected')
    parser.add_argument('--max_tagSNP', type=int, default=1000000)
    parser.add_argument('--max_coverage_regional', type=float, default=.2)

    d_args = vars(parser.parse_args())

    return d_args


if __name__ == '__main__':
    main()
