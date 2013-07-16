#!/software/python-3.3.2/bin/python3

import argparse

def main():

    needle, haystack = func_argparse()

    s = ''
    with open(haystack) as text:
        for chrom,pos in KnuthMorrisPratt(text, needle):
            print(chrom,pos+1)
            s += '%i %i\n' %(chrom,pos+1)
    with open('kmp.out','w') as file:
        file.write(s)

    return


def build_shifts(pattern):

    # build table of shift amounts
    shifts = [1] * (len(pattern) + 1)
    shift = 1
    for pos in range(len(pattern)):
        while shift <= pos and pattern[pos] != pattern[pos-shift]:
            shift += shifts[pos-shift]
        shifts[pos+1] = shift

    return shifts


def func_argparse():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        '--string','--haystack',
        dest='haystack',
        help='file/string to be searched through',
        required = True,
        )

    parser.add_argument(
        '--substring','--needle',
        dest='needle',
        help='substring to search for',
        required = True,
        )

    haystack = vars(parser.parse_args())['haystack']
    needle = vars(parser.parse_args())['needle']

    return needle, haystack


def KnuthMorrisPratt(text, pattern):

    # This function is heavily based on:
    # Knuth-Morris-Pratt string matching
    # David Eppstein, UC Irvine, 1 Mar 2002

    '''Yields all starting positions of copies of the pattern in the text.
Calling conventions are similar to string.find, but its arguments can be
lists or iterators, not just strings, it returns all matches, not just
the first one, and it does not need the whole text in memory at once.
Whenever it yields, it will have read the text exactly up to and including
the match that caused the yield.'''

    # allow indexing into pattern and protect against change during yield
    pattern = list(pattern)

    # build table of shift amounts
    shifts = [1] * (len(pattern) + 1)
    shift = 1
    for pos in range(len(pattern)):
        while shift <= pos and pattern[pos] != pattern[pos-shift]:
            shift += shifts[pos-shift]
        shifts[pos+1] = shift

    # do the actual search
    chrom = 0
    while True:
        c = text.read(1)
        if not c:
            break
        ## next chromosome
        if c == '>':
             text.readline()
             chrom += 1
             startPos = 0
             matchLen = 0
             continue
        if c == '\n':
             continue
        while matchLen == len(pattern) or \
              matchLen >= 0 and pattern[matchLen] != c:
            startPos += shifts[matchLen]
            matchLen -= shifts[matchLen]
        matchLen += 1
        if matchLen == len(pattern):
            yield chrom, startPos

if __name__ == '__main__':
    main()
