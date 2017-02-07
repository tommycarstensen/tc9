## Adapted from https://gist.github.com/atdt/875e0dba6a15e3fa6018

# Python implementation of Aho-Corasick string matching
#
# Alfred V. Aho and Margaret J. Corasick, "Efficient string matching: an aid to
# bibliographic search", CACM, 18(6):333-340, June 1975.
#
# <http://xlinux.nist.gov/dads//HTML/ahoCorasick.html>
#
# Copyright (C) 2015 Ori Livneh <ori@wikimedia.org>
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
FAIL = -1

import argparse
import itertools

def main():

    args = parse_args()

    print(args.needles)
    with open(args.haystack) as haystack:
        results = aho_corasick(haystack, args.needles)
    print(results)

    return


def parse_args():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        '--string','--haystack', '--ref',
        dest='haystack',
        help='file/string to be searched through',
        required = True,
        )

    parser.add_argument(
        '--substrings', '--needles',
        dest='needles',
        help='substring(s) to search for',
        required = True,
        nargs='+',
        )

    args = parser.parse_args()

    return args


def aho_corasick(haystack, keywords):

    transitions = {}
    outputs = {}
    fails = {}

    new_state = 0

    for keyword in keywords:
        state = 0

        for j, char in enumerate(keyword):
            res = transitions.get((state, char), FAIL)
            if res == FAIL:
                break
            state = res

        for char in keyword[j:]:
            new_state += 1
            transitions[(state, char)] = new_state
            state = new_state

        outputs[state] = [keyword]

    queue = []
    for (from_state, char), to_state in transitions.items():
        if from_state == 0 and to_state != 0:
            queue.append(to_state)
            fails[to_state] = 0

    while queue:
        r = queue.pop(0)
        for (from_state, char), to_state in transitions.items():
            if from_state == r:
                queue.append(to_state)
                state = fails[from_state]

                while True:
                    res = transitions.get((state, char), state and FAIL)
                    if res != FAIL:
                        break
                    state = fails[state]

                failure = transitions.get((state, char), state and FAIL)
                fails[to_state] = failure
                outputs.setdefault(to_state, []).extend(
                    outputs.get(failure, []))

    state = 0
    results = []
    i = 0
    CHROM = 0
    ## This should be preceded by a loop over chromosomes
    ## to prevent the unlikely event of string matching across chromosomes.
    for char in itertools.chain.from_iterable(haystack):
        ## next chromosome
        if char == '>':
             CHROM += 1
             POS = 0
             continue
        if char == '\n':
             continue
        i += 1
        ## Position on chromosome.
        POS += 1
#        if not POS % 10000:
#            print(CHROM, POS)
        while True:
            res = transitions.get((state, char), state and FAIL)
            if res != FAIL:
                state = res
                break
            state = fails[state]

        for match in outputs.get(state, ()):
#            pos = i - len(match) + 1
            print(CHROM, POS - len(match) + 1, match)
            results.append((POS - len(match) + 1, match))

    return results


if __name__ == '__main__':
    main()
