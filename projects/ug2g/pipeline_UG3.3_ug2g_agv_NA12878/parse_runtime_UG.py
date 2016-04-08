import glob
import re
import datetime
import time

format_strptime = '%b %d %H:%M:%S %Y'
for file in glob.glob('LSF/UnifiedGenotyper/*/*.out'):
    t1 = t2 = nt = nct = 0
    with open(file) as f:
        for line in f:
            pattern = '.*--nt (\d{1,2})'
            match = re.match(pattern, line)
            if match:
                nt = match.group(1)
            pattern = '.*--nct (\d{1,2})'
            match = re.match(pattern, line)
            if match:
                nct = match.group(1)
                print('{}\t{}\t{}'.format(nt, nct, delta.seconds))
                t1 = t2 = nt = nct = 0
            pattern = 'Started at [A-Z]\w\w ([A-Z]\w\w .\d \d\d:\d\d:\d\d 201\d)'
            match = re.match(pattern, line)
            if match:
                t1 = datetime.datetime.strptime(match.group(1), format_strptime)
            pattern = 'Results reported at [A-Z]\w\w ([A-Z]\w\w .\d \d\d:\d\d:\d\d 201\d)'
            match = re.match(pattern, line)
            if match:
                t2 = datetime.datetime.strptime(match.group(1), format_strptime)
                delta = t2-t1
