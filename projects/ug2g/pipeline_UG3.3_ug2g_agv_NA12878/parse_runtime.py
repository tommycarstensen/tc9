import glob
import re
import datetime
import time
import os

format_strptime = '%b %d %H:%M:%S %Y'
sum_delta_days = 0
for step in 'VariantRecalibrator ApplyRecalibration beagle UnifiedGenotyper'.split():
    for file in glob.glob('LSF/{}/*/*.out'.format(step)):
        t1 = t2 = 0
        delta_max = 0
        with open(file) as f:
            for line in f:
                pattern = 'Started at [A-Z]\w\w ([A-Z]\w\w .\d \d\d:\d\d:\d\d 201\d)'
                match = re.match(pattern, line)
                if match:
                    t1 = datetime.datetime.strptime(match.group(1), format_strptime)
                pattern = 'Results reported at [A-Z]\w\w ([A-Z]\w\w .\d \d\d:\d\d:\d\d 201\d)'
                match = re.match(pattern, line)
                if match:
                    t2 = datetime.datetime.strptime(match.group(1), format_strptime)
                    delta = t2-t1
                    delta_days = delta.days+delta.seconds/(24*60*60)
                    sum_delta_days += delta_days
    print(step, sum_delta_days)
