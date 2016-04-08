import glob
import re
import datetime
import time

format_strptime = '%b %d %H:%M:%S %Y'
for file in glob.glob('LSFtest/SHAPEIT.20.*out'):
    with open(file) as f:
        for line in f:
            pattern = 'Started at [A-Z]\w\w ([A-Z]\w\w \d{1,2} \d\d:\d\d:\d\d 201\d)'
            match = re.match(pattern, line)
            if match:
                break
        t1 = datetime.datetime.strptime(match.group(1), format_strptime)
        line = f.readline()
        pattern = 'Results reported at [A-Z]\w\w ([A-Z]\w\w \d{1,2} \d\d:\d\d:\d\d 201\d)'
        match = re.match(pattern, line)
        t2 = datetime.datetime.strptime(match.group(1), format_strptime)
        delta = t2-t1
        pattern = 'LSFtest/SHAPEIT.20.(\d{1,2}).\d{1,2}.out'
        match = re.match(pattern, file)
        threads = int(match.group(1))
        print('{}\t{}'.format(threads, delta.seconds))
