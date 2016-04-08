import glob
import re
import datetime
import time
import os

format_strptime = '%b %d %H:%M:%S %Y'
print('nthreads CPUhours MaxMemGB CPUhours*MaxMemGB file')
for file in glob.glob('test/LSF/*.out'):
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
            if 'Max Memory' in line:
                MaxMemoryGB = float(line.split()[3])/1024
                CPUhours = delta.days*24+delta.seconds/(60*60)
                if delta.days > 0 or delta.seconds > 180:
                    print(nthreads, round(CPUhours,2), round(MaxMemoryGB,2), round(CPUhours*MaxMemoryGB,2), file, sep='\t')
            if match:
                t2 = datetime.datetime.strptime(match.group(1), format_strptime)
                delta = t2-t1
                if delta.days == 0 and delta.seconds < 180:
                    continue
#                    print(delta.days*24+int(delta.seconds/(60*60)))
                    pass
                nthreads =os.path.basename(file).split('.')[0]
#                print(nthreads, delta.days*24*60*60+delta.seconds, file, sep='\t')
                if delta.days*24*60*60+delta.seconds > delta_max:
                    delta_max = delta.days*24*60*60+delta.seconds
