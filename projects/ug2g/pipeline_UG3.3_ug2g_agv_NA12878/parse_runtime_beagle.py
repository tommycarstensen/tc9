import glob
import re
import datetime
import time
import os

format_strptime = '%b %d %H:%M:%S %Y'
for file in glob.glob('LSF/beagle/*/*.out'):
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
                if delta.seconds > 60:
                    print(delta.days*24+int(delta.seconds/(60*60)))
                    pass
                if delta.days*24*60*60+delta.seconds > delta_max:
                    delta_max = delta.days*24*60*60+delta.seconds
    if delta_max < 3600 and not line.strip():
        if os.path.isfile('out_{}.vcf.gz.tbi'.format(file[4:-4])):
            os.remove(file)
        elif not os.path.isfile('out_{}.vcf.gz'.format(file[4:-4])):
            continue
        else:
            print(file)
            print('out_{}.log'.format(file[4:-4]))
            print('out_{}.vcf.gz'.format(file[4:-4]))
            print('vcf', (time.time()-os.path.getmtime(
                'out_{}.vcf.gz'.format(file[4:-4])))/3600)
            print('log', (time.time()-os.path.getmtime(
                'out_{}.log'.format(file[4:-4])))/3600)
            print(os.path.getsize('out_{}.vcf.gz'.format(file[4:-4])))
            print(os.path.isfile('out_{}.log'.format(file[4:-4])))
        if not any([
            os.path.isfile('out_{}.log'.format(file[4:-4])),
            os.path.isfile('out_{}.vcf.gz'.format(file[4:-4])),
            ]):
            print('removing', file)
            os.remove(file)
