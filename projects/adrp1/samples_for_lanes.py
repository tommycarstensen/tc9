import xlrd
import random
import math
from collections import Counter

path = 'RNA_samples_randomisation_plan_and_plate_layout_050116.xlsx'
xl_workbook = xlrd.open_workbook(path)

xl_sheet = xl_workbook.sheet_by_name('sheet 1')

print('#', xl_sheet.cell(0, 0))
print('#', xl_sheet.cell(0, 1))
print('#', xl_sheet.cell(0, 2))
print('#', xl_sheet.cell(0, 3))

pops = ('MKK', 'LWK', 'YRI', 'ESN', 'GWD', 'MSL',)

n_samples = 579 + 21  # 579 + 21 YRI replicates
n_hb = n_Hologic_batches = math.ceil(n_samples / 48)
n_cb = n_Coriell_batches = math.ceil(n_samples / 100)
n_wells = 96
n_plates = n_Sanger_plates = math.ceil(n_samples / n_wells)
n_average_samples_per_lane = 6  # 12 samples per 2 lanes
n_samples_per_lane = 12
print('#n_Hologic_batches', n_Hologic_batches)
print('#n_Coriell_batches', n_Coriell_batches)
print('#n_Sanger_plates', n_Sanger_plates)
print()

d = {hb: {cb: [] for cb in range(1, n_cb+1)} for hb in range(1, n_hb+1)}
d_ID2pop = {}
d_plate = {
    plate: {
        'init': Counter(), 'curr': Counter()} for plate in range(
            1, n_plates+1)}

for row_idx in range(1, xl_sheet.nrows):
    ID = Coriell_ID = xl_sheet.cell(row_idx, 0).value
    pop = Population = xl_sheet.cell(row_idx, 1).value
    cb = Coriell_Batch = xl_sheet.cell(row_idx, 2).value
    hb = Hologic_Batch = xl_sheet.cell(row_idx, 3).value
    well = xl_sheet.cell(row_idx, 5).value
    d[hb][cb].append(ID)
    d[hb][cb].append(ID)
    d_ID2pop[ID] = pop
    plate = ((hb-1)//2)+1
    if plate == 6 and int(well[1:]) >= 10:
        plate = 7
    d_plate[plate]['init'][pop] += 2
    d_plate[plate]['curr'][pop] += 2

for plate in range(1, n_plates+1):
    n_samples = int(sum(d_plate[plate]['init'].values()) / 2)
    n_lanes = int(n_samples / n_average_samples_per_lane)
    for lane in range(1, n_lanes+1):
        selected_pops = []
        for hb in range(plate*2-1, min(plate*2, n_hb)+1):
            for cb in range(1, n_cb+1):
                ## Count the populations in the remaining set of samples.
                c = Counter([d_ID2pop[ID] for ID in d[hb][cb]])
                ## (Re)set counter.
                i = 0
                ## Calculate the current maximal proportion of each population on the plate.
                frac_max = max([d_plate[plate]['curr'][_]/d_plate[plate]['init'][_] for _ in pops if d_plate[plate]['init'][_] != 0])
                while True:
                    ## Randomly choose an ID.
                    ID = random.choice(d[hb][cb])
                    ## Get the pop the ID belongs to.
                    pop = d_ID2pop[ID]
                    ## Calculate the proportion of the population on the plate
                    ## to which the randomly selected sample belongs to.
                    frac = d_plate[plate]['curr'][pop]/d_plate[plate]['init'][pop]
                    if frac == frac_max:
                        break
                    ## Break if we tried a thousand random choices
                    ## and this sample belongs to the most abundant
                    ## of the remaining populations
                    ## for the particular combination of hb and cb.
                    if pop == max(c) and i > 1000:
                        break
                    i += 1
                d_plate[plate]['curr'][pop] -= 1
                d[hb][cb].remove(ID)
                ## Print lane for each ID.
                print(ID, d_ID2pop[ID], cb, hb, lane, sep='\t')
                selected_pops.append(d_ID2pop[ID])
##        ## Print a summary for the lane.
##        print(
##            lane+(plate-1)*int(n_wells/n_average_samples_per_lane),
##            selected_pops.count('ESN'), selected_pops.count('GWD'),
##            selected_pops.count('LWK'), selected_pops.count('MKK'),
##            selected_pops.count('MSL'), selected_pops.count('YRI'),
##            sep='\t',
##            )
