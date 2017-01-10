import xlrd
import random
import math
##from collections import Counter

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
print('#n_Hologic_batches', n_Hologic_batches)
print('#n_Coriell_batches', n_Coriell_batches)
print('#n_Sanger_plates', n_Sanger_plates)
print()

d = {hb: {cb: [] for cb in range(1, n_cb+1)} for hb in range(1, n_hb+1)}
d_ID2pop = {}

for row_idx in range(1, xl_sheet.nrows):
    ID = Coriell_ID = xl_sheet.cell(row_idx, 0).value
    pop = Population = xl_sheet.cell(row_idx, 1).value
    cb = Coriell_Batch = xl_sheet.cell(row_idx, 2).value
    hb = Hologic_Batch = xl_sheet.cell(row_idx, 3).value
    d[hb][cb].append(ID)
    d[hb][cb].append(ID)
    d_ID2pop[ID] = pop

for plate in range(1, n_plates+1):
    if plate == n_plates:
        n_lanes = int((n_samples % n_wells) / n_average_samples_per_lane)
    else:
        n_lanes = int(n_wells / n_average_samples_per_lane)
    for lane in range(1, n_lanes+1):
        for hb in range(plate*2-1, min(plate*2, n_hb)+1):
            for cb in range(1, n_cb+1):
                ID = random.choice(d[hb][cb])
                d[hb][cb].remove(ID)
                print(ID, d_ID2pop[ID], cb, hb, lane, sep='\t')
