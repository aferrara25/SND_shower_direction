import os, sys
from glob import glob

min = -1
max = -1
tb = "testbeam"
ti18 = "physics"



if len(sys.argv) < 2 :
    print(f'Usage: python write_args.py data_path min_run max_run append?. \nUse -1 for min and max to write all folder info.')
    sys.exit(2)

if len(sys.argv) > 3 :
    min = int(sys.argv[2])
    max = int(sys.argv[3])

data_path = sys.argv[1]

if tb in data_path : bool = 'true'
elif ti18 in data_path : bool = 'false'

runlist = glob(f"{data_path}/run*")

#print(runlist, flush = True)

if not bool:    runnumbers = [int((run.split("_")[len(run.split("_"))-1]).lstrip('0') )for run in runlist]
else :          runnumbers = [int(run.split("_")[len(run.split("_"))-1]) for run in runlist]

runfiles = [len( glob(f"{data_path}/run_*{str(run)}/sndsw_raw-*")) for run in runnumbers]




if min != -1 and max != -1:
    begin = runnumbers.index(min)
    end = runnumbers.index(max)+1

    runN = runnumbers[begin:end]
    runF = runfiles[begin:end]

else :
    runN = runnumbers
    runF = runfiles

if len(runN) != len(runF) : 
    print('something is wrong')
    sys.exit(1)

if len(sys.argv) == 5: option = sys.argv[4]
else : option = 'w'

with open(f'./prova.txt', option) as f:
    for i in range(len(runN)): 
        f.write(f'{runN[i]} {runF[i]} {bool}\n')
    