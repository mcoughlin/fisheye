
import os, glob

params = {}
params["datadisk"]='/lsst/all-sky-ASIVA/ALOPachon'

eventtimes = []

folders = glob.glob(os.path.join(params["datadisk"],"run_*"))
for folder in folders:
    fitsfolders = glob.glob(os.path.join(folder,"1*"))
    for fitsfolder in fitsfolders:

        fitsfolderSplit = fitsfolder.split("/")
        thistime = fitsfolderSplit[-1]

        eventtimes.append(thistime)

for eventtime in eventtimes:
    system_command = "python run_ASIVA.py --doPlotFits -t %s -o %s"%(eventtime,eventtime)
    os.system(system_command)


