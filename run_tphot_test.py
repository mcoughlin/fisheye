
import os
import numpy as np

import matplotlib
#matplotlib.rc('text', usetex=True)
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 16})
import matplotlib.pyplot as plt

outputdir = "/lsst/home/coughlin/allsky/data/ut012914/tphottest" 
plotdir = os.path.join(outputdir,"plots")
if not os.path.isdir(plotdir):
    os.mkdir(plotdir)

aprads = np.arange(1,20)
low_fluxes = []
high_fluxes = []
data_out = {}

for aprad in aprads: 
    file = "%s/files/%d.tph"%(outputdir,aprad)
    tphot_command = "tphot /lsst/home/coughlin/allsky/data/FITS/ut012914/M/ut012914.0300.long.M.fits -bias 2048 -rad 4 -sig 2 -min 50 -aprad %d -okfit 1 -chin 100000 -npar 7 -out %s"%(aprad,file)
    #os.system(tphot_command)

    data_out_tph = np.loadtxt(file)

    data_out["%d"%aprad] = {}
    data_out["%d"%aprad]["data"] = data_out_tph
    data_out["%d"%aprad]["flux"] = data_out_tph[:,7]
    data_out["%d"%aprad]["snr"] = data_out_tph[:,7] / data_out_tph[:,8]
 
    high_fluxes.append(data_out["%d"%aprad]["flux"][3])
    low_fluxes.append(data_out["%d"%aprad]["flux"][1])

    print data_out_tph[1000,0], data_out_tph[1000,1], len(data_out["%d"%aprad]["flux"])

print stop

for i in xrange(len(data_out["1"]["flux"])):

    snrs = []
    fluxes = []
    for key in data_out.iterkeys():
        fluxes.append(data_out[key]["flux"][i])
        snrs.append(data_out[key]["snr"][i])

    if snrs[0] < 50:
        continue

    plotName = os.path.join(plotdir,'fluxes_%d.png'%i)
    plt.plot(aprads,fluxes,"*")
    plt.show()
    plt.savefig(plotName,dpi=200)
    plt.close('all')

plotName = os.path.join(plotdir,'high_fluxes.png')
plt.plot(aprads,high_fluxes,"*")
plt.show()
plt.savefig(plotName,dpi=200)
plt.close('all')

plotName = os.path.join(plotdir,'low_fluxes.png')
plt.plot(aprads,low_fluxes,"*")
plt.show()
plt.savefig(plotName,dpi=200)
plt.close('all')

