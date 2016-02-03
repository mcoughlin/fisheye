#!/usr/bin/python

"""
% startnight.py

Written by Chris Stubbs and Michael Coughlin (michael.coughlin@ligo.org)

This program runs the nightly Cerro Pachon All-Sky Camera Project camera.

"""

import os, time, glob
import numpy as np


# =============================================================================
#
#                                    MAIN
#
# =============================================================================

if __name__=="__main__":

    FITSfolder = "/scratch/coughlin/allsky/data/" 
    folders = glob.glob(os.path.join(FITSfolder,"ut*"))
    folders = sorted(folders)

    filters = ["M","B","G","R"]

    for filter in filters:

        data = np.array([])
        for folder in folders:
 
            file = "%s/fisheyeplots/%s/mags.txt"%(folder,filter)

            if os.path.isfile(file):
                data_out = np.loadtxt(file)
                if len(data_out) == 0:
                    continue
 
                if len(data) == 0:
                    data = data_out
                else:
                    data = np.concatenate((data,data_out))
        data = data[:,:2]
        new_array = [tuple(row) for row in data]
        data = np.unique(new_array)
        data = data[data[:,0].argsort()]

        catalogfile = "/home/coughlin/git-repo/fisheye/catalog/%s.txt"%filter
        f = open(catalogfile,'w')
        for row in data:
            f.write('%.5f %.5f\n'%(row[0],row[1]))
        f.close() 
 
