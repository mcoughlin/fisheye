#!/usr/bin/python

"""
% startnight.py

Written by Chris Stubbs and Michael Coughlin (michael.coughlin@ligo.org)

This program runs the nightly Cerro Pachon All-Sky Camera Project camera.

"""

import os, time, glob
import numpy as np
import subprocess

# =============================================================================
#
#                                    MAIN
#
# =============================================================================

if __name__=="__main__":

    allskyFolder = "/lsst/all-sky"
    FITSfolder = "/lsst/all-sky/FITS" 
    folders = glob.glob(os.path.join(allskyFolder,"ut*"))
    folders = sorted(folders)

    flags = "--doMakeFits"
    n = 0

    for folder in folders:

        folderName = folder.split("/")[-1]

        folderDir = os.path.join(FITSfolder,folderName)
        if os.path.isdir(folderDir): 
            continue

        n = n + 1
        if n>2:
            continue

        fisheye_command = "python run_fisheye.py -f %s -o %s %s"%(folderName,folderName,flags)
        #os.system(fisheye_command)
        print fisheye_command

        proc = subprocess.Popen(fisheye_command, shell=True,
             stdin=None, stdout=None, stderr=None, close_fds=True)
        #os.spawnlp(os.P_NOWAIT, fisheye_command)
      
        #fisheye_command = "python"
        #fisheye_args = "run_fisheye.py -f %s -o %s %s"%(folderName,folderName,flags)

        #os.spawnlp(os.P_NOWAIT, fisheye_command,fisheye_args)

