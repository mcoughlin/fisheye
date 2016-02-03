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

    FITSfolder = "/scratch/coughlin/allsky/data" 
    folders = glob.glob(os.path.join(FITSfolder,"ut*"))
    folders = sorted(folders)

    flags = "--doCombineNight"

    for folder in folders:

        folderName = folder.split("/")[-1]

        file = "/lsst/home/coughlin/git-repo/fisheye/catalognight/M/%s.txt"%folderName
        if os.path.isfile(file):
            continue

        fisheye_command = "python run_fisheye.py -o %s %s"%(folderName,flags)
        os.system(fisheye_command)

