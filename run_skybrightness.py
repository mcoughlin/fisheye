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

    FITSfolder = "/lsst/all-sky/FITS" 
    folders = glob.glob(os.path.join(FITSfolder,"ut*"))
    folders = sorted(folders)

    flags = "--doSkyBrightness"

    for folder in folders:

        folderName = folder.split("/")[-1]
        fisheye_command = "python run_fisheye.py -f %s -o %s %s"%(folderName,folderName,flags)
        os.system(fisheye_command)

