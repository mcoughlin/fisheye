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

    outFolder = 
    skybrightnessfiles = []

    for folder in folders:
        filename = os.path.join(folder,'skybrightnessplots/skybrightness.txt')
        folderName = folder.split("/")[-1]

        print filename
        print folderName
        skybrightnessfiles.append(filename)



