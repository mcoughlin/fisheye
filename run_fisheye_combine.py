#!/usr/bin/python

"""
% startnight.py

Written by Chris Stubbs and Michael Coughlin (michael.coughlin@ligo.org)

This program runs the nightly Cerro Pachon All-Sky Camera Project camera.

"""

from __future__ import division

import os, sys, pickle, math, optparse, glob, time, subprocess
from datetime import date, datetime
import numpy as np
import scipy.ndimage
import astropy
import aplpy

from scipy.stats import mode

import matplotlib
#matplotlib.rc('text', usetex=True)
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 16})
import matplotlib.pyplot as plt

# =============================================================================
#
#                               DEFINITIONS
#
# =============================================================================

def parse_commandline():
    """
    Parse the options given on the command-line.
    """
    parser = optparse.OptionParser()

    parser.add_option("-f", "--folderName", help="folder name",default="ut")
    parser.add_option("-o", "--outputFolder", help="output folder",default="fisheye_combined")
    parser.add_option("-v", "--verbose", action="store_true", default=False,
                      help="Run verbosely. (Default: False)")

    parser.add_option("--doCompareSkyBrightness",  action="store_true", default=False)
    parser.add_option("--doCompareFisheye",  action="store_true", default=False)

    opts, args = parser.parse_args()

    return opts

def compareskybrightness(params):

    filters = params["filters"]
    colors = params["colors"]

    for filter,color in zip(filters,colors):

        outpath = os.path.join(params["skybrightnessplotpath"],filter)
        if not os.path.isdir(outpath):
            os.mkdir(outpath)
        moviedir = os.path.join(outpath,"movie")
        if not os.path.isdir(moviedir):
            os.mkdir(moviedir)

        dirpath = os.path.join(params["dirpath"],"%s*"%params["folderName"])
        folders = glob.glob(dirpath)
        folders = sorted(folders)

        n=1
        for folder in folders:
            file = os.path.join(folder,"skybrightnessplots/%s/skybrightness.png"%filter)
            if not os.path.isfile(file):
                continue

            n = n + 1
            filename = os.path.join(moviedir,"fishy-%04d.png"%n)
            cp_command = "cp %s %s"%(file,filename)
            os.system(cp_command)

        moviefiles = os.path.join(moviedir,"fishy-%04d.png")
        filename = os.path.join(moviedir,"long.mpg")
        ffmpeg_command = 'ffmpeg -an -y -r 20 -i %s -b:v %s %s'%(moviefiles,'5000k',filename)
        os.system(ffmpeg_command)
        filename = os.path.join(moviedir,"long.gif")
        ffmpeg_command = 'ffmpeg -an -y -r 20 -i %s -b:v %s %s'%(moviefiles,'5000k',filename)
        os.system(ffmpeg_command)

        rm_command = "rm %s/*.png"%(moviedir)
        os.system(rm_command)

def comparefisheye(params):

    filters = params["filters"]
    colors = params["colors"]

    for filter in filters:

        plotDir = os.path.join(params["timeplotpath"],filter)
        if not os.path.isdir(plotDir):
            os.mkdir(plotDir)

        dirpath = os.path.join(params["dirpath"],"%s*"%params["folderName"])
        folders = glob.glob(dirpath)
        folders = sorted(folders)

        data_all = []
        data_date = []
        data_ave = []

        for folder in folders:

            magsDir = os.path.join(folder,"fisheyeplots/%s"%filter)
            filename = os.path.join(magsDir,"mags.txt")
            if not os.path.isfile(filename):
                continue

            data_out = np.loadtxt(filename)
            if len(data_out) == 0:
                continue

            if len(data_all) == 0:
                data_all = data_out
            else:
                data_all = np.vstack((data_all,data_out))

            if len(data_out) < 10:
                continue

            folderName = folder.split("/")[-1]
            data_date.append(folderName)
            data_ave.append(np.median(data_out[:,2]))

        plotName = os.path.join(plotDir,'std.png')
        unique_mag = np.unique(data_all[:,0])
        for i in xrange(len(unique_mag)):
            index = unique_mag[i]

            indexes = np.where(data_all[:,0] == index)[0]

            #if len(indexes) < 10:
            #    continue

            plt.semilogy(index,np.std(data_all[indexes,2]),'*')

        plt.xlim([0,360])
        #plt.ylim([11,15])
        plt.xlabel('')
        plt.ylabel('std(m - minst)')
        #cb = plt.colorbar()
        #cb.set_label('minst')
        #cb.set_label('minst - m')
        plt.show()
        plt.savefig(plotName,dpi=200)
        plt.close('all')

        plotName = os.path.join(plotDir,'mag.png')
        plt.plot(data_all[:,0],data_all[:,2],'*')
        plt.xlim([0,360])
        plt.ylim([11,15])
        plt.xlabel('')
        plt.ylabel('m - minst')
        #cb = plt.colorbar()
        #cb.set_label('minst')
        #cb.set_label('minst - m')
        plt.show()
        plt.savefig(plotName,dpi=200)
        plt.close('all')

        plotName = os.path.join(plotDir,'ave.png')
        fig = plt.gcf()
        fig.set_size_inches(10,12)
        plt.plot(data_ave,'*')
        plt.xlabel('date')
        plt.ylabel('m - minst')
        locs, labels = plt.xticks()
        locs = np.arange(0,len(data_ave))
        plt.xticks(locs, data_date, rotation=90)
        plt.show()
        plt.savefig(plotName,dpi=200)
        plt.close('all')

def initialize(params):

    # set up directory names and data paths
    params["datadisk"]='/lsst/all-sky'
    params["dirpath"]='/lsst/home/coughlin/allsky/data'
    params["dirname"]=params["outputFolder"]

    # create tonight's directory
    #params["datapath"]=os.path.join(params["datadisk"],params["dirname"])
    #if not os.path.isdir(params["datapath"]):
    #    os.mkdir(params["datapath"])

    # create tonight's directory
    params["dirpathname"]=os.path.join(params["dirpath"],params["dirname"])
    if not os.path.isdir(params["dirpathname"]):
        os.mkdir(params["dirpathname"])

    # create tonight's directory
    params["expathname"]=os.path.join(params["dirpathname"],"extractor")
    if not os.path.isdir(params["expathname"]):
        os.mkdir(params["expathname"])

    # create tonight's directory
    params["imagepathname"]=os.path.join(params["dirpath"],params["dirname"])
    if not os.path.isdir(params["imagepathname"]):
        os.mkdir(params["imagepathname"])

    params["fitsplotpath"]=os.path.join(params["dirpathname"],"fitsplots")
    if not os.path.isdir(params["fitsplotpath"]):
        os.mkdir(params["fitsplotpath"])

    params["fisheyeplotpath"]=os.path.join(params["dirpathname"],"fisheyeplots")
    if not os.path.isdir(params["fisheyeplotpath"]):
        os.mkdir(params["fisheyeplotpath"])

    params["cloudsplotpath"]=os.path.join(params["dirpathname"],"cloudsplots")
    if not os.path.isdir(params["cloudsplotpath"]):
        os.mkdir(params["cloudsplotpath"])

    params["timeplotpath"]=os.path.join(params["dirpathname"],"timeplots")
    if not os.path.isdir(params["timeplotpath"]):
        os.mkdir(params["timeplotpath"])

    params["skybrightnessplotpath"]=os.path.join(params["dirpathname"],"skybrightnessplots")
    if not os.path.isdir(params["skybrightnessplotpath"]):
        os.mkdir(params["skybrightnessplotpath"])


    # prefix for images, gets framecounter.cr2 appended
    params["imageprefix"]="%s."%params["dirname"]
    # interval after end of last image before starting next one, in seconds
    params["pauseinterval"]=1

    # fix this later
    #if ($spaceleft>20); then
    #       mail -s "allsky camera lots space!" stubbs@physics.harvard.edu
    #fi

    params["filters"] = ["M","B","G","R"]
    params["colors"] = ["bw","g1","b","g2"]

    return params

# =============================================================================
#
#                                    MAIN
#
# =============================================================================

if __name__=="__main__":

    opts = parse_commandline()

    outputFolder = opts.outputFolder
    folderName = opts.folderName
    params = {}

    params["outputFolder"] = outputFolder
    params["folderName"] = folderName

    params = initialize(params)

    # Compare fisheye fits across nights
    params["doCompareFisheye"] = opts.doCompareFisheye
    # Compare sky brightness across nights
    params["doCompareSkyBrightness"] = opts.doCompareSkyBrightness

    if params["doCompareSkyBrightness"]:
        compareskybrightness(params)
    if params["doCompareFisheye"]:
        comparefisheye(params)

