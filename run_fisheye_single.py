#!/usr/bin/python

"""
% run_fisheye.py

Written by Chris Stubbs and Michael Coughlin (michael.coughlin@ligo.org)

This program analyzes data from the nightly Cerro Pachon All-Sky Camera Project camera.

"""

import os, sys, pickle, math, optparse, glob, time, subprocess
import fisheye_single
import fisheye

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

    parser.add_option("-f", "--file", help="file name",\
        default="/Volumes/ALLSKY1/ut042914/ut042914.0500.long.cr2.gz")
    #    default="/Volumes/2TB/ut062114/ut062114.0074.long.cr2.gz")
    parser.add_option("-u", "--user", help="username",default="coughlin")
    #parser.add_option("-o", "--outputFolder", help="output folder",default="webserver")
    parser.add_option("-o", "--outputFolder", help="output folder",default="ws")
    parser.add_option("-v", "--verbose", action="store_true", default=False,
                      help="Run verbosely. (Default: False)")

    parser.add_option("--doMakeFits",  action="store_true", default=False)
    parser.add_option("--doClouds",  action="store_true", default=False)
    parser.add_option("--doFisheye",  action="store_true", default=False)
    parser.add_option("--doSkyBrightness",  action="store_true", default=False)
    parser.add_option("--doCopyFiles",  action="store_true", default=False)
    parser.add_option("--doCopyFilesSingle",  action="store_true", default=False)

    parser.add_option("--doGetLatest",  action="store_true", default=False)

    opts, args = parser.parse_args()

    return opts

def getLatest(params):

    params["datadisk"]='/Volumes/2TB/'

    folders = sorted(glob.glob(os.path.join(params["datadisk"],"ut*")))
    folder = folders[-1]
    files = sorted(glob.glob(os.path.join(folder,"*.gz")))
    file = files[-1]
     
    params["file"] = file

    return params

def initialize(params):

    # set up directory names and data paths
    params["datadisk"]='/Users/christopherstubbs/allsky'
    params["datadisk"]='/Volumes/2TB/'
    params["catalogpath"] = "/Users/christopherstubbs/git-repo/fisheye/catalog"
    params["dirname"]=params["outputFolder"]

    # create tonight's directory
    params["datapath"]=os.path.join(params["datadisk"],params["dirname"])
    if not os.path.isdir(params["datapath"]):
        os.mkdir(params["datapath"])

    fileSplit = params["file"].split("/")[-1].split(".")
    fileCombine = "%s_%s"%(fileSplit[0],fileSplit[1])

    # create tonight's directory
    params["datapath"]=os.path.join(params["datapath"],fileCombine)
    if not os.path.isdir(params["datapath"]):
        os.mkdir(params["datapath"])

    params["fitspath"]=params["datapath"]
    params["dirpathname"]=params["datapath"]
    params["folderName"] = fileCombine

    params["fisheyeplotpath"]=os.path.join(params["datapath"],"fisheyeplots")
    if not os.path.isdir(params["fisheyeplotpath"]):
        os.mkdir(params["fisheyeplotpath"])

    params["cloudsplotpath"]=os.path.join(params["datapath"],"cloudsplots")
    if not os.path.isdir(params["cloudsplotpath"]):
        os.mkdir(params["cloudsplotpath"])

    params["skybrightnessplotpath"]=os.path.join(params["dirpathname"],"skybrightnessplots")
    if not os.path.isdir(params["skybrightnessplotpath"]):
        os.mkdir(params["skybrightnessplotpath"])

    params["maxframes"] = 10000

    cp_command = "cp %s %s"%(params["file"],params["datapath"])
    os.system(cp_command)

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
    file = opts.file
    user = opts.user

    params = {}
    params["outputFolder"] = outputFolder
    params["file"] = file
    params["user"] = user

    # Create fits files from .gz files
    params["doMakeFits"] = opts.doMakeFits
    # Do analysis where we look for clouds
    params["doClouds"] = opts.doClouds
    # Calculate sky brightness
    params["doSkyBrightness"] = opts.doSkyBrightness
    # Run fisheye code, which converts from pixel x,y to RA and Dec
    params["doFisheye"] = opts.doFisheye
    # Copy files to NCSA
    params["doCopyFiles"] = opts.doCopyFiles
    # Copy files to NCSA
    params["doCopyFilesSingle"] = opts.doCopyFilesSingle
    # Get Latest file
    params["doGetLatest"] = opts.doGetLatest

    if params["doGetLatest"]:
        params = getLatest(params)

    params = initialize(params)

    if params["doMakeFits"]:
        fisheye_single.makefits(params)
    if params["doFisheye"]:
        fisheye_single.fisheye(params)
    if params["doSkyBrightness"]:
        fisheye_single.skybrightness(params)
    if params["doClouds"]:
        fisheye_single.clouds(params)
    if params["doCopyFiles"]:
        fisheye_single.copy_files(params)
    if params["doCopyFilesSingle"]:
        fisheye_single.copy_files_single(params)
