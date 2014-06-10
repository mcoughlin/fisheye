#!/usr/bin/python

"""
% run_fisheye.py

Written by Chris Stubbs and Michael Coughlin (michael.coughlin@ligo.org)

This program analyzes data from the nightly Cerro Pachon All-Sky Camera Project camera.

"""

import os, sys, pickle, math, optparse, glob, time, subprocess
import fisheye, photodiode

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

    parser.add_option("-m", "--maxframes", help="max frames",default=100000,type=int)
    parser.add_option("-u", "--user", help="username",default="coughlin")
    parser.add_option("-f", "--folderName", help="folder name",default="ut041414")
    parser.add_option("-o", "--outputFolder", help="output folder",default="ut041414")
    parser.add_option("-v", "--verbose", action="store_true", default=False,
                      help="Run verbosely. (Default: False)")

    parser.add_option("--doMakeFits",  action="store_true", default=False)
    parser.add_option("--doMakeMovie",  action="store_true", default=False)
    parser.add_option("--doClouds",  action="store_true", default=False)
    parser.add_option("--doFisheye",  action="store_true", default=False)
    parser.add_option("--doCombineFisheye",  action="store_true", default=False)
    parser.add_option("--doPhotometry",  action="store_true", default=False)
    parser.add_option("--doSkyBrightness",  action="store_true", default=False)
    parser.add_option("--doXY2Sky",  action="store_true", default=False)
    parser.add_option("--doHealpix",  action="store_true", default=False)
    parser.add_option("--doPhotodiode",  action="store_true", default=False)

    opts, args = parser.parse_args()

    return opts

def initialize(params):

    # set up directory names and data paths
    params["datadisk"]='/lsst/all-sky'
    params["fitsdisk"]='/lsst/all-sky/FITS'    
    params["fitsdisk"]='/lsst/home/coughlin/allsky/data/FITS'
    params["photodiodedisk"]='/lsst/home/coughlin/allsky/data/PD'
    params["dirpath"]='/lsst/home/%s/allsky/data'%params["user"]
    params["dirname"]=params["outputFolder"]

    # create tonight's directory
    params["datapath"]=os.path.join(params["datadisk"],params["dirname"])
    if not os.path.isdir(params["datapath"]):
        os.mkdir(params["datapath"])

    params["fitspath"]=os.path.join(params["fitsdisk"],params["dirname"])
    if not os.path.isdir(params["fitspath"]):
        os.mkdir(params["fitspath"])

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

    params["xy2skyplotpath"]=os.path.join(params["dirpathname"],"xy2skyplots")
    if not os.path.isdir(params["xy2skyplotpath"]):
        os.mkdir(params["xy2skyplotpath"])

    params["skybrightnessplotpath"]=os.path.join(params["dirpathname"],"skybrightnessplots")
    if not os.path.isdir(params["skybrightnessplotpath"]):
        os.mkdir(params["skybrightnessplotpath"])

    params["healpixpath"]=os.path.join(params["dirpathname"],"healpixplots")
    if not os.path.isdir(params["healpixpath"]):
        os.mkdir(params["healpixpath"])

    params["photodiodepath"]=os.path.join(params["dirpathname"],"photodiodeplots")
    if not os.path.isdir(params["photodiodepath"]):
        os.mkdir(params["photodiodepath"])

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
    maxframes = opts.maxframes
    user = opts.user

    params = {}
    params["folderName"] = folderName
    params["outputFolder"] = outputFolder
    params["maxframes"] = maxframes
    params["user"] = user

    # Create fits files from .gz files
    params["doMakeFits"] = opts.doMakeFits
    # Do analysis where we look for clouds
    params["doClouds"] = opts.doClouds
    # Calculate sky brightness
    params["doSkyBrightness"] = opts.doSkyBrightness
    # Do analysis where we convert x,y to RA/Dec
    params["doXY2Sky"] = opts.doXY2Sky
    # Run fisheye code, which converts from pixel x,y to RA and Dec
    params["doFisheye"] = opts.doFisheye
    # Track stars across the night
    params["doCombineFisheye"] = opts.doCombineFisheye
    # Create movie of fits images
    params["doMakeMovie"] = opts.doMakeMovie
    # Run tphot / source extractor on the images
    params["doPhotometry"] = opts.doPhotometry
    # Convert images to healpix maps
    params["doHealpix"] = opts.doHealpix
    # Retrieve photodiode data
    params["doPhotodiode"] = opts.doPhotodiode

    params = initialize(params)

    if params["doMakeFits"]:
        fisheye.makefits(params)
    if params["doMakeMovie"]:
        fisheye.makemovie(params)
    if params["doPhotometry"]:
        fisheye.photometry(params)
        fisheye.photometry_tphot(params)
    if params["doFisheye"]:
        fisheye.fisheye(params)
    if params["doSkyBrightness"]:
        fisheye.skybrightness(params)
    if params["doCombineFisheye"]:
        fisheye.combinefisheye(params)
    if params["doClouds"]:
        fisheye.clouds(params)
    if params["doXY2Sky"]:
        fisheye.xy2sky(params)
    if params["doHealpix"]:
        fisheye.healpix(params)
    if params["doPhotodiode"]:
        photodiode.photodiode(params)

