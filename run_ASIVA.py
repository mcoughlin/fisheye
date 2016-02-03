#!/usr/bin/python

"""
% run_fisheye.py

Written by Chris Stubbs and Michael Coughlin (michael.coughlin@ligo.org)

This program analyzes data from the nightly Cerro Pachon All-Sky Camera Project camera.

"""

import os, sys, pickle, math, optparse, glob, time, subprocess
import asiva

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

    parser.add_option("-t", "--eventtime", help="time",default="1392172512")
    parser.add_option("-u", "--user", help="username",default="coughlin")
    parser.add_option("-o", "--outputFolder", help="output folder",default="1392172512")
    parser.add_option("-v", "--verbose", action="store_true", default=False,
                      help="Run verbosely. (Default: False)")

    parser.add_option("--doPlotFits",  action="store_true", default=False)
    parser.add_option("--doMakeMovie",  action="store_true", default=False)

    opts, args = parser.parse_args()

    return opts

def initialize(params):

    # set up directory names and data paths
    params["datadisk"]='/lsst/all-sky-ASIVA/ALOPachon'
    params["fitsdisk"]='/scratch/%s/allsky/data/ASIVA'%params["user"]
    params["codepath"]='/lsst/home/coughlin/git-repo/fisheye'

    params["dirpath"]='/lsst/home/%s/allsky/data/ASIVA'%params["user"]
    params["dirpath"]='/scratch/%s/allsky/data/ASIVA'%params["user"]
    params["dirname"]=params["outputFolder"]

    params["fitspath"]=os.path.join(params["fitsdisk"],params["dirname"])
    if not os.path.isdir(params["fitspath"]):
        os.mkdir(params["fitspath"])

    # create tonight's directory
    params["dirpathname"]=os.path.join(params["dirpath"],params["dirname"])
    if not os.path.isdir(params["dirpathname"]):
        os.mkdir(params["dirpathname"])

    # prefix for images, gets framecounter.cr2 appended
    params["imageprefix"]="%s."%params["dirname"]
    # interval after end of last image before starting next one, in seconds
    params["pauseinterval"]=1

    # fix this later
    #if ($spaceleft>20); then
    #       mail -s "allsky camera lots space!" stubbs@physics.harvard.edu
    #fi

    params["filters"] = ["filt1","filt2"]
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
    user = opts.user
    eventtime = opts.eventtime

    params = {}
    params["outputFolder"] = outputFolder
    params["user"] = user
    params["eventtime"] = eventtime

    params = initialize(params)

    # Create fits files from .gz files
    params["doPlotFits"] = opts.doPlotFits
    # Do analysis where we look for clouds
    params["doMakeMovie"] = opts.doMakeMovie

    if params["doPlotFits"]:
        asiva.plotfits(params)
    if params["doMakeMovie"]:
        asiva.makemovie(params)
