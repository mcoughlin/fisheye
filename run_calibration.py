#!/usr/bin/python

"""
% run_fisheye.py

Written by Chris Stubbs and Michael Coughlin (michael.coughlin@ligo.org)

This program analyzes data from the nightly Cerro Pachon All-Sky Camera Project camera.

"""

import os, sys, pickle, math, optparse, glob, time, subprocess

import numpy as np
import scipy.ndimage
import astropy
import aplpy

from scipy.stats import mode
#from scipy.interpolate import griddata

import matplotlib
#matplotlib.rc('text', usetex=True)
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 16})
import matplotlib.pyplot as plt

from matplotlib.mlab import griddata


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
    parser.add_option("-e", "--exposureLength", help="exposure length",default="10")
    parser.add_option("-u", "--user", help="username",default="coughlin")
    parser.add_option("-f", "--folderName", help="folder name",default="d")
    parser.add_option("-o", "--outputFolder", help="output folder",default="d")
    parser.add_option("-v", "--verbose", action="store_true", default=False,
                      help="Run verbosely. (Default: False)")

    parser.add_option("--doCalibration",  action="store_true", default=False)

    opts, args = parser.parse_args()

    return opts

def plotfits(fitsfile,hdulist,plotdir):

    if not os.path.isdir(plotdir):
        os.mkdir(plotdir)

    plotName = os.path.join(plotdir,'fits.png')
    gc = aplpy.FITSFigure(fitsfile)
    gc.show_grayscale()
    gc.tick_labels.set_font(size='small')
    gc.save(plotName)
    gc.close()

    hist_vals = hdulist[0].data

    plotName = os.path.join(plotdir,'hist.png')
    bins = np.arange(-50,50,1)
    hist, bins = np.histogram(hist_vals, bins=bins)
    width = 0.7 * (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2
    plt.bar(center, hist, align='center', width=width)
    plt.show()
    plt.savefig(plotName,dpi=200)
    plt.close('all')

    hist_vals = hdulist[0].data[:100,:100]

    plotName = os.path.join(plotdir,'hist_100.png')
    bins = np.arange(-50,50,1)
    hist, bins = np.histogram(hist_vals, bins=bins)
    width = 0.7 * (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2
    plt.bar(center, hist, align='center', width=width)
    plt.show()
    plt.savefig(plotName,dpi=200)
    plt.close('all')

def calibration(params):

    filters = params["filters"]
    colors = params["colors"]

    for filter in filters:
        #outpath = os.path.join(params["dirpathname"],filter)
        outpath = os.path.join(params["datapath"],filter)

        plotdir = os.path.join(outpath,"plots")
        if not os.path.isdir(plotdir):
            os.mkdir(plotdir)

        files = glob.glob(os.path.join(outpath,"*.fits.old"))
        files = sorted(files)

        hdulistave = []

        n = 0.0
        for fitsfile in files:
        #for fitsfile in []:

            try:
                hdulist = astropy.io.fits.open(fitsfile)
            except:
                continue
            fitsSplit = fitsfile.split("/")

            n = n + 1.0

            #fitshdr_command = "fitshdr %s"%(fitsfile)
            #p = subprocess.Popen(fitshdr_command.split(" "),stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            #output, errors = p.communicate()

            bias = 2047.228

            hdulist[0].data = hdulist[0].data - bias

            if len(hdulistave) == 0:
                hdulistave = hdulist
                #hdulistmax = hdulist
                #hdulistmin = hdulist
                #hdulistave = [1]

            else:
                hdulistave[0].data = hdulistave[0].data.copy() + hdulist[0].data.copy()
                #hdulistmax[0].data = np.maximum(hdulistmax[0].data,hdulist[0].data.copy())

                #hdulistmin[0].data = np.minimum(hdulistmin[0].data,hdulist[0].data.copy())

            plotdir = os.path.join(outpath,"plots/%s"%fitsSplit[-1])
            plotfits(fitsfile,hdulist,plotdir)

            print plotdir
            print stop

        #hdulistave[0].data = hdulistave[0].data / n
        #hdulistmax[0].data = hdulistmax[0].data / n
        #hdulistmin[0].data = hdulistmin[0].data / n

        fitsfileave = os.path.join(plotdir,'ave.fits')
        #hdulistave.writeto(fitsfileave,output_verify='warn',clobber=True)
        plotdir = os.path.join(outpath,"plots/ave")
        #plotfits(fitsfileave,hdulistave,plotdir)

        fitsfilemax = os.path.join(plotdir,'max.fits')
        #hdulistmax.writeto(fitsfilemax,output_verify='warn',clobber=True)
        plotdir = os.path.join(outpath,"plots/max")
        #plotfits(fitsfilemax,hdulistmax,plotdir)

        fitsfilemin = os.path.join(plotdir,'min.fits')
        #hdulistmin.writeto(fitsfilemin,output_verify='warn',clobber=True)
        plotdir = os.path.join(outpath,"plots/min")
        #plotfits(fitsfilemin,hdulistmin,plotdir)

        print stop

        xmed = np.median(hdulistave[0].data,axis=0)
        ymed = np.median(hdulistave[0].data,axis=1)

        #xmed = np.mean(hdulistave[0].data,axis=0)
        #ymed = np.mean(hdulistave[0].data,axis=1)

        plotName = os.path.join(plotdir,'xmed.png')
        plt.plot(xmed,'*')
        plt.ylim([-3,3])
        plt.show()
        plt.savefig(plotName,dpi=200)
        plt.close('all')

        plotName = os.path.join(plotdir,'ymed.png')
        plt.plot(ymed,'*')
        plt.ylim([-3,3])
        plt.show()
        plt.savefig(plotName,dpi=200)
        plt.close('all')

        print "%s: %d"%(filter,len(files))
        continue

        #for fitsfile1 in files:
        for fitsfile1 in []:
            for fitsfile2 in files:

                hdulist1 = astropy.io.fits.open(fitsfile1)
                fitsSplit1 = fitsfile1.split("/")

                hdulist2 = astropy.io.fits.open(fitsfile2)
                fitsSplit2 = fitsfile2.split("/")

                if fitsSplit1[-1] == fitsSplit2[-1]:
                    continue

                hdulist1 = astropy.io.fits.open(fitsfile1)
                hdulist2 = astropy.io.fits.open(fitsfile2)

                fitsfilesubtract = os.path.join(plotdir,'subtract.fits')
                hdulist = astropy.io.fits.open(fitsfile1)
                hdulist[0].data = hdulist2[0].data - hdulist1[0].data 
                hdulist.writeto(fitsfilesubtract,output_verify='warn',clobber=True)

                plotdir = os.path.join(outpath,"plots/%s_%s"%(fitsSplit1[-1],fitsSplit2[-1]))
                plotfits(fitsfilesubtract,hdulistsubtract,plotdir)
            
def initialize(params):

    # set up directory names and data paths
    params["datadisk"]='/home/coughlin/allsky/calibdata/'
    params["dirpath"]='/lsst/home/%s/allsky/calibdata'%params["user"]
    params["dirname"]=params["outputFolder"]

    # create tonight's directory
    params["datapath"]=os.path.join(params["datadisk"],params["dirname"])
    if not os.path.isdir(params["datapath"]):
        os.mkdir(params["datapath"])

    # create tonight's directory
    params["datapath"]=os.path.join(params["datapath"],"%.5f"%eval(params["exposureLength"]))
    if not os.path.isdir(params["datapath"]):
        os.mkdir(params["datapath"])

    # create tonight's directory
    params["dirpathname"]=os.path.join(params["dirpath"],params["dirname"])
    if not os.path.isdir(params["dirpathname"]):
        os.mkdir(params["dirpathname"])

    # create tonight's directory
    params["dirpathname"]=os.path.join(params["dirpathname"],"%.5f"%eval(params["exposureLength"]))
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
    exposureLength = opts.exposureLength

    params = {}
    params["folderName"] = folderName
    params["outputFolder"] = outputFolder
    params["maxframes"] = maxframes
    params["user"] = user
    params["exposureLength"] = exposureLength

    # Create fits files from .gz files
    params["doCalibration"] = opts.doCalibration

    params = initialize(params)

    if params["doCalibration"]:
        calibration(params)
