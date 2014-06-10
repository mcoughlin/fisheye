#!/usr/bin/python

"""
% run_fisheye.py

Written by Chris Stubbs and Michael Coughlin (michael.coughlin@ligo.org)

This program analyzes data from the nightly Cerro Pachon All-Sky Camera Project camera.

"""

import os, sys, pickle, math, optparse, glob, time, subprocess
import numpy as np
import calendar

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

    parser.add_option("-u", "--user", help="username",default="coughlin")
    parser.add_option("-s", "--site", help="site",default="Arica")
    parser.add_option("-v", "--verbose", action="store_true", default=False,
                      help="Run verbosely. (Default: False)")

    parser.add_option("--doAeronet",  action="store_true", default=False)

    opts, args = parser.parse_args()

    return opts

def aeronet(params):

    files = glob.glob(os.path.join(params["datadisk"],"*%s.lev*"%params["site"]))
    files = sorted(files)

    timestamps = []
    data = []

    for file in files:
        lines = [line.strip() for line in open(file)]

        for line in lines:
            lineSplit = line.split(",")
            if len(lineSplit) < 10:
                continue
            if lineSplit[0] == "Date(dd-mm-yy)":
                continue

            day = lineSplit[0]
            daytime = lineSplit[1]
            strtimestamp = "%s %s"%(day.replace(":","/"),daytime)

            timestamp = time.strptime(strtimestamp, "%d/%m/%Y %H:%M:%S")  
            #timestamps.append(timestamp)
            timestamps.append(calendar.timegm(timestamp))

            vals = lineSplit[3:-2]
            for i in xrange(len(vals)):
                if vals[i] == "N/A":
                    vals[i] = float('nan') 
                else:
                    vals[i] = float(vals[i])
        
            if len(data) == 0:
                data = vals
            else:
                data = np.vstack((data,vals))      

    data = data.T
    data[np.isnan(data)] = 0

    indexes = []
    for i in xrange(len(vals)):
        if np.sum(data[i,:]) == 0:
            continue    
        indexes.append(i)

    # Date(dd-mm-yy),Time(hh:mm:ss),Julian_Day,AOT_1640,AOT_1020,AOT_870,AOT_675,AOT_667,AOT_555,AOT_551,AOT_532,AOT_531,AOT_500,AOT_490,AOT_443,AOT_440,AOT_412,AOT_380,AOT_340,Water(cm),%TripletVar_1640,%TripletVar_1020,%TripletVar_870,%TripletVar_675,%TripletVar_667,%TripletVar_555,%TripletVar_551,%TripletVar_532,%TripletVar_531,%TripletVar_500,%TripletVar_490,%TripletVar_443,%TripletVar_440,%TripletVar_412,%TripletVar_380,%TripletVar_340,%WaterError,440-870Angstrom,380-500Angstrom,440-675Angstrom,500-870Angstrom,340-440Angstrom,440-675Angstrom(Polar),Last_Processing_Date(dd/mm/yyyy),Solar_Zenith_Angle

    indexes = np.arange(32,39)
    data = data[indexes,:]

    x = timestamps
    y = np.arange(len(data))
    ylabels = ["440-870","380-500","440-675","500-870","340-440","440-675",""]

    X, Y = np.meshgrid(x,y)
 
    plotName = os.path.join(params["aeronetplotpath"],'wavelength.png')
    fig = plt.gcf()
    fig.set_size_inches(14,18)

    plt.pcolor(X,Y,data)
    #plt.xlim([0,2897])
    #plt.ylim([0,1935])
    #plt.title("average loss")
    locs, labels = plt.xticks()
    #locs = np.arange(0,len(data_ave))
    plt.xticks(rotation=90)
    locs, labels = plt.yticks()
    plt.yticks(y+0.5, ylabels)

    cb = plt.colorbar()
    cb.set_label('amp')
    plt.show()
    plt.savefig(plotName,dpi=200)
    plt.close('all')

def initialize(params):

    # set up directory names and data paths
    params["datadisk"]='/lsst/home/%s/Aeronet/data'%params["user"]
    params["dirpath"]='/lsst/home/%s/Aeronet'%params["user"]
    params["dirname"]=params["site"]

    # create tonight's directory
    params["dirpathname"]=os.path.join(params["dirpath"],params["dirname"])
    if not os.path.isdir(params["dirpathname"]):
        os.mkdir(params["dirpathname"])

    params["aeronetplotpath"]=os.path.join(params["dirpathname"],"aeronetplots")
    if not os.path.isdir(params["aeronetplotpath"]):
        os.mkdir(params["aeronetplotpath"])

    # prefix for images, gets framecounter.cr2 appended
    params["imageprefix"]="%s."%params["dirname"]
    # interval after end of last image before starting next one, in seconds
    params["pauseinterval"]=1

    # fix this later
    #if ($spaceleft>20); then
    #       mail -s "allsky camera lots space!" stubbs@physics.harvard.edu
    #fi

    return params

# =============================================================================
#
#                                    MAIN
#
# =============================================================================

if __name__=="__main__":

    opts = parse_commandline()

    site = opts.site
    user = opts.user

    params = {}
    params["site"] = site
    params["user"] = user

    # Plot Aeronet data
    params["doAeronet"] = opts.doAeronet

    params = initialize(params)

    if params["doAeronet"]:
        aeronet(params)
