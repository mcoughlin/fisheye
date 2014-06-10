#!/usr/bin/python

"""
% photodiode.py

Written by Jamison Frost Burke and Michael Coughlin (michael.coughlin@ligo.org)

This program runs the nightly Cerro Pachon All-Sky Camera Project camera.

"""

from __future__ import division

import os, sys, pickle, math, optparse, glob, time, subprocess
from datetime import date, datetime, timedelta
import numpy as np
import scipy.ndimage
import astropy
import aplpy

from astropy.time import Time

import healpy as hp

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

def readPhotoDiode(fn):

        arr = np.loadtxt(fn, dtype=[('date', 'S8'), ('time', 'S8'), ('vals', '8f')])

        fmt = '%y/%m/%d %H:%M:%S'
        t = [datetime.strptime('%s %s' % (row['date'], row['time']), fmt) for row in arr]

        return t, arr['vals']

def photodiode(params):

    files = glob.glob(os.path.join(params["photodiodedisk"],'*.txt'))
    mjdobss = []
    vals = []

    datestr = params["folderName"].replace("ut","")
    year = int("20%s"%datestr[4:6])
    day = int(datestr[2:4]) 
    month = int(datestr[0:2])

    day1 = datetime(year, month, day) - timedelta(0.5)
    day2 = day1 + timedelta(1)

    day1mjd = Time(day1, scale='utc').mjd
    day2mjd = Time(day2, scale='utc').mjd

    for file in files:
        times, values = readPhotoDiode(file)
        values = values.T

        t = Time(times, scale='utc')

        if len(mjdobss) == 0:
            mjdobss = t.mjd
            vals = values
        else:
            mjdobss = np.hstack((mjdobss,t.mjd))
            vals = np.hstack((vals,values))

    indexes = np.intersect1d(np.where(mjdobss >= day1mjd)[0],np.where(mjdobss <= day2mjd)[0])
    mjdobss = mjdobss[indexes]
    vals = vals[:,indexes]

    indexes = np.argsort(mjdobss)
    mjdobss = mjdobss[indexes]
    vals = vals[:,indexes]

    data = {}
    filters = {'R':0,'Z':2,'Y':4}
    for filter,index in filters.items():
        plotDir = os.path.join(params["photodiodepath"],filter)
        if not os.path.isdir(plotDir):
            os.mkdir(plotDir)

        plotName = os.path.join(plotDir,'photodiode.png')
        plt.plot(mjdobss,vals[index,:],'b*')
        plt.xlabel('time')
        plt.ylabel(filter)
        plt.title(params["folderName"])
        plt.show()
        plt.savefig(plotName,dpi=200)
        plt.close('all')

        f = open(os.path.join(plotDir,'photodiode.txt'),'w')
        for mjdobs,sky in zip(mjdobss,vals[index,:]):
            f.write("%.10f %.10e\n"%(mjdobs,sky))
        f.close()

        data[filter] = {}
        data[filter]["mjdobs"] = mjdobss
        data[filter]["skybrightness"] = vals[index,:]

    f = open(os.path.join(params["photodiodepath"],'photodiode.txt'),'w')
    for mjdobs,r,z,y in zip(data["R"]["mjdobs"],data["R"]["skybrightness"],data["Z"]["skybrightness"],data["Y"]["skybrightness"]):
        f.write("%.10f %.10e %.10e %.10e\n"%(mjdobs,r,z,y))
    f.close()

    plotName = os.path.join(params["photodiodepath"],'photodiode.png')
    plt.plot(data["R"]["mjdobs"],data["R"]["skybrightness"],'r*',label='R')
    plt.plot(data["Z"]["mjdobs"],data["Z"]["skybrightness"],'g*',label='Z')
    plt.plot(data["Y"]["mjdobs"],data["Y"]["skybrightness"],'b*',label='Y')
    plt.xlabel('time')
    plt.title(params["folderName"])
    plt.legend()
    plt.show()
    plt.savefig(plotName,dpi=200)
    plt.close('all')


