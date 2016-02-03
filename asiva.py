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
import ephem

import healpy as hp

import matplotlib
#matplotlib.rc('text', usetex=True)
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 16})
matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib.mlab import griddata

# =============================================================================
#
#                               DEFINITIONS
#
# =============================================================================

def plotfits(params):

    folders = glob.glob(os.path.join(params["datadisk"],"run_*"))
    for folder in folders:
        fitsfolders = glob.glob(os.path.join(folder,"1*"))
        for fitsfolder in fitsfolders:

            fitsfolderSplit = fitsfolder.split("/")
            thistime = fitsfolderSplit[-1]

            if thistime == params["eventtime"]:
                thisfitsfolder = fitsfolder
            

    plotdir = params["fitspath"]
    fitsfiles = glob.glob(os.path.join(thisfitsfolder,"*open*.fits"))
    for fitsfile in fitsfiles:

        hdulist1 = astropy.io.fits.open(fitsfile)
        hdulist2 = astropy.io.fits.open(fitsfile.replace("open","closed"))

        fitsSplit = fitsfile.split("/")
        eventtime = int(fitsSplit[-2])
        utctime = datetime.fromtimestamp(eventtime)

        data1 = hdulist1[0].data
        data2 = hdulist2[0].data
        data = data1 - data2

        data = data[0,:,:]
        data = np.abs(data)

        a,b = data.shape
        xmid = np.floor(b/2.0) - 10
        ymid = np.floor(a/2.0) + 20

        thetas = (2*np.pi/360.0)*np.arange(360)
        r = 185
        xs = r*np.cos(thetas) + xmid
        ys = r*np.sin(thetas) + ymid

        figurename = fitsSplit[-1]
        figurenameSplit = figurename.split("_")
        filt = figurenameSplit[1]

        vmin = 5000
        vmax = 10000

        plotName = os.path.join(plotdir,'%s.png'%filt)
        plt.figure()
        plt.imshow(data,vmin=vmin,vmax=vmax)
        #plt.plot(xs,ys,'r')
        plt.xlim([0,b])
        plt.ylim([0,a])
        plt.title(utctime)
        plt.show()
        plt.savefig(plotName,dpi=200)
        plt.close('all')

        xi = np.arange(b) - xmid
        yi = np.arange(a) - ymid
        XI,YI = np.meshgrid(xi,yi)
        RI = np.sqrt(XI**2 + YI**2)

        #indexes = np.where(RI[:] > r)[0]
        data_cut = data
        data_cut[RI > r] = 0.0
        data_cut = np.abs(data_cut)

        plotName = os.path.join(plotdir,'%s_cut.png'%filt)
        plt.figure()
        plt.imshow(data_cut,vmin=vmin,vmax=vmax)
        #plt.plot(xs,ys,'r')
        plt.xlim([0,b])
        plt.ylim([0,a])
        plt.title(utctime)
        plt.show()
        plt.savefig(plotName,dpi=200)
        plt.close('all')

        LAT=-30.67833
        LNG=-70.73669
        #HA=0.015
        #DEC=-27.092
        #AZ=95.082
        SCALE=0.04705

        xx = XI*SCALE + LNG
        yy = YI*SCALE + LAT
        zz = data_cut

        lng_min = np.min(xx)
        lng_max = np.max(xx)
        lat_min = np.min(yy)
        lat_max = np.max(yy)

        xx = xx.ravel()
        yy = yy.ravel()
        zz = zz.ravel()

        print xx.shape, yy.shape, zz.shape

        xi = np.linspace(-180,180,361)
        yi = np.linspace(-90,90,181)
        zi = griddata(xx,yy,zz,xi,yi,interp='nn')
        zi = np.array(zi)

        zi[np.isnan(zi)] = 0

        z_min = np.min(zi)
        z_max = np.max(zi)

        plotName = os.path.join(plotdir,'%s_pcolor.png'%filt)
        plt.pcolor(xi, yi, zi, cmap='RdBu', vmin=z_min, vmax=z_max)
        # plot data points.
        #plt.xlim([-180,180])
        #plt.ylim([-90,90])
        plt.title(utctime)
        plt.xlim([lng_min,lng_max])
        plt.ylim([lat_min,lat_max])
        cb = plt.colorbar()
        #cb.set_label('vals')
        #cb.set_label('(minst - m) - m0')
        plt.show()
        plt.savefig(plotName,dpi=200)
        plt.close('all')

        nside = 512
        print "Pixel area: %.4f square degrees" % hp.nside2pixarea(nside, degrees=True)

        phi = xx*2*np.pi/360.0
        theta = np.pi/2.0 - yy*2*np.pi/360.0
        indexes = hp.ang2pix(nside, theta, phi)
        healpix_map = np.zeros(hp.nside2npix(nside), dtype=np.double)
        healpix_map[indexes] = zz

        plotName = os.path.join(plotdir,'%s_mollview.png'%filt)
        plt.figure()
        hp.mollview(healpix_map)
        plt.show()
        plt.savefig(plotName,dpi=200)
        plt.close('all')

        hpfits = os.path.join(plotdir,'%s.fits'%filt)
        hp.write_map(hpfits, healpix_map)


def makemovie(params):

    filters = params["filters"]

    for filter in filters:

        #outpath = os.path.join(params["dirpathname"],filter)
        outpath = os.path.join(params["dirpathname"],filter)
        if not os.path.isdir(outpath):
            os.mkdir(outpath)
        moviedir = os.path.join(outpath,"movie")
        if not os.path.isdir(moviedir):
            os.mkdir(moviedir)

        n=1
        folders = sorted(glob.glob(os.path.join(params["fitsdisk"],'1*')))
        for folder in folders:
            file = os.path.join(folder,'%s_cut.png'%filter)
            if not os.path.isfile(file):
                continue
            n = n + 1
            filename = os.path.join(moviedir,"asiva-%04d.png"%n)
            cp_command = "cp %s %s"%(file,filename)
            os.system(cp_command)

            print file

        moviefiles = os.path.join(moviedir,"asiva-%04d.png")
        filename = os.path.join(moviedir,"asiva.mpg")
        ffmpeg_command = 'ffmpeg -an -y -r 20 -i %s -b:v %s %s'%(moviefiles,'5000k',filename)
        os.system(ffmpeg_command)
        filename = os.path.join(moviedir,"asiva.gif")
        ffmpeg_command = 'ffmpeg -an -y -r 20 -i %s -b:v %s %s'%(moviefiles,'5000k',filename)
        os.system(ffmpeg_command)

        rm_command = "rm %s/*.png"%(moviedir)
        os.system(rm_command)

