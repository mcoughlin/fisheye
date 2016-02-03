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

from scipy.stats import mode
#from scipy.interpolate import griddata
from scipy.optimize import curve_fit
from scipy.spatial import Delaunay

from skimage.transform import radon, rescale
from skimage.feature import blob_dog, blob_log, blob_doh

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

import utils

import ephem
from utils import calcDist_haversine, mjd2djd

from grabESO import grabESO

# =============================================================================
#
#                               DEFINITIONS
#
# =============================================================================

def esocompare(params):

    filters = params["filters"]
    colors = params["colors"]

    skybrightnessfile = os.path.join(params["skybrightnessplotpath"],'skybrightness.txt')
    #nightData = readChris('../data/'+front+ut+tag)

    def readChris(filename):
        # Looks like file goes Date M,R,G,B
        names=['time','M','R','G','B']
        types=[float]*5
        data = np.loadtxt(filename, dtype=zip(names,types))
        return data

    nightData = readChris(skybrightnessfile)

    sbs=[]
    lsstObs = ephem.Observer()
    lsstObs.lat = -0.527868529 #radians of '-30:14:40.7'
    lsstObs.lon = -1.2348102646986 #radians of '-70:44:57.9'
    lsstObs.elevation = 2662.75 #meters
    # Set the target alt and az
    targetAlt = 90.
    targetAz = 0.
    airmass = 1./np.cos(np.radians(90.-targetAlt) )
    esoobj = grabESO(restoreFile='')

    # Convert from MJD to DJD.
    nightData['time'] = mjd2djd(nightData['time'])

    # Compute when the 1/3, 2/3, 3/3 points of the night are.
    lsstObs.date = nightData['time'][0]
    sun = ephem.Sun(lsstObs)
    midnight = lsstObs.next_antitransit(sun)
    lsstObs.date = midnight
    nightStart = lsstObs.previous_setting(sun)
    nightEnd = lsstObs.next_rising(sun)
    night1_3 = (nightEnd-nightStart)/3. + nightStart
    night2_3 = (nightEnd-nightStart)/3.*2 + nightStart
    sunAlts = []
    moonAlts = []
    moonSeps = []

    for time in nightData['time']:

        # In the first third of the night
        nightTime = 1
        # In the 2nd third of the night
        if time > night1_3:
            nightTime = 2
        # In the final third of the night
        if time > night2_3:
            nightTime = 3

        lsstObs.date = time

        # Compute the ecliptic lat and lon for the pointing
        targetRA, targetDec = lsstObs.radec_of(targetAz, np.radians(targetAlt))
        ecl = ephem.Ecliptic(targetRA, targetDec)
        eclLon = np.degrees(ecl.lon)
        eclLat = np.degrees(ecl.lat)
        # Make eclLon go between -180 and 180
        #if eclLon > 180:
        #    eclLon = eclLon - 360.
        # Make eclLon go between -180 and 180, changing where zero is
        #eclLon = eclLon - 180.

        sun = ephem.Sun(lsstObs)
        moon = ephem.Moon(lsstObs)
        moon_sun_sep = np.degrees(ephem.separation(moon,sun) )
        moon_target_sep = np.degrees(ephem.separation((moon.az,moon.alt ), (np.radians(targetAz),np.radians(targetAlt) )))
        esoobj.values['SKYMODEL.MOON.SUN.SEP'] = moon_sun_sep
        esoobj.values['SKYMODEL.MOON.TARGET.SEP'] =  moon_target_sep
        esoobj.values['SKYMODEL.MOON.ALT'] = np.degrees(moon.alt)
        esoobj.values["SKYMODEL.TIME"] = nightTime
        # Turn off zodiacal light until I understand the coord conversion
        esoobj.values["SKYMODEL.INCL.ZODIACAL"] = 'N'
    #    esoobj.values['SKYMODEL.ECL.LON'] = eclLon
    #    esoobj.values['SKYMODEL.ECL.LAT'] = eclLat
        sbs.append(esoobj.query(esoobj.values))
        sunAlts.append(sun.alt)
        moonAlts.append(moon.alt)
        moonSeps.append(moon_target_sep)

    esoobj.saveValues()

    cannonBands = {'R':'R', 'G':'V','B':'B'}
    colorcodes = ['r','g','b']

    # Set up figure for
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    ax1.set_xlabel('DJD')
    ax1.set_ylabel('mag/sq arcsec')
    #ax1.set_title('UT '+ut)

    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111)
    ax2.set_xlabel('sun alt (degrees)')
    ax2.set_ylabel('mag/sq arcsec')
    #ax2.set_title('UT '+ut)

    fig3 = plt.figure()
    ax3 = fig3.add_subplot(111)
    ax3.plot(nightData['time'], moonSeps)
    ax3.set_xlabel('DJD')
    #ax3.set_ylabel('Moon-target sep (degrees)')

    for key,colorc in zip(cannonBands.keys(), colorcodes):
        band = key
        vals = []
        for dic in sbs:
            vals.append(dic[cannonBands[key]])
        vals = np.array(vals)
        # Let's find the zeropoint from the middle of the night
        good = np.where( (nightData['time'] > night1_3) & (nightData['time'] < night2_3) )[0]
        zp = np.median(-2.5*np.log10(nightData[band][good]) - vals[good])


        ax1.plot(nightData['time'], -2.5*np.log10(nightData[band])-zp, colorc+'o')
        ax1.plot(nightData['time'], vals, 'k')

        ax2.plot(np.degrees(sunAlts), -2.5*np.log10(nightData[band])-zp, colorc)
        ax2.plot(np.degrees(sunAlts), vals, colorc+'--')

    # Flip the y axis limits
    ax1.set_ylim( ax1.get_ylim()[::-1] )
    ax2.set_ylim( ax2.get_ylim()[::-1] )

    plotName = os.path.join(params["esocomparepath"],"SB.png")
    fig1.savefig(plotName)
    plotName = os.path.join(params["esocomparepath"],"sunAlt.png")
    fig2.savefig(plotName)
    plotName = os.path.join(params["esocomparepath"],"moonSeps.png")
    fig3.savefig(plotName)

def healpix(params):

    filters = params["filters"]
    colors = params["colors"]
    data = {}

    nside = 512
    print "Pixel area: %.4f square degrees" % hp.nside2pixarea(nside, degrees=True)

    for filter in filters:

        baseplotDir = os.path.join(params["healpixpath"],filter)
        if not os.path.isdir(baseplotDir):
            os.mkdir(baseplotDir)

        outpath = os.path.join(params["dirpathname"],filter)
        outpath = os.path.join(params["fitspath"],filter)

        files = glob.glob(os.path.join(outpath,"*.long.%s.fits"%filter))
        files = sorted(files)

        magsDir = os.path.join(params["fisheyeplotpath"],filter)
        filename = os.path.join(magsDir,"mags.txt")
        mags_data = np.loadtxt(filename)

        ra_min = np.min(mags_data[:,0]) * (2*np.pi/360.0)
        ra_max = np.max(mags_data[:,0]) * (2*np.pi/360.0)
        dec_min = np.min(mags_data[:,1]) * (2*np.pi/360.0)
        dec_max = np.max(mags_data[:,1]) * (2*np.pi/360.0)

        phi_min = ra_min
        phi_max = ra_max
        theta_min = 0.5*np.pi - dec_max
        theta_max = 0.5*np.pi - dec_min

        runNumbers = []
        mjdobss = []
        skybrightness = []

        for file in files:

            folderName = file.split('.')[-4]
            runNumber = int(folderName)

            if runNumber > params["maxframes"]:
                continue

            hdulist = astropy.io.fits.open(file)

            npix = hp.nside2npix(nside)
            theta, phi = hp.pix2ang(nside, np.arange(npix))

            plotDir = os.path.join(baseplotDir,folderName)

            #if os.path.isdir(plotDir):
            #    continue

            if not os.path.isdir(plotDir):
                os.mkdir(plotDir)

            data_flat = np.ndarray.flatten(hdulist[0].data)
            image_array = hdulist[0].data

            image_array = np.flipud(image_array)
            #theta = np.linspace(0,np.pi, num=image_array.shape[0])[:, None]
            #theta = np.linspace((50.0/90.0)*np.pi/2,np.pi, num=image_array.shape[0])[:, None]
            #phi = np.linspace(-np.pi, np.pi, num=image_array.shape[1])
            theta = np.linspace(theta_min,theta_max, num=image_array.shape[0])[:, None]
            phi = np.linspace(phi_min, phi_max, num=image_array.shape[1])
            pix = hp.ang2pix(nside, theta, phi)
            healpix_map = np.zeros(hp.nside2npix(nside), dtype=np.double)
            healpix_map[pix] = image_array
 
            plotName = os.path.join(plotDir,'mollview.png')
            hp.mollview(healpix_map, min=2048, max=2100, title=folderName,xsize=2000)
            plt.show()
            plt.savefig(plotName,dpi=200)
            plt.close('all')

            plotName = os.path.join(plotDir,'gnomview.png')
            hp.gnomview(healpix_map, min=2048, max=2100, title=folderName,xsize=2000)
            plt.show()
            plt.savefig(plotName,dpi=200)
            plt.close('all')

            hpfits = os.path.join(plotDir,'hpfits.fits')
            hp.write_map(hpfits, healpix_map)

            #print plotName
            #print stop

        moviedir = os.path.join(baseplotDir,"movie")
        if not os.path.isdir(moviedir):
            os.mkdir(moviedir)

        folders = sorted(glob.glob(os.path.join(baseplotDir,'*')))
        n=1
        for folder in folders:
            n = n + 1
            file = os.path.join(folder,"mollview.png")
            filename = os.path.join(moviedir,"fishy-%04d.png"%n)
            cp_command = "cp %s %s"%(file,filename)
            os.system(cp_command)

        moviefiles = os.path.join(moviedir,"fishy-%04d.png")
        filename = os.path.join(moviedir,"mollview.mpg")
        ffmpeg_command = 'ffmpeg -an -y -r 20 -i %s -b:v %s %s'%(moviefiles,'5000k',filename)
        os.system(ffmpeg_command)
        filename = os.path.join(moviedir,"mollview.gif")
        ffmpeg_command = 'ffmpeg -an -y -r 20 -i %s -b:v %s %s'%(moviefiles,'5000k',filename)
        os.system(ffmpeg_command)

        rm_command = "rm %s/*.png"%(moviedir)
        os.system(rm_command)

        n=1
        for folder in folders:
            n = n + 1
            file = os.path.join(folder,"gnomview.png")
            filename = os.path.join(moviedir,"fishy-%04d.png"%n)
            cp_command = "cp %s %s"%(file,filename)
            os.system(cp_command)

        moviefiles = os.path.join(moviedir,"fishy-%04d.png")
        filename = os.path.join(moviedir,"gnomview.mpg")
        ffmpeg_command = 'ffmpeg -an -y -r 20 -i %s -b:v %s %s'%(moviefiles,'5000k',filename)
        os.system(ffmpeg_command)
        filename = os.path.join(moviedir,"gnomview.gif")
        ffmpeg_command = 'ffmpeg -an -y -r 20 -i %s -b:v %s %s'%(moviefiles,'5000k',filename)
        os.system(ffmpeg_command)

        rm_command = "rm %s/*.png"%(moviedir)
        os.system(rm_command)

def xy2sky(params):

    filters = params["filters"]
    colors = params["colors"]

    for filter in filters:

        baseplotDir = os.path.join(params["xy2skyplotpath"],filter)
        if not os.path.isdir(baseplotDir):
            os.mkdir(baseplotDir)

        #outpath = os.path.join(params["dirpathname"],filter)
        outpath = os.path.join(params["fitspath"],filter)
        files = glob.glob(os.path.join(outpath,"*.long.%s.fish"%filter))
        files = sorted(files)

        for file in files:

            folderName = file.split('.')[-4]
            runNumber = int(folderName)
            if runNumber > params["maxframes"]:
                continue

            fileNameSplit = file.split("/")[-1].replace("fish","fits")
            fitsfile = os.path.join(outpath,fileNameSplit)

            hdulist = astropy.io.fits.open(fitsfile)
            (a,b) = hdulist[0].data.shape

            LAT=-30.67833
            LNG=-70.73669

            lines = [line.strip() for line in open(file)]
            lineSplit = [x for x in lines[0].split(" ") if x != '']
            HA = float(lineSplit[2])
            DEC = float(lineSplit[4])
            AZ = float(lineSplit[6])
            SCALE = float(lineSplit[8])
            QUAD = float(lineSplit[10])
            CUBE = float(lineSplit[12])
            CX = float(lineSplit[14])
            CY = float(lineSplit[16])

            fitshdr_command = "fitshdr %s"%(fitsfile)
            p = subprocess.Popen(fitshdr_command.split(" "),stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            output, errors = p.communicate()

            xi = np.linspace(0,2896,2897)
            yi = np.linspace(0,1934,1935)

            xx, yy = np.meshgrid(xi, yi)
            xx = xx.flatten()
            yy = yy.flatten()
            zz = hdulist[0].data.flatten()

            ras = []
            decs = []

            for line in output.split("\n"):
                lineSplit = line.split(" ")
                lineSplit = [x for x in lineSplit if x != '']
                if len(lineSplit) == 0:
                    continue

                if lineSplit[0] == "MJD-OBS":
                    utcobs=float(lineSplit[2])
                    gmst = np.mod(280.460618+360.985647366*(utcobs-51544.5),360.0)
                    LST = np.mod(gmst+LNG,360.0)

                    for x in xi:
                        xxyyfile = file.replace("fish","xxyy")
                        f = open(xxyyfile,'w')
                        for i in xrange(len(yi)):
                            f.write('%.5f %.5f\n'%(x,yi[i]))
                        f.close()

                        fisheye_command = "awk '{print $1,$2}' %s | fisheye -ha %.5f -dec %.5f -az %.5f -scale %.5f -quad %.5f -cube %.5f -cx %.5f -cy %.5f -lat %.5f -xy2sky"%(xxyyfile,HA,DEC,AZ,SCALE,QUAD,CUBE,CX,CY,LAT)

                        output = os.popen(fisheye_command).readlines()
                        for line in output:
                            lineStrip = line.strip()
                            lineSplit = lineStrip.split(" ")
                            lineSplit = [x for x in lineSplit if x != '']
                     
                            ras.append(float(lineSplit[2]))
                            decs.append(float(lineSplit[3]))

            plotDir = os.path.join(baseplotDir,folderName)

            if not os.path.isdir(plotDir):
                os.mkdir(plotDir)

            xi = np.linspace(-180,180,3601)
            yi = np.linspace(-90,90,1801)
            zi = griddata(ras,decs,zz,xi,yi,interp='nn')
            zi = np.array(zi)

            zi[np.isnan(zi)] = 0

            z_min = np.min(zi)
            z_max = np.max(zi)

            #contour_levels = np.arange(-2.0, 2.0, 0.1)
            plotName = os.path.join(plotDir,'pcolor_ra_dec.png')
            plt.pcolor(xi, yi, zi, cmap='RdBu', vmin=z_min, vmax=z_max)
            # plot data points.
            plt.xlim([-180,180])
            plt.ylim([-90,90])
            plt.title(runNumber)
            cb = plt.colorbar()
            cb.set_label('vals')
            #cb.set_label('(minst - m) - m0')
            plt.show()
            plt.savefig(plotName,dpi=200)
            plt.close('all')

def combineclouds(params):

    filters = params["filters"]
    colors = params["colors"]

    for filter in filters:
        #outpath = os.path.join(params["dirpathname"],filter)
        outpath = os.path.join(params["fitspath"],filter)

        baseplotDir = os.path.join(params["cloudsplotpath"],filter)
        if not os.path.isdir(baseplotDir):
            os.mkdir(baseplotDir)

        folders = glob.glob(os.path.join(baseplotDir,"*"))
        folders = sorted(folders)        

        runNumbers = []
        medianVals = []
        stdVals = []
        magPerc = []
        mjds = []
        nstars = []

        magThresh = -1        
        bins = np.arange(-4,1.1,0.1)
        histMatrix = []

        fisheyeDir = os.path.join(params["fisheyeplotpath"],filter)
        fileName = os.path.join(fisheyeDir,'nstars.txt')
        nstarsData = np.loadtxt(fileName)

        for folder in folders:

            folderSplit = folder.split('/')
            folderName = folderSplit[-1]

            try:
                runNumber = int(folderName)
            except:
                print "%s is not a number..."%folderName 
                continue

            #if not runNumber == 898:
            #   continue

            #if not runNumber == 700:
            #   continue

            if runNumber > params["maxframes"]:
                continue

            plotDir = os.path.join(baseplotDir,folderName)

            interpMethod = 0
            if interpMethod:

               fitsfile = os.path.join(plotDir,'clouds.fits')

               if not os.path.isfile(fitsfile):
                   continue

               hdulist = astropy.io.fits.open(fitsfile)

               scidata = hdulist[0].data
               im = scidata

               hist_vals = im.ravel()
               hist_vals = hist_vals[~np.isnan(hist_vals)]

               hist_vals_median = np.median(hist_vals)
               hist_vals_std = np.std(hist_vals)

               indexes = np.where(hist_vals < magThresh)[0]
               perc = 100.0 * len(indexes) / len(hist_vals)

               magPerc.append(perc)

               runNumbers.append(runNumber)
               medianVals.append(hist_vals_median)
               stdVals.append(hist_vals_std)

               index = np.argmin(np.absolute(nstarsData[:,0] - runNumber))
               mjds.append(nstarsData[index,1])
               nstars.append(nstarsData[index,2])

               hist, bins = np.histogram(hist_vals, bins=bins)
               hist = 100 * hist / hist.sum()

            else:

               trianglesFile = os.path.join(plotDir,'triangles.txt')

               if not os.path.isfile(trianglesFile):
                   continue

               try:
                   data_out = np.loadtxt(trianglesFile)
               except:
                   continue
               hist_vals = data_out[:,0]
               weights = data_out[:,-1]                    

               hist_vals_median = np.median(hist_vals)
               hist_vals_std = np.std(hist_vals)

               indexes = np.where(hist_vals < magThresh)[0]
               perc = 100.0 * len(indexes) / len(hist_vals)

               magPerc.append(perc)

               runNumbers.append(runNumber)
               medianVals.append(hist_vals_median)
               stdVals.append(hist_vals_std)

               index = np.argmin(np.absolute(nstarsData[:,0] - runNumber))
               mjds.append(nstarsData[index,1])
               nstars.append(nstarsData[index,2])

               hist, bins = np.histogram(hist_vals, bins=bins, weights=weights)
               hist = 100 * hist / hist.sum()

            if len(histMatrix) == 0:
                histMatrix = hist
            else:
                histMatrix = np.vstack([histMatrix,hist])

        tt = np.array(mjds)
        if len(tt) == 0:
            print "No available data for %s"%filter
            continue

        plotDir = os.path.join(baseplotDir,"combine")
        if not os.path.isdir(plotDir):
            os.mkdir(plotDir)

        histFile = os.path.join(plotDir,'hist.txt')
        f = open(histFile,'w')
        for i in xrange(len(tt)):
            f.write('%.10f'%tt[i])
            for j in xrange(len(bins)-1):
                f.write(' %.10f '%histMatrix[i,j])
            f.write('\n')
        f.close()

        percFile = os.path.join(plotDir,'perc.txt')
        f = open(percFile,'w')
        for i in xrange(len(tt)):
            f.write('%.10f %.10f\n'%(tt[i],magPerc[i]))
        f.close()

        tt = tt - tt[0]
        tt = tt * 24

        plotName = os.path.join(plotDir,'timeseries.png')
        plt.errorbar(tt, medianVals, yerr=stdVals)
        plt.xlabel('Hours since twlight')
        plt.ylabel('delta magnitude')
        plt.ylim([0,100])
        plt.show()
        plotName = os.path.join(plotDir,'timeseries.png')
        plt.savefig(plotName,dpi=200)
        plotName = os.path.join(plotDir,'timeseries.eps')
        plt.savefig(plotName,dpi=200)
        plt.close('all')

        plt.plot(tt, magPerc)
        plt.xlabel('Hours since twlight')
        plt.ylabel('percentage greater than 1 magnitude')
        plt.show()
        plotName = os.path.join(plotDir,'perc.png')
        plt.savefig(plotName,dpi=200)
        plotName = os.path.join(plotDir,'perc.eps')
        plt.savefig(plotName,dpi=200)
        plt.close('all')

        data = histMatrix
        mags = bins[:-1]

        X,Y = np.meshgrid(tt,mags)
        ax = plt.subplot(111)
        im = plt.pcolor(X,Y,data.T, cmap=plt.cm.jet,vmin=0,vmax=10)
        plt.xlabel("Hours since twlight")
        plt.ylabel("delta magnitude")
        cbar=plt.colorbar()
        cbar.set_label(r'PDF')
        plt.xlim([np.min(tt),np.max(tt)])
        plt.ylim([np.min(mags),np.max(mags)])
        plt.grid()
        plt.show()
        plotName = os.path.join(plotDir,'tf.png')
        plt.savefig(plotName,dpi=200)
        plotName = os.path.join(plotDir,'tf.eps')
        plt.savefig(plotName,dpi=200)
        plt.close('all')

        low = data.min()/2
        high = data.max() * 2
        nbins = 50
        bins = np.linspace(low, high, num=nbins+1)        
        bins,specvar = utils.spectral_histogram(data,bins=bins)
        bins = bins[:-1]
        percentile = 10
        spectral_variation_10per = utils.spectral_percentiles(specvar,bins,percentile)
        percentile = 50
        spectral_variation_50per = utils.spectral_percentiles(specvar,bins,percentile)
        percentile = 90
        spectral_variation_90per = utils.spectral_percentiles(specvar,bins,percentile)

        plotName = os.path.join(plotDir,'specvar.png')
        X,Y = np.meshgrid(mags,bins)
        ax = plt.subplot(111)
        #im = plt.pcolor(X,Y,np.transpose(spectral_variation_norm), cmap=plt.cm.jet)

        masked_array = np.ma.array(specvar, mask=np.isnan(specvar))
        cmap = plt.cm.jet
        cmap.set_bad('w',1.)

        #im = plt.pcolor(X,Y,specvar, cmap=plt.cm.jet,vmin=0,vmax=100)
        im = plt.pcolor(X,Y,masked_array, cmap=cmap,vmin=0,vmax=50)
        #ax.set_xscale('log')
        #ax.set_yscale('log')
        plt.plot(mags,spectral_variation_10per,'w',label='10')
        plt.plot(mags,spectral_variation_50per,'w',label='50')
        plt.plot(mags,spectral_variation_90per,'w',label='90')
        #plt.loglog(fl,low,'k-.',linewidth=3)
        #plt.loglog(fh,high,'k-.',label='LNM/HNM',linewidth=3)
        plt.xlim([-3,2])
        plt.ylim([0,50]) 
        #plt.ylim([1e-10,1e-5])
        plt.xlabel("Delta magnitude")
        plt.ylabel("PDF")
        cbar=plt.colorbar()
        cbar.set_label(r'Percentile')
        plt.grid()
        plt.show()
        plt.savefig(plotName,dpi=200)
        plt.close('all')

def clouds(params):

    filters = params["filters"]
    colors = params["colors"]

    #filters = [filters[-1]]

    for filter in filters:
        #outpath = os.path.join(params["dirpathname"],filter)
        outpath = os.path.join(params["fitspath"],filter)

        baseplotDir = os.path.join(params["cloudsplotpath"],filter)
        if not os.path.isdir(baseplotDir):
            os.mkdir(baseplotDir)

        magsDir = os.path.join(params["fisheyeplotpath"],filter)
        filename = os.path.join(magsDir,"mags.txt")

        #filename = os.path.join(params["catalogpath"],"%s.txt"%filter)
        mags_data = np.loadtxt(filename)

        ra_min = np.min(mags_data[:,0]) * (2*np.pi/360.0)
        ra_max = np.max(mags_data[:,0]) * (2*np.pi/360.0)
        dec_min = np.min(mags_data[:,1]) * (2*np.pi/360.0)
        dec_max = np.max(mags_data[:,1]) * (2*np.pi/360.0)

        phi_min = ra_min
        phi_max = ra_max
        theta_min = 0.5*np.pi - dec_max
        theta_max = 0.5*np.pi - dec_min

        plotName = os.path.join(baseplotDir,'mags_hist.png')
        bins = np.arange(10,20,0.1)
        hist, bins = np.histogram(mags_data[:,2], bins=bins)
        width = 0.7 * (bins[1] - bins[0])
        center = (bins[:-1] + bins[1:]) / 2
        plt.bar(center, hist, align='center', width=width)
        plt.show()
        plt.savefig(plotName,dpi=200)
        plt.close('all')

        plotName = os.path.join(baseplotDir,'mags_slope_hist.png')
        bins = np.arange(-1,1,0.1)
        hist, bins = np.histogram(mags_data[:,3], bins=bins)
        width = 0.7 * (bins[1] - bins[0])
        center = (bins[:-1] + bins[1:]) / 2
        plt.bar(center, hist, align='center', width=width)
        plt.show()
        plt.savefig(plotName,dpi=200)
        plt.close('all')

        mags_data_slope_median = np.median(mags_data[:,3])
        print mags_data_slope_median

        files = glob.glob(os.path.join(outpath,"*.long.%s.fish"%filter))
        files = sorted(files)

        zi_ave = []

        #for file in []:
        for file in files:

            folderName = file.split('.')[-4]
            runNumber = int(folderName)

            #if not runNumber == 898:
            #   continue

            #if not runNumber == 110:
            #   continue

            if runNumber > params["maxframes"]:
                continue

            fitsfile = file.replace("fish","fits")
            fileprefix = file.replace(".fits","")
            outputfile = file.replace("fish","mch")

            if not os.path.isfile(outputfile):
                continue

            data_out = np.loadtxt(outputfile,comments='RA')
            if len(data_out) == 0:
                continue

            #if not runNumber == 400:
            #    continue 

            data_all = data_out

            indexes = np.where(~np.isnan(data_all[:,17]))[0]
            data_all = data_all[indexes,:]
            magnitudes = data_all[:,10] - data_all[:,17]

            zenith = (90 - data_all[:,15]) * 2 * np.pi / 360.0
            secant = (1.002432 * np.cos(zenith)**2 + 0.148386 * np.cos(zenith) + 0.0096467) / (np.cos(zenith)**3 + 0.149864 * np.cos(zenith)**2 + 0.0102963 * np.cos(zenith) + 0.00303978)

            plotDir = os.path.join(baseplotDir,folderName)

            #if os.path.isdir(plotDir):
            #    continue

            if not os.path.isdir(plotDir):
                os.mkdir(plotDir)

            plotName = os.path.join(plotDir,'mag_x_y.png')
            plt.scatter(data_all[:,3],data_all[:,4],s=10,c=magnitudes,vmin=12,vmax=15,edgecolor='none')
            plt.xlim([0,2897])
            plt.ylim([0,1935])
            plt.xlabel('x')
            plt.ylabel('y')
            plt.title(runNumber)
            cb = plt.colorbar()
            #cb.set_label('minst')
            cb.set_label('minst - m')
            plt.show()
            plt.savefig(plotName,dpi=200)
            plt.close('all')

            LAT=-30.67833
            LNG=-70.73669
            #HA=0.015
            #DEC=-27.092
            #AZ=95.082
            #SCALE=0.04705

            lines = [line.strip() for line in open(file)]
            lineSplit = [x for x in lines[0].split(" ") if x != '']
            HA = float(lineSplit[2])
            DEC = float(lineSplit[4])
            AZ = float(lineSplit[6])
            SCALE = float(lineSplit[8])
            QUAD = float(lineSplit[10])
            CUBE = float(lineSplit[12])
            CX = float(lineSplit[14])
            CY = float(lineSplit[16])

            fitshdr_command = "fitshdr %s"%(fitsfile)
            p = subprocess.Popen(fitshdr_command.split(" "),stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            output, errors = p.communicate()
            for line in output.split("\n"):
                lineSplit = line.split(" ")
                lineSplit = [x for x in lineSplit if x != '']
                if len(lineSplit) == 0:
                    continue

                if lineSplit[0] == "MJD-OBS":
                    utcobs=float(lineSplit[2])
                    gmst = np.mod(280.460618+360.985647366*(utcobs-51544.5),360.0)
                    LST = np.mod(gmst+LNG,360.0)

            mags = []
            magnitudes = []
            xs = []
            ys = []

            for i in xrange(len(mags_data[:,0])):

                ra = mags_data[i,0]
                dec = mags_data[i,1]

                index1 = np.where(ra == data_all[:,0])[0]
                index2 = np.where(dec == data_all[:,1])[0]

                index = index1

                zenith = (90 - data_all[index,15]) * 2 * np.pi / 360.0
                secant = 1/np.cos(zenith)
                ha = LST - ra


                fisheye_command = "echo '%.5f %.5f' | fisheye -ha %.5f -dec %.5f -az %.5f -scale %.5f -quad %.5f -cube %.5f -cx %.5f -cy %.5f -lat %.5f -sky2xy"%(ha,dec,HA,DEC,AZ,SCALE,QUAD,CUBE,CX,CY,LAT)
                output = os.popen(fisheye_command).readlines()[0].strip()
                outputSplit = output.split(" ")
                outputSplit = [x for x in outputSplit if x != '']

                x = float(outputSplit[0])
                y = float(outputSplit[1])

                if np.sqrt(((2897.0/2) - x)**2 + ((1935.0/2) - y)**2) > 1100:
                    continue

                xs.append(x)
                ys.append(y)

                if len(data_all[index,2]) == 0:
                    #mags.append(magnitudes[i])
                    mags.append(-1000)
                    magnitudes.append(-1000)
                else:

                    mag_estimate = mags_data[index,2][0] + secant*mags_data[index,3][0]
                    #mag_estimate = mags_data[index,2][0] + secant*mags_data_slope_median
                    mags.append(mag_estimate[0])
                    magnitude = data_all[index,10] - data_all[index,17]
                    magnitudes.append(magnitude[0])

            xs = np.array(xs)
            ys = np.array(ys)

            #mags = []
            #magnitudes = data_all[:,10] - data_all[:,17]
            #mags = data_all[:,10] - data_all[:,17] + 1
            #mags = mags.T
            #xs = data_all[:,3]
            #ys = data_all[:,4]
            #xs = np.array(xs)
            #ys = np.array(ys)

            for i in xrange(len(data_all[:,0])):
                continue
                ra = data_all[i,0]
                dec = data_all[i,1]
                zenith = (90 - data_all[i,15]) * 2 * np.pi / 360.0
                secant = 1/np.cos(zenith)

                ha = LST - ra
                index1 = np.where(ra == mags_data[:,0])[0]
                index2 = np.where(dec == mags_data[:,1])[0]

                index = index1

                if len(mags_data[index,2]) == 0:
                    #mags.append(magnitudes[i])
                    mags.append(-1000)
                else:
 
                    mag_estimate = mags_data[index,2][0] + secant*mags_data[index,3][0]
                    #mag_estimate = mags_data[index,2][0] + secant*mags_data_slope_median
                    mags.append(mag_estimate)

            #data_all = data_all_new
            #magnitudes = data_all[:,10] - data_all[:,17]
            #print mags
            magnitudes = np.array(magnitudes)
            mags = np.array(mags).T
            diffs = magnitudes - mags

            diffs[mags==-1000] = -4

            if np.sum(diffs) == 0.0:
                continue

            if len(diffs) < 3:
                continue

            colorRange = 2.0

            plotName = os.path.join(plotDir,'diff_x_y.png')
            plt.scatter(xs,ys,s=10,c=diffs,vmin=-colorRange,vmax=colorRange,edgecolor='none')
            plt.xlim([0,2897])
            plt.ylim([0,1935])
            plt.xlabel('x')
            plt.ylabel('y')
            plt.title(runNumber)
            cb = plt.colorbar()
            #cb.set_label('minst')
            cb.set_label('(minst - m) - m0')
            plt.show()
            plt.savefig(plotName,dpi=200)
            plt.close('all')

            plotName = os.path.join(plotDir,'sky_x_y.png')
            plt.scatter(data_all[:,3],data_all[:,4],s=10,c=data_all[:,19],vmin=0,vmax=10,edgecolor='none')
            plt.xlim([0,2897])
            plt.ylim([0,1935])
            plt.xlabel('x')
            plt.ylabel('y')
            plt.title(runNumber)
            cb = plt.colorbar()
            #cb.set_label('minst')
            cb.set_label('sky')
            plt.show()
            plt.savefig(plotName,dpi=200)
            plt.close('all')

            xi = np.linspace(0,2897,1000)
            yi = np.linspace(0,1935,1000)

            #xi = np.linspace(0,3000,1000)
            #yi = np.linspace(-200,2000,1000)

            zi = griddata(xs.T,ys.T,diffs,xi,yi,interp='nn')
            zi = np.array(zi)
            zi_copy = zi.copy()
            zi_copy[np.isnan(zi_copy)] = 0

            if len(zi_ave) == 0:
                zi_ave = zi_copy
            else:
                zi_ave = zi_ave + zi_copy

            hist_vals = zi.ravel()
            hist_vals = hist_vals[~np.isnan(hist_vals)]
 
            bins = np.arange(-3,3,0.1)
            plotName = os.path.join(plotDir,'hist.png')
            hist, bins = np.histogram(hist_vals, bins=bins)

            width = 0.7 * (bins[1] - bins[0])
            center = (bins[:-1] + bins[1:]) / 2
            plt.bar(center, hist, align='center', width=width)
            plt.show()
            plt.savefig(plotName,dpi=200)
            plt.close('all')

            zi_ra = griddata(data_all[:,3],data_all[:,4],data_all[:,0],xi,yi,interp='nn')
            zi_ra = np.array(zi_ra)
            contour_levels_ra = np.arange(0.0, 360.0, 60)

            zi_dec = griddata(data_all[:,3],data_all[:,4],data_all[:,1],xi,yi,interp='nn')
            zi_dec = np.array(zi_dec)
            contour_levels_dec = np.arange(-90, 90.0, 10)

            contour_levels = np.arange(-4.0, 1.0, 0.5)
            plotName = os.path.join(plotDir,'contour_x_y.png')
            plt.scatter(xs,ys,s=10,c=diffs,zorder=10,vmax=colorRange, vmin=-colorRange)
            #plt.scatter(data_all[:,3],data_all[:,4],s=10,c=diffs,zorder=10,vmax=2, vmin=-2,edgecolor='none')
            plt.contour(xi,yi,zi,contour_levels,linewidths=0.5,colors='k')
            plt.contourf(xi,yi,zi,contour_levels,cmap=plt.cm.rainbow,vmax=colorRange, vmin=-colorRange)

            plt.xlabel('PIXEL')
            plt.ylabel('PIXEL')
            cb = plt.colorbar()
            #cb.set_label('minst')
            cb.set_label('Change in magnitude')            

            CS_ra = plt.contour(xi,yi,zi_ra,contour_levels_ra,linewidths=0.5,colors='k')
            plt.clabel(CS_ra, fontsize=9, inline=1)
            CS_dec = plt.contour(xi,yi,zi_dec,contour_levels_dec,linewidths=0.5,colors='k')
            plt.clabel(CS_dec, fontsize=9, inline=1)

            #cb = plt.colorbar() # draw colorbar
            # plot data points.
            plt.xlim([0,2897])
            plt.ylim([0,1935])
            plt.title(runNumber)
            plt.xlim([0,3000])
            plt.ylim([0,3000])
            plt.show()
            plt.savefig(plotName,dpi=200)
            plt.close('all')

            points = np.vstack([xs,ys]).T
            tri = Delaunay(points)
            plotName = os.path.join(plotDir,'Delaunay.png')

            theta = np.linspace(0,2*np.pi,100)
            cradiusx = 1100.0
            cradiusy = 1100.0
            xradius = (2897.0/2) + cradiusx*np.cos(theta)
            yradius = (1935.0/2) - cradiusy*np.sin(theta)

            trianglesFile = os.path.join(plotDir,'triangles.txt')
            f = open(trianglesFile,'w')

            xtris = []
            ytris = []
            trivals =[]
            areas = []
            for indexes in tri.simplices:
                trivals.append(np.max(diffs[indexes]))
                xtris.append(xs[indexes])
                ytris.append(ys[indexes])

                a = np.array([xs[indexes][0],ys[indexes][0]])
                b = np.array([xs[indexes][1],ys[indexes][1]])
                c = np.array([xs[indexes][2],ys[indexes][2]])

                ab=a-b
                ac=a-c
                area=np.abs(ab[0]*ac[1] - ab[1]*ac[0])/2
                areas.append(area)

                f.write('%.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f\n'%(\
                    np.max(diffs[indexes]),\
                    xs[indexes][0],xs[indexes][1],xs[indexes][2],\
                    ys[indexes][0],ys[indexes][1],ys[indexes][2],\
                    area))

            plt.triplot(points[:,0], points[:,1], tri.simplices.copy())
            plt.plot(points[:,0], points[:,1], 'o')

            plt.plot(xradius,yradius,'k--')

            plt.show()
            plt.savefig(plotName,dpi=200)
            plt.close('all')

            plotName = os.path.join(plotDir,'triangles.png')
            fig, ax = plt.subplots()
            range = np.linspace(0, 1, 100)
            colorsall = cm.rainbow(range)
            
            patches = []
            colors = []
            fig, ax = plt.subplots()
            plt.scatter(xs,ys,s=10,c=diffs,zorder=0,vmax=1, vmin=-4,alpha=1)
            for xtips,ytips,trival in zip(xtris,ytris,trivals):
                thisrange = (trival - -4.0)/6.0
                index = np.argmin(np.absolute(range-thisrange))
                color = colorsall[index]
                colors.append(color)
                plt.fill(xtips,ytips,
                facecolor=color,alpha=1, edgecolor='black')
                #polygon = Polygon(np.vstack((xtips,ytips)).T, True)
                #patches.append(polygon)
            cb = plt.colorbar()
            cb.set_label(r'Change in magnitude')

            plt.xlabel('PIXEL')
            plt.ylabel('PIXEL')
            plt.xlim([0,3000])
            plt.ylim([0,2000])
            plt.title(runNumber)
 
            plt.show()
            plt.savefig(plotName,dpi=200)
            plt.close('all')

            #print plotName
            #print stop

            nside = 512
            image_array = zi

            image_array = np.flipud(image_array)
            #theta = np.linspace(0,np.pi, num=image_array.shape[0])[:, None]
            #theta = np.linspace((50.0/90.0)*np.pi/2,np.pi, num=image_array.shape[0])[:, None]
            #phi = np.linspace(-np.pi, np.pi, num=image_array.shape[1])
            theta = np.linspace(theta_min,theta_max, num=image_array.shape[0])[:, None]
            phi = np.linspace(phi_min, phi_max, num=image_array.shape[1])
            pix = hp.ang2pix(nside, theta, phi)
            healpix_map = np.zeros(hp.nside2npix(nside), dtype=np.double)
            healpix_map[pix] = image_array

            plotName = os.path.join(plotDir,'mollview.png')
            hp.mollview(healpix_map, min=-colorRange, max=colorRange, title=folderName,xsize=2000)
            plt.show()
            plt.savefig(plotName,dpi=200)
            plt.close('all')

            plotName = os.path.join(plotDir,'gnomview.png')
            hp.gnomview(healpix_map, min=-colorRange, max=colorRange, title=folderName,xsize=2000)
            plt.show()
            plt.savefig(plotName,dpi=200)
            plt.close('all')

            hpfits = os.path.join(plotDir,'hpfits.fits')
            hp.write_map(hpfits, healpix_map)

            fitsfile = os.path.join(plotDir,'clouds.fits')
            hdulist = astropy.io.fits.PrimaryHDU(zi)
            hdulist.writeto(fitsfile,output_verify='warn',clobber=True)

            zi = griddata(data_all[:,3],data_all[:,4],data_all[:,19],xi,yi,interp='nn')
            zi = np.array(zi)

            contour_levels = np.arange(5, 10.0, 0.1)
            plotName = os.path.join(plotDir,'skycontour_x_y.png')
            plt.scatter(data_all[:,3],data_all[:,4],s=10,c=data_all[:,19],zorder=10,vmax=10, vmin=5)
            #plt.scatter(data_all[:,3],data_all[:,4],s=10,c=diffs,zorder=10,vmax=2, vmin=-2,edgecolor='none')
            plt.contour(xi,yi,zi,contour_levels,linewidths=0.5,colors='k')
            plt.contourf(xi,yi,zi,contour_levels,cmap=plt.cm.rainbow)
            #cb = plt.colorbar() # draw colorbar
            # plot data points.
            plt.xlim([0,2897])
            plt.ylim([0,1935])
            plt.xlim([0,3000])
            plt.ylim([0,3000])
            plt.title(runNumber)
            cb = plt.colorbar()
            #cb.set_label('minst')
            cb.set_label('Sky brightness')
            plt.show()
            plt.savefig(plotName,dpi=200)
            plt.close('all')

            zi = griddata(data_all[:,3],data_all[:,4],data_all[:,0],xi,yi,interp='nn')
            zi = np.array(zi)

            contour_levels = np.arange(0, 360.0, 60)
            plotName = os.path.join(plotDir,'ra.png')
            plt.scatter(data_all[:,3],data_all[:,4],s=10,c=data_all[:,0],zorder=10,vmax=360, vmin=0)
            #plt.scatter(data_all[:,3],data_all[:,4],s=10,c=diffs,zorder=10,vmax=2, vmin=-2,edgecolor='none')
            plt.contour(xi,yi,zi,contour_levels,linewidths=0.5,colors='k')
            #plt.contourf(xi,yi,zi,contour_levels,cmap=plt.cm.rainbow)
            #for contour_level in contour_levels:
            #    print contour_level
            #    plt.contour(xi,yi,zi,[contour_level],linewidths=0.5,colors='k')
            #cb = plt.colorbar() # draw colorbar
            # plot data points.
            plt.xlim([0,2897])
            plt.ylim([0,1935])
            plt.title(runNumber)
            #cb = plt.colorbar()
            #cb.set_label('minst')
            #cb.set_label('(minst - m) - m0')
            plt.show()
            plt.savefig(plotName,dpi=200)
            plt.close('all')

            zi = griddata(data_all[:,3],data_all[:,4],data_all[:,1],xi,yi,interp='nn')
            zi = np.array(zi)

            contour_levels = np.arange(-90, 90.0, 5)
            plotName = os.path.join(plotDir,'dec.png')
            #plt.scatter(data_all[:,3],data_all[:,4],s=10,c=data_all[:,19],zorder=10,vmax=10, vmin=0)
            #plt.scatter(data_all[:,3],data_all[:,4],s=10,c=diffs,zorder=10,vmax=2, vmin=-2,edgecolor='none')
            plt.contour(xi,yi,zi,contour_levels,linewidths=0.5,colors='k')
            plt.contourf(xi,yi,zi,contour_levels,cmap=plt.cm.rainbow)
            #cb = plt.colorbar() # draw colorbar
            # plot data points.
            plt.xlim([0,2897])
            plt.ylim([0,1935])
            plt.title(runNumber)
            cb = plt.colorbar()
            #cb.set_label('minst')
            #cb.set_label('(minst - m) - m0')
            plt.show()
            plt.savefig(plotName,dpi=200)
            plt.close('all')

            #p = cs.collections[0].get_paths()[0]
            #v = p.vertices
            #x = v[:,0]
            #y = v[:,1]

        zi_ave_norm = np.array(zi_ave) / float(len(files))

        plotDir = os.path.join(baseplotDir,"ave")
        if not os.path.isdir(plotDir):
            os.mkdir(plotDir)

        if True:
            contour_levels = np.arange(-1, 1, 0.1)
            plotName = os.path.join(plotDir,'contour_x_y.png')
            plt.contour(xi,yi,zi_ave_norm,contour_levels,linewidths=0.5,colors='k')
            plt.contourf(xi,yi,zi_ave_norm,contour_levels,cmap=plt.cm.rainbow)
            cb = plt.colorbar() # draw colorbar
            # plot data points.
            plt.xlim([0,2897])
            plt.ylim([0,1935])
            plt.title("average loss")
            #cb = plt.colorbar()
            #cb.set_label('minst')
            #cb.set_label('(minst - m) - m0')
            plt.show()
            plt.savefig(plotName,dpi=200)
            plt.close('all')

        moviedir = os.path.join(baseplotDir,"movie")
        if not os.path.isdir(moviedir):
            os.mkdir(moviedir)

        folders = sorted(glob.glob(os.path.join(baseplotDir,'*')))
        n=1
        for folder in folders:
            n = n + 1
            file = os.path.join(folder,"diff_x_y.png")
            filename = os.path.join(moviedir,"fishy-%04d.png"%n)
            cp_command = "cp %s %s"%(file,filename)
            os.system(cp_command)

        moviefiles = os.path.join(moviedir,"fishy-%04d.png")
        filename = os.path.join(moviedir,"diff.mpg")
        ffmpeg_command = 'ffmpeg -an -y -r 20 -i %s -b:v %s %s'%(moviefiles,'5000k',filename)
        os.system(ffmpeg_command)
        filename = os.path.join(moviedir,"diff.gif")
        ffmpeg_command = 'ffmpeg -an -y -r 20 -i %s -b:v %s %s'%(moviefiles,'5000k',filename)
        os.system(ffmpeg_command)

        rm_command = "rm %s/*.png"%(moviedir)
        os.system(rm_command)

        n=1
        for folder in folders:
            n = n + 1
            file = os.path.join(folder,"contour_x_y.png")
            filename = os.path.join(moviedir,"fishy-%04d.png"%n)
            cp_command = "cp %s %s"%(file,filename)
            os.system(cp_command)

        moviefiles = os.path.join(moviedir,"fishy-%04d.png")
        filename = os.path.join(moviedir,"contour.mpg")
        ffmpeg_command = 'ffmpeg -an -y -r 20 -i %s -b:v %s %s'%(moviefiles,'5000k',filename)
        os.system(ffmpeg_command)
        filename = os.path.join(moviedir,"contour.gif")
        ffmpeg_command = 'ffmpeg -an -y -r 20 -i %s -b:v %s %s'%(moviefiles,'5000k',filename)
        os.system(ffmpeg_command)

        rm_command = "rm %s/*.png"%(moviedir)
        os.system(rm_command)

        n=1
        for folder in folders:
            n = n + 1
            file = os.path.join(folder,"skycontour_x_y.png")
            filename = os.path.join(moviedir,"fishy-%04d.png"%n)
            cp_command = "cp %s %s"%(file,filename)
            os.system(cp_command)

        moviefiles = os.path.join(moviedir,"fishy-%04d.png")
        filename = os.path.join(moviedir,"skycontour.mpg")
        ffmpeg_command = 'ffmpeg -an -y -r 20 -i %s -b:v %s %s'%(moviefiles,'5000k',filename)
        os.system(ffmpeg_command)
        filename = os.path.join(moviedir,"skycontour.gif")
        ffmpeg_command = 'ffmpeg -an -y -r 20 -i %s -b:v %s %s'%(moviefiles,'5000k',filename)
        os.system(ffmpeg_command)

        rm_command = "rm %s/*.png"%(moviedir)
        os.system(rm_command)

        n=1
        for folder in folders:
            n = n + 1
            file = os.path.join(folder,"mollview.png")
            filename = os.path.join(moviedir,"fishy-%04d.png"%n)
            cp_command = "cp %s %s"%(file,filename)
            os.system(cp_command)

        moviefiles = os.path.join(moviedir,"fishy-%04d.png")
        filename = os.path.join(moviedir,"mollview.mpg")
        ffmpeg_command = 'ffmpeg -an -y -r 20 -i %s -b:v %s %s'%(moviefiles,'5000k',filename)
        os.system(ffmpeg_command)
        filename = os.path.join(moviedir,"mollview.gif")
        ffmpeg_command = 'ffmpeg -an -y -r 20 -i %s -b:v %s %s'%(moviefiles,'5000k',filename)
        os.system(ffmpeg_command)

        rm_command = "rm %s/*.png"%(moviedir)
        os.system(rm_command)

        n=1
        for folder in folders:
            n = n + 1
            file = os.path.join(folder,"triangles.png")
            filename = os.path.join(moviedir,"fishy-%04d.png"%n)
            cp_command = "cp %s %s"%(file,filename)
            os.system(cp_command)

        moviefiles = os.path.join(moviedir,"fishy-%04d.png")
        filename = os.path.join(moviedir,"triangles.mpg")
        ffmpeg_command = 'ffmpeg -an -y -r 20 -i %s -b:v %s %s'%(moviefiles,'5000k',filename)
        os.system(ffmpeg_command)
        filename = os.path.join(moviedir,"triangles.gif")
        ffmpeg_command = 'ffmpeg -an -y -r 20 -i %s -b:v %s %s'%(moviefiles,'5000k',filename)
        os.system(ffmpeg_command)

        rm_command = "rm %s/*.png"%(moviedir)
        os.system(rm_command)

def getVariableStars():

    variableStars = []
    variableStarsFile = "/lsst/home/coughlin/git-repo/fisheye/variable_stars.dat"
    lines = [line.strip() for line in open(variableStarsFile)]
    for line in lines:
        lineSplit = line.split("|")
        ra = lineSplit[4]
        raSplit = ra.split(" ")
        if raSplit[0] == "":
            continue
        raSplit = [float(x) for x in ra.split(" ")]
        ra = (raSplit[0] + raSplit[1]/60.0 + raSplit[2]/3600.0) * 360.0 / 24.0
        dec = lineSplit[5]
        decSplit = [float(x) for x in dec.split(" ")]
        dec = (decSplit[0] + decSplit[1]/60.0 + decSplit[2]/3600.0)

        variableStars.append((ra,dec))
    variableStars = np.array(variableStars)
 
    return variableStars

def m_func(x_4d, alpha, beta, gamma, m_0):
    m_corr = x_4d[:,0] - (alpha + beta*x_4d[:,1])*x_4d[:,2] - gamma * x_4d[:,3] + m_0
    #m_corr = x_4d[:,0] - (alpha + beta*x_4d[:,1])*x_4d[:,2] - gamma * 0 + m_0
    return m_corr

def combinenight(params):

    filters = params["filters"]
    colors = params["colors"]

    for filter in filters:
    #for filter in []:
        #outpath = os.path.join(params["dirpathname"],filter)
        outpath = os.path.join(params["fitspath"],filter)

        data_all = []
        data_all_tph = []

        files = glob.glob(os.path.join(outpath,"*.long.%s.fits"%filter))
        #files = glob.glob(os.path.join(outpath,"*.short.%s.fish"%filter))
        files = sorted(files)

        runNumbers = []

        filename = os.path.join(params["catalognightpath"],"%s/%s.txt"%(filter,params["dirname"]))
        fid = open(filename,'w')

        for file in files:

            runNumber = int(file.split('.')[-4])
            if runNumber > params["maxframes"]:
                continue

            fileprefix = file.replace(".fits","")
            outputfile = file.replace("fits","mch")
            tphfile = file.replace("fits","tph")

            if not os.path.isfile(outputfile):
                continue

            fitshdr_command = "fitshdr %s"%(file)
            p = subprocess.Popen(fitshdr_command.split(" "),stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            output, errors = p.communicate()

            for line in output.split("\n"):
                lineSplit = line.split(" ")
                lineSplit = [x for x in lineSplit if x != '']
                if len(lineSplit) == 0:
                    continue

                if lineSplit[0] == "BIAS":
                    bias=float(lineSplit[2])
                elif lineSplit[0] == "EXPTIME":
                    exptime=float(lineSplit[2])
                elif lineSplit[0] == "MJD-OBS":
                    mjdobs=float(lineSplit[2])
                elif lineSplit[0] == "NAXIS1":
                    nx=float(lineSplit[2])
                elif lineSplit[0] == "NAXIS2":
                    ny=float(lineSplit[2])

            data_out = np.loadtxt(outputfile,comments='RA')
            lines = [line.rstrip() for line in open(outputfile)] 
            lines=lines[1:] 

            for line in lines:
                fid.write('%.10f %s\n'%(mjdobs,line))
        fid.close()

def combinefisheye(params):

    filters = params["filters"]
    colors = params["colors"]

    for filter in filters:
    #for filter in []:
        #outpath = os.path.join(params["dirpathname"],filter)
        outpath = os.path.join(params["fitspath"],filter)

        data_all = []
        data_all_tph = []

        files = glob.glob(os.path.join(outpath,"*.long.%s.fits"%filter))
        #files = glob.glob(os.path.join(outpath,"*.short.%s.fish"%filter))
        files = sorted(files)

        runNumbers = []
        nstars = []
        mjds = []

        minmags = []
        maxmags = []
        medmags = []

        minsky = []
        maxsky = []
        medsky = []

        moonposx = []
        moonposy = []
        angles = []

        for file in files:

            runNumber = int(file.split('.')[-4])
            if runNumber > params["maxframes"]:
                continue

            fileprefix = file.replace(".fits","")
            outputfile = file.replace("fits","mch")
            tphfile = file.replace("fits","tph")

            if not os.path.isfile(outputfile):
                continue
            data_out = np.loadtxt(outputfile,comments='RA')

            if not os.path.isfile(tphfile):
                continue
            data_out_tph = np.loadtxt(tphfile)

            fitshdr_command = "fitshdr %s"%(file)
            p = subprocess.Popen(fitshdr_command.split(" "),stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            output, errors = p.communicate()

            for line in output.split("\n"):
                lineSplit = line.split(" ")
                lineSplit = [x for x in lineSplit if x != '']
                if len(lineSplit) == 0:
                    continue

                if lineSplit[0] == "BIAS":
                    bias=float(lineSplit[2])
                elif lineSplit[0] == "EXPTIME":
                    exptime=float(lineSplit[2])
                elif lineSplit[0] == "MJD-OBS":
                    mjdobs=float(lineSplit[2])
                elif lineSplit[0] == "NAXIS1":
                    nx=float(lineSplit[2])
                elif lineSplit[0] == "NAXIS2":
                    ny=float(lineSplit[2])

            runNumbers.append(runNumber)
            nstars.append(len(data_out))
            mjds.append(mjdobs)

            indexes = np.where(~np.isnan(data_out[:,17]))[0]
            magnitudes = data_out[indexes,17]
            magnitudes = sorted(magnitudes)

            index10 = int(np.floor(0.1*len(magnitudes)))
            index50 = int(np.floor(0.5*len(magnitudes)))
            index90 = int(np.floor(0.9*len(magnitudes)))

            minmags.append(magnitudes[index10])
            maxmags.append(magnitudes[index90])
            medmags.append(magnitudes[index50])

            indexes = np.where(~np.isnan(data_out[:,19]))[0]
            sky = data_out[indexes,19]
            sky = sorted(sky)

            index10 = int(np.floor(0.1*len(sky)))
            index50 = int(np.floor(0.5*len(sky)))
            index90 = int(np.floor(0.9*len(sky)))

            minsky.append(sky[index10])
            maxsky.append(sky[index90])
            medsky.append(sky[index50])

            #minmags.append(np.min(magnitudes))
            #maxmags.append(np.max(magnitudes))
            #medmags.append(np.median(magnitudes))

            index = np.argmax(data_out_tph[:,7])
            moonposx.append(data_out_tph[index,0])
            moonposy.append(data_out_tph[index,1]) 

            if len(data_out) == 0:
                continue            
            
            addruns = runNumber * np.ones((len(data_out),1))
            data_out = np.hstack((data_out,addruns))

            if len(data_all) == 0:
                data_all = data_out
                data_all_tph = data_out_tph
            else:
                data_all = np.vstack((data_all,data_out))
                data_all_tph = np.vstack((data_all_tph,data_out_tph))

            moon = ephem.Moon()
            t = astropy.time.Time(mjdobs,format='mjd',scale='utc')
            utc = t.unix
            utc_timestamp = datetime.utcfromtimestamp(utc)
            utc_date = utc_timestamp.strftime('%Y/%m/%d %H:%M:%S')
            moon.compute(utc_date)

            ra_moon = (360/(2*np.pi))*float(repr(moon.ra))
            dec_moon = (180/np.pi)*float(repr(moon.dec))

            ra1 = float(repr(moon.ra))
            ra2 = data_out[:,0] * ((2*np.pi)/360)
            d1 = float(repr(moon.dec))
            d2 = data_out[:,1] * ((2*np.pi)/360)

            cosA = np.sin(d1)*np.sin(d2) + np.cos(d1)*np.cos(d2)*np.cos(ra1-ra2)
            angle = np.arccos(cosA)*(360/(2*np.pi))
            angle = angle.T

            if len(angles) == 0:
                angles = angle
            else:
                angles = np.hstack((angles,angle))

            if runNumber == 850:
                angles_850 = angle
                data_out_850 = data_out
            if runNumber == 950:
                angles_950 = angle
                data_out_950 = data_out

            #if len(data_out) < 2000:
            #    print file
            #    print len(data_out)

        plotDir = os.path.join(params["fisheyeplotpath"],filter)
        if not os.path.isdir(plotDir):
            os.mkdir(plotDir)

        if len(data_all) == 0:
            continue

        tt = np.array(mjds)
        if len(tt) == 0:
            print "No available data for %s"%filter
            continue

        tt = tt - tt[0]
        tt = tt * 24

        plotName = os.path.join(plotDir,'nstars.png')
        plt.plot(tt,nstars,'*')
        plt.xlabel('Hours since twilight')
        plt.ylabel('number of stars')
        plt.show()
        plt.savefig(plotName,dpi=200)
        plt.close('all')

        plotName = os.path.join(plotDir,'mags.png')
        plt.plot(tt,minmags,'r*')
        plt.plot(tt,medmags,'g*')
        plt.plot(tt,maxmags,'b*')
        plt.xlabel('Hours since twilight')
        plt.ylabel('magnitudes [10th,50th,90th]')
        plt.show()
        plt.savefig(plotName,dpi=200)
        plt.close('all')

        plotName = os.path.join(plotDir,'sky.png')
        plt.plot(tt,minsky,'r*')
        plt.plot(tt,maxsky,'g*')
        plt.plot(tt,medsky,'b*')
        plt.xlabel('Hours since twilight')
        plt.ylabel('sky [10th,50th,90th]')
        plt.show()
        plt.savefig(plotName,dpi=200)
        plt.close('all')

        fileName = os.path.join(plotDir,'nstars.txt')
        f = open(fileName,'w')
        for runNumber,mjd,nstar in zip(runNumbers,mjds,nstars):
            f.write('%d %.10f %d\n'%(runNumber,mjd,nstar))
        f.close()

        indexes = np.where(~np.isnan(data_all[:,17]))[0]
        data_all = data_all[indexes,:]
        magnitudes = data_all[:,10] - data_all[:,17]
        angles = angles[indexes]

        #indexes = np.where(data_all[:,18] < 0.1)[0]
        #data_all = data_all[indexes,:]
        #magnitudes = data_all[:,10] - data_all[:,17]
        #magnitudes = data_all[:,10]
        #magnitudes = data_all[:,0] 

        plotName = os.path.join(plotDir,'angles.png')
        plot = plt.plot(angles,data_all[:,19],'k*')
        #plt.xlim([0,2897])
        #plt.ylim([0,1935])
        plt.xlabel('Angle [deg]')
        plt.ylabel('Sky brightness')
        #cb = plt.colorbar()
        #cb.set_label('minst')
        #cb.set_label('')
        plt.show()
        plt.savefig(plotName,dpi=200)
        plt.close('all')

        plotName = os.path.join(plotDir,'angles_850_950.png')
        plot = plt.plot(angles_850,data_out_850[:,19],'k*')
        plot = plt.plot(angles_950,data_out_950[:,19],'m^')
        #plt.xlim([0,2897])
        #plt.ylim([0,1935])
        plt.xlabel('Angle [deg]')
        plt.ylabel('Sky brightness')
        #cb = plt.colorbar()
        #cb.set_label('minst')
        #cb.set_label('')
        plt.show()
        plt.savefig(plotName,dpi=200)
        plt.close('all')

        print plotName
        print stop

        zenith = (90 - data_all[:,15]) * 2 * np.pi / 360.0
        secant = (1.002432 * np.cos(zenith)**2 + 0.148386 * np.cos(zenith) + 0.0096467) / (np.cos(zenith)**3 + 0.149864 * np.cos(zenith)**2 + 0.0102963 * np.cos(zenith) + 0.00303978)
        secant = 1/np.cos(zenith)

        plotName = os.path.join(plotDir,'mag_x_y.png')
        plot = plt.scatter(data_all[:,11],data_all[:,12],s=10,c=magnitudes,vmin=12,vmax=15,edgecolor='none')
        plt.xlim([0,2897])
        plt.ylim([0,1935])
        plt.xlabel('x')
        plt.ylabel('y')
        cb = plt.colorbar()
        #cb.set_label('minst')
        cb.set_label('minst - m')
        plt.show()
        plt.savefig(plotName,dpi=200)
        plt.close('all')

        #data_all[:,17] = data_all[:,17] + 18

        magnitudes = magnitudes[data_all[:,0].argsort()]
        data_all = data_all[data_all[:,0].argsort(),:]
        ar,return_index,return_inverse = np.unique(data_all[:,0],return_index=True, return_inverse=True)
        indexes = np.where(return_inverse == mode(return_inverse)[0])[0]
        data_all_cut = data_all[indexes,:]
        magnitudes_cut = magnitudes[indexes]

        indexes = data_all_cut[:,-1].argsort()
        data_all_cut = data_all_cut[indexes,:]
        magnitudes_cut = magnitudes_cut[indexes]

        radecdiff = np.sqrt((data_all[:,0] - data_all_cut[0,0])**2 + (data_all[:,1] - data_all_cut[0,1])**2)
        ar,return_index,return_inverse = np.unique(radecdiff,return_index=True, return_inverse=True)

        #indexes = np.where(return_inverse != mode(return_inverse)[0])[0]
        #indexes = np.where(return_inverse != return_inverse[0])[0]
        indexes = np.where(radecdiff != 0)[0]
        data_all_cut_2 = data_all[indexes,:]
        magnitudes_cut_2 = magnitudes[indexes]  
        ar,return_index,return_inverse = np.unique(data_all_cut_2[:,0],return_index=True, return_inverse=True)
        indexes = np.where(return_inverse == mode(return_inverse)[0])[0] 
        data_all_cut_2 = data_all_cut_2[indexes,:]
        magnitudes_cut_2 = magnitudes_cut_2[indexes]

        #data_all_cut = data_all_cut[100:500,:]
        #magnitudes_cut = magnitudes_cut[100:500]

        #indexes = np.where(return_inverse == index_2)[0]
        #print indexes

        #data_all_cut_2 = data_all[indexes,:]
        #magnitudes_cut_2 = magnitudes[indexes]

        xerror = data_all_cut[:,3] - data_all_cut[:,11] 
        yerror = data_all_cut[:,4] - data_all_cut[:,12]
        rerror = np.sqrt(xerror**2 + yerror**2)

        plotName = os.path.join(plotDir,'mag_xerror_max.png')
        plt.plot(xerror,magnitudes_cut,'*')
        plt.xlabel('x error')
        plt.ylabel('minst - m')
        plt.show()
        plt.savefig(plotName,dpi=200)
        plt.close('all')

        plotName = os.path.join(plotDir,'mag_xerror_max.png')
        plt.plot(yerror,magnitudes_cut,'*')
        plt.xlabel('y error')
        plt.ylabel('minst - m')
        plt.show()
        plt.savefig(plotName,dpi=200)
        plt.close('all')

        plotName = os.path.join(plotDir,'mag_rerror_max.png')
        plt.plot(rerror,magnitudes_cut,'*')
        plt.xlim([0,2])
        plt.xlabel('r error')
        plt.ylabel('minst - m')
        plt.show()
        plt.savefig(plotName,dpi=200)
        plt.close('all')

        plotName = os.path.join(plotDir,'mag_x_y_max.png')
        plt.scatter(data_all_cut[:,11],data_all_cut[:,12],s=10,c=magnitudes_cut,edgecolor='none')
        plt.xlabel('x')
        plt.ylabel('y')
        cb = plt.colorbar()
        cb.set_label('minst - m')
        plt.show()
        plt.savefig(plotName,dpi=200)
        plt.close('all')

        plotName = os.path.join(plotDir,'mag_sky.png')
        plt.plot(data_all[:,19],magnitudes,'*')
        plt.xlabel('sky')
        plt.ylabel('minst - m')
        plt.show()
        plt.savefig(plotName,dpi=200)
        plt.close('all')

        plotName = os.path.join(plotDir,'mag_time.png')
        plt.plot(data_all[:,21],magnitudes,'*')
        plt.xlabel('run numbers')
        plt.ylabel('minst - m')
        plt.show()
        plt.savefig(plotName,dpi=200)
        plt.close('all')

        plotName = os.path.join(plotDir,'mag_time_max.png')
        plt.plot(data_all_cut[:,21],magnitudes_cut,'*')
        plt.xlabel('run numbers')
        plt.ylabel('minst - m')
        plt.show()
        plt.savefig(plotName,dpi=200)
        plt.close('all')

        plotName = os.path.join(plotDir,'mag_time_error_max.png')
        plt.errorbar(data_all_cut[:,21],magnitudes_cut,yerr=data_all_cut[:,18])
        plt.xlabel('run numbers')
        plt.ylabel('minst - m')
        plt.ylim([11.5,14.0])
        plt.show()
        plt.savefig(plotName,dpi=200)
        plt.close('all')

        plotName = os.path.join(plotDir,'mag_sky_max.png')
        plt.plot(data_all_cut[:,19],magnitudes_cut,'*')
        plt.xlabel('sky')
        plt.ylabel('minst - m')
        plt.show()
        plt.savefig(plotName,dpi=200)
        plt.close('all')

        plotName = os.path.join(plotDir,'mag_time_max_2.png')
        plt.plot(data_all_cut[:,21],magnitudes_cut,'b*')
        plt.plot(data_all_cut_2[:,21],magnitudes_cut_2,'r*')
        plt.xlabel('run numbers')
        plt.ylabel('minst - m')
        plt.show()
        plt.savefig(plotName,dpi=200)
        plt.close('all')

        plotName = os.path.join(plotDir,'flux_chin.png')
        plt.semilogx(data_all_tph[:,7],data_all_tph[:,-1],'*')
        plt.xlabel('flux')
        plt.ylabel('chi/n')
        plt.show()
        plt.savefig(plotName,dpi=200)
        plt.close('all')

        plotName = os.path.join(plotDir,'flux_major.png')
        plt.semilogx(data_all_tph[:,7],data_all_tph[:,9],'*')
        plt.xlabel('flux')
        plt.ylabel('major')
        plt.show()
        plt.savefig(plotName,dpi=200)
        plt.close('all')

        plotName = os.path.join(plotDir,'flux_minor.png')
        plt.semilogx(data_all_tph[:,7],data_all_tph[:,10],'*')
        plt.xlabel('flux')
        plt.ylabel('minor')
        plt.show()
        plt.savefig(plotName,dpi=200)
        plt.close('all')

        plotName = os.path.join(plotDir,'major_minor.png')
        plt.scatter(data_all_tph[:,9],data_all_tph[:,10],s=10,c=np.log10(data_all_tph[:,7]),edgecolor='none',vmin=2,vmax=5)
        #plt.plot(data_all_tph[:,9],data_all_tph[:,10],'*')
        plt.xlabel('major')
        plt.ylabel('minor')
        plt.show()
        plt.savefig(plotName,dpi=200)
        plt.close('all')

        plotName = os.path.join(plotDir,'major_skyval.png')
        ax = plt.gca()
        plt.scatter(data_all_tph[:,9],data_all_tph[:,4],s=10,c=np.log10(data_all_tph[:,7]),edgecolor='none',vmin=2,vmax=5)
        #plt.plot(data_all_tph[:,9],data_all_tph[:,10],'*')
        ax.set_yscale('log')
        ax.set_xscale('log')
        plt.xlabel('major')
        plt.ylabel('skyval')
        plt.show()
        plt.savefig(plotName,dpi=200)
        plt.close('all')

        plotName = os.path.join(plotDir,'flux_skyfit.png')
        plt.loglog(data_all_tph[:,7],data_all_tph[:,6],'*')
        plt.xlabel('flux')
        plt.ylabel('skyfit')
        plt.show()
        plt.savefig(plotName,dpi=200)
        plt.close('all')

        z = np.polyfit(data_all_cut[:,16], magnitudes_cut, 1)
        p = np.poly1d(z)
        scatter = np.absolute(p(data_all_cut[:,16])-magnitudes_cut)
        scatter_sorted = np.sort(scatter)
        max_error = scatter_sorted[int(0.9*len(scatter_sorted))]
        indexes = np.where(scatter<=max_error)[0]
        z = np.polyfit(data_all_cut[indexes,16], magnitudes_cut[indexes], 1)
        p = np.poly1d(z)
        xp = np.linspace(np.min(data_all_cut[:,16]),np.max(data_all_cut[:,16]),1000)
        xp = np.linspace(1e2,1e4,1000)

        r_at_0 = p(0)

        plotName = os.path.join(plotDir,'mag_r_max.png')
        #plt.semilogx(data_all_cut[:,16],data_all_cut[:,17]+data_all_cut[:,20],'*')
        plt.semilogx(data_all_cut[:,16],magnitudes_cut,'*')
        plt.semilogx(xp, p(xp),'b-')
        plt.xlabel('r')
        plt.ylabel('minst - m')
        title_text = "minst - m at r=0: %.5f"%r_at_0
        plt.title(title_text)
        plt.show()
        plt.savefig(plotName,dpi=200)
        plt.close('all')

        plotName = os.path.join(plotDir,'mag_dm_max.png')
        #plt.semilogx(data_all_cut[:,16],data_all_cut[:,17]+data_all_cut[:,20],'*')
        plt.plot(data_all_cut[:,18],magnitudes_cut,'*')
        plt.xlabel('dm')
        plt.ylabel('minst - m')
        plt.show()
        plt.savefig(plotName,dpi=200)
        plt.close('all')

        zenith = (90 - data_all_cut[:,15]) * 2 * np.pi / 360.0
        secant = (1.002432 * np.cos(zenith)**2 + 0.148386 * np.cos(zenith) + 0.0096467) / (np.cos(zenith)**3 + 0.149864 * np.cos(zenith)**2 + 0.0102963 * np.cos(zenith) + 0.00303978) 
        #secant = 1/np.cos(zenith)

        #secant = 1/np.cos(zenith * 2 * np.pi / 360.0)
        z = np.polyfit(secant, magnitudes_cut, 1)
        #z = np.polyfit(secant, data_all_cut[:,17]+data_all_cut[:,20], 1)
        p = np.poly1d(z)
        scatter = np.absolute(p(secant)-magnitudes_cut)
        #scatter = np.absolute(p(secant)-data_all_cut[:,17]-data_all_cut[:,20])
        scatter_sorted = np.sort(scatter)
        max_error = scatter_sorted[int(0.9*len(scatter_sorted))]
        indexes = np.where(scatter<=max_error)[0]
        z = np.polyfit(secant[indexes], magnitudes_cut[indexes], 1)
        p = np.poly1d(z)
        xp = np.linspace(1,3,100)

        secant_at_0 = p(0)

        plotName = os.path.join(plotDir,'mag_secant_max.png')
        plt.plot(secant,magnitudes_cut,'*')
        #plt.plot(secant,data_all_cut[:,17]+data_all_cut[:,20],'*')
        plt.plot(xp, p(xp),'b-')
        plt.xlim([1,3])
        plt.xlabel('secant')
        plt.ylabel('minst - m')
        title_text = "minst - m at secant=0: %.5f"%secant_at_0
        plt.title(title_text)
        plt.show()
        plt.savefig(plotName,dpi=200)
        plt.close('all')

        plotName = os.path.join(plotDir,'mag_modx_max.png')
        plt.plot(np.mod(data_all_cut[:,11],1),magnitudes_cut+secant,'*')
        plt.xlim([-0.05,0.95])
        plt.xlabel('mod(x,1)')
        plt.ylabel('minst - m')
        plt.show()
        plt.savefig(plotName,dpi=200)
        plt.close('all')

        plotName = os.path.join(plotDir,'mag_mody_max.png')
        plt.plot(np.mod(data_all_cut[:,12],1),magnitudes_cut+secant,'*')
        plt.xlim([-0.05,0.95])
        plt.xlabel('mod(y,1)')
        plt.ylabel('minst - m')
        plt.show()
        plt.savefig(plotName,dpi=200)
        plt.close('all')

        plotName = os.path.join(plotDir,'mag_modx_mody_max.png')
        plt.scatter(np.mod(data_all_cut[:,11],1),np.mod(data_all_cut[:,12],1),s=10,c=magnitudes_cut+secant,edgecolor='none')
        plt.xlim([-0.05,0.95])
        plt.ylim([-0.05,0.95])
        plt.xlabel('mod(x,1)')
        plt.ylabel('mod(y,1)')
        cb = plt.colorbar()
        cb.set_label('minst - m')
        plt.show()
        plt.savefig(plotName,dpi=200)
        plt.close('all')

        m_inst = data_all_cut[:,17]
        b_minus_v = data_all_cut[:,8]
        airmass = secant
        f_of_r = data_all_cut[:,16]
        ydata = np.zeros(data_all_cut[:,17].shape)

        alpha = 0
        beta = 0
        gamma = 0
        m_0 = -4
        p0 = (alpha,beta,gamma,m_0)

        x_4d = np.vstack((m_inst,b_minus_v, airmass-airmass, f_of_r)).T
        popt, pcov = curve_fit(m_func, x_4d, ydata, p0=p0)
        y_fit = m_func(x_4d, popt[0], popt[1], popt[2], popt[3])

        plotName = os.path.join(plotDir,'mag_sec_max_fit.png')
        plt.plot(secant,y_fit,'*')
        #plt.xlim([-0.05,0.95])
        #plt.xlabel('mod(y,1)')
        #plt.ylabel('minst - m')
        plt.show()
        plt.savefig(plotName,dpi=200)
        plt.close('all')

        plotName = os.path.join(plotDir,'mag_time_max_fit.png')
        plt.plot(data_all_cut[:,21],y_fit,'*')
        #plt.xlim([-0.05,0.95])
        #plt.xlabel('mod(y,1)')
        #plt.ylabel('minst - m')
        plt.show()
        plt.savefig(plotName,dpi=200)
        plt.close('all')
 
        unique_return_inverse = np.unique(return_inverse)

        plotName = os.path.join(plotDir,'mag_r.png')
        for i in xrange(len(unique_return_inverse)):
            index = unique_return_inverse[i]

            indexes = np.where(return_inverse == index)[0]

            if len(indexes) < 10:
                continue

            data_all_cut = data_all[indexes,:]
            magnitudes_cut = magnitudes[indexes]

            plt.semilogx(data_all_cut[:,16],magnitudes_cut,'*')

        plt.xlabel('r')
        plt.ylabel('minst - m')
        plt.show()
        plt.savefig(plotName,dpi=200)
        plt.close('all')

        rmserrors = []

        variableStars = getVariableStars()

        filename = os.path.join(plotDir,"mags.txt")
        fid = open(filename,"w+")

        plotName = os.path.join(plotDir,'mag_secant.png')
        for i in xrange(len(unique_return_inverse)):
            index = unique_return_inverse[i]

            indexes = np.where(return_inverse == index)[0]

            if len(indexes) < 10:
                continue

            data_all_cut = data_all[indexes,:]
            magnitudes_cut = magnitudes[indexes]
            indexes = np.where(~np.isnan(data_all_cut[:,15]))[0]
            data_all_cut = data_all_cut[indexes,:]         
            magnitudes_cut = magnitudes_cut[indexes]

            varDiff = np.sqrt((variableStars[:,0] - data_all_cut[0,0])**2 + \
                (variableStars[:,1] - data_all_cut[0,1])**2)

            if np.min(varDiff) < 0.1:
                continue

            zenith = (90 - data_all_cut[:,15]) * 2 * np.pi / 360.0
            secant = (1.002432 * np.cos(zenith)**2 + 0.148386 * np.cos(zenith) + 0.0096467) / (np.cos(zenith)**3 + 0.149864 * np.cos(zenith)**2 + 0.0102963 * np.cos(zenith) + 0.00303978)
            secant = 1/np.cos(zenith)

            #secant = 1/np.cos(zenith * 2 * np.pi / 360.0)
            z = np.polyfit(secant, magnitudes_cut, 1)
            p = np.poly1d(z)
            scatter = np.absolute(p(secant)-magnitudes_cut)
            scatter_sorted = np.sort(scatter)
            max_error = scatter_sorted[int(0.9*len(scatter_sorted))]

            indexes = np.where(scatter<=max_error)[0]
            z = np.polyfit(secant[indexes], magnitudes_cut[indexes], 1)
            p = np.poly1d(z)
            xp = np.linspace(1,3,100)

            #fit m_com | r=0 = m_inst - (alpha + beta (m_g - m_z))*airmass + gamma*f(r) - m_effect (-14)
            #magnitudes = data_all[:,17]

            plt.plot(secant,magnitudes_cut,'*')
            #plt.plot(xp, p(xp),'b-')

            #print magnitudes_cut
            #print p(magnitudes_cut)
            rmserror = np.mean(np.absolute(magnitudes_cut - p(secant)))
            #if rmserror > 1.5:
            #    print data_all_cut[0,0],data_all_cut[0,1],p
            #    print magnitudes_cut

            #p = np.median(magnitudes_cut)

            rmserrors.append(rmserror)
            #plt.plot(data_all_cut[0,0],rmserror,'*')

            m_inst = data_all_cut[:,17]
            b_minus_v = data_all_cut[:,8]
            airmass = secant
            f_of_r = data_all_cut[:,16]
            ydata = np.zeros(data_all_cut[:,17].shape)

            alpha = 0
            beta = 0
            gamma = 0
            m_0 = -4

            p0 = (alpha,beta,gamma,m_0)
            x_4d = np.vstack((m_inst,b_minus_v, airmass, f_of_r)).T
            popt, pcov = curve_fit(m_func, x_4d, ydata, p0=p0)

            print popt

            fid.write("%.5f %.5f %.10e %.10e\n"%(data_all_cut[0,0],data_all_cut[0,1],p[0],p[1]))

        fid.close()

        plt.xlabel('secant')
        plt.ylabel('minst')
        plt.xlim([1,3])
        plt.show()
        plt.savefig(plotName,dpi=200)
        plt.close('all')

    plotDir = os.path.join(params["fisheyeplotpath"],"colorcompare")
    if not os.path.isdir(plotDir):
        os.mkdir(plotDir)
    for filter1 in filters:
        for filter2 in filters:

            magsDir = os.path.join(params["fisheyeplotpath"],filter1)
            filename = os.path.join(magsDir,"mags.txt")

            if not os.path.isfile(filename):
                continue

            mags_data_1 = np.loadtxt(filename)

            if len(mags_data_1) == 0:
                continue

            magsDir = os.path.join(params["fisheyeplotpath"],filter2)
            filename = os.path.join(magsDir,"mags.txt")

            if not os.path.isfile(filename):
                continue

            mags_data_2 = np.loadtxt(filename)

            if len(mags_data_2) == 0:
                continue

            data1 = []
            data2 = []

            for i in xrange(len(mags_data_1)):

                ra = mags_data_1[i,0]
                dec = mags_data_1[i,1]
                index1 = np.where(ra == mags_data_2[:,0])[0]
                index2 = np.where(dec == mags_data_2[:,1])[0]

                if len(index1) == 0:
                    continue

                index = index1[0]

                #if not len(mags_data_2[index,2]) == 0:
 
                data1.append(mags_data_1[i,2])
                data2.append(mags_data_2[index,2])

            plotName = os.path.join(plotDir,'%s_%s.png'%(filter1,filter2))
            plt.plot(data1,data2,'*')
            plt.xlabel(filter1)
            plt.ylabel(filter2)
            plt.xlim([13,15])
            plt.ylim([13,15])
            plt.show()
            plt.savefig(plotName,dpi=200)
            plt.close('all')

def airplane(params):

    filters = params["filters"]
    colors = params["colors"]
    data = {}

    for filter in filters:
        outpath = os.path.join(params["dirpathname"],filter)
        outpath = os.path.join(params["fitspath"],filter)

        files = glob.glob(os.path.join(outpath,"*.long.%s.median.fits"%filter))
        files = sorted(files)

        runNumbers = []
        mjdobss = []
        skybrightness = []

        baseplotDir = os.path.join(params["airplanepath"],filter)
        if not os.path.isdir(baseplotDir):
            os.mkdir(baseplotDir)

        for file in files:

            folderName = file.split('.')[-5]
            runNumber = int(folderName)

            #if not runNumber == 898:
            #    continue

            if not runNumber == 700:
                continue

            if runNumber > params["maxframes"]:
                continue

            fitshdr_command = "fitshdr %s"%(file)
            p = subprocess.Popen(fitshdr_command.split(" "),stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            output, errors = p.communicate()

            for line in output.split("\n"):
                lineSplit = line.split(" ")
                lineSplit = [x for x in lineSplit if x != '']
                if len(lineSplit) == 0:
                    continue

                if lineSplit[0] == "BIAS":
                    bias=float(lineSplit[2])
                elif lineSplit[0] == "EXPTIME":
                    exptime=float(lineSplit[2])
                elif lineSplit[0] == "MJD-OBS":
                    mjdobs=float(lineSplit[2])
                elif lineSplit[0] == "NAXIS1":
                    nx=float(lineSplit[2])
                elif lineSplit[0] == "NAXIS2":
                    ny=float(lineSplit[2])

            hdulist = astropy.io.fits.open(file)
            scidata = hdulist[0].data
            im = scidata

            #theta = np.linspace(0., 180., max(im.shape), endpoint=True)
            theta = np.linspace(0.0,180.,181)
            sinogram = radon(im, theta=theta, circle=False)

            plotDir = os.path.join(baseplotDir,folderName)
            if not os.path.isdir(plotDir):
                os.mkdir(plotDir)

            plotName = os.path.join(plotDir,'radon.png')
            plt.imshow(sinogram, cmap=plt.cm.Greys_r,
                extent=(0, 180, 0, sinogram.shape[0]), aspect='auto')
            plt.xlabel('Projection angle [deg]')
            plt.ylabel('Projection position [pixels]')
            plt.show()
            plt.savefig(plotName,dpi=200)
            plt.close('all')

            blobs = blob_doh(im, max_sigma=10, threshold=.000001)

            plotName = os.path.join(plotDir,'blobs.png')
            fig, ax = plt.subplots(1, 1)
            plt.imshow(im, cmap=plt.cm.Greys_r,
                extent=(0, sinogram.shape[0], 0, sinogram.shape[1])) #aspect='auto')
            color = 'r'
            print blobs
            for blob in blobs:
                y, x, r = blob
                print y, x, r
                c = plt.Circle((x, y), r, color=color, linewidth=2, fill=True)
                fig.gca().add_artist(c)
                #ax.add_patch(c)
            #plt.xlim([x - 10,x + 10])
            #plt.ylim([y - 10,y + 10])
            #plt.xlabel('Projection angle [deg]')
            #plt.ylabel('Projection position [pixels]')
            plt.show()
            plt.savefig(plotName,dpi=200)
            plt.close('all')

            print plotName
            print stop

            runNumbers.append(runNumber)
            skybrightness.append(skyval)
            mjdobss.append(mjdobs)

        plotName = os.path.join(plotDir,'skybrightness.png')
        plt.semilogy(runNumbers,skybrightness,'*')
        plt.xlabel('run number')
        plt.ylabel('skybrightness')
        plt.title(params["outputFolder"])
        plt.show()
        plt.savefig(plotName,dpi=200)
        plt.close('all')

        f = open(os.path.join(plotDir,'skybrightness.txt'),'w')
        for mjdobs,sky in zip(mjdobss,skybrightness):
            f.write("%.10f %.10e\n"%(mjdobs,sky))
        f.close()

        data[filter] = {}
        data[filter]["mjdobs"] = mjdobss
        data[filter]["skybrightness"] = skybrightness

    f = open(os.path.join(params["skybrightnessplotpath"],'skybrightness.txt'),'w')
    for mjdobs,m,r,g,b in zip(data["M"]["mjdobs"],data["M"]["skybrightness"],data["R"]["skybrightness"],data["G"]["skybrightness"],data["B"]["skybrightness"]):
        f.write("%.10f %.10e %.10e %.10e %.10e\n"%(mjdobs,m,r,g,b))
    f.close()

def skybrightness(params):

    filters = params["filters"]
    colors = params["colors"]
    data = {}

    for filter in filters:
        outpath = os.path.join(params["dirpathname"],filter)
        outpath = os.path.join(params["fitspath"],filter)

        files = glob.glob(os.path.join(outpath,"*.long.%s.median.fits"%filter))
        files = sorted(files)

        runNumbers = []
        mjdobss = []
        skybrightness = []

        for file in files:

            runNumber = int(file.split('.')[-5])
            if runNumber > params["maxframes"]:
                continue

            fitshdr_command = "fitshdr %s"%(file)
            p = subprocess.Popen(fitshdr_command.split(" "),stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            output, errors = p.communicate()

            for line in output.split("\n"):
                lineSplit = line.split(" ")
                lineSplit = [x for x in lineSplit if x != '']
                if len(lineSplit) == 0:
                    continue

                if lineSplit[0] == "BIAS":
                    bias=float(lineSplit[2])
                elif lineSplit[0] == "EXPTIME":
                    exptime=float(lineSplit[2])
                elif lineSplit[0] == "MJD-OBS":
                    mjdobs=float(lineSplit[2])
                elif lineSplit[0] == "NAXIS1":
                    nx=float(lineSplit[2])
                elif lineSplit[0] == "NAXIS2":
                    ny=float(lineSplit[2])

            hdulist = astropy.io.fits.open(file)
            scidata = hdulist[0].data
            im = scidata

            # ny is actually nx, vice versa
            index_x = int(np.ceil(ny/2))
            index_y = int(np.ceil(nx/2))

            skyval = (im[index_x,index_y]-bias)/exptime
            runNumbers.append(runNumber)
            skybrightness.append(skyval)
            mjdobss.append(mjdobs)

        plotDir = os.path.join(params["skybrightnessplotpath"],filter)
        if not os.path.isdir(plotDir):
            os.mkdir(plotDir)

        plotName = os.path.join(plotDir,'skybrightness.png')
        plt.semilogy(runNumbers,skybrightness,'*')
        plt.xlabel('run number')
        plt.ylabel('skybrightness')
        plt.title(params["outputFolder"])
        plt.show()
        plt.savefig(plotName,dpi=200)
        plt.close('all')

        f = open(os.path.join(plotDir,'skybrightness.txt'),'w')
        for mjdobs,sky in zip(mjdobss,skybrightness):
            f.write("%.10f %.10e\n"%(mjdobs,sky))
        f.close()

        data[filter] = {}
        data[filter]["mjdobs"] = mjdobss
        data[filter]["skybrightness"] = skybrightness

    f = open(os.path.join(params["skybrightnessplotpath"],'skybrightness.txt'),'w')
    for mjdobs,m,r,g,b in zip(data["M"]["mjdobs"],data["M"]["skybrightness"],data["R"]["skybrightness"],data["G"]["skybrightness"],data["B"]["skybrightness"]):
        f.write("%.10f %.10e %.10e %.10e %.10e\n"%(mjdobs,m,r,g,b))
    f.close()

def fisheye(params):

    filters = params["filters"]
    colors = params["colors"]

    for filter in filters:
        outpath = os.path.join(params["dirpathname"],filter)
        outpath = os.path.join(params["fitspath"],filter)

        files = glob.glob(os.path.join(outpath,"*.%s.fits"%filter))
        files = sorted(files)

        for file in files:

            runNumber = int(file.split('.')[-4])
            if runNumber > params["maxframes"]:
                continue

            #if not runNumber == 300:
            #    continue

            fileprefix = file.replace(".fits","")
            outputfile = file.replace("fits","fish")
            if os.path.isfile(outputfile):
                continue

            fitshdr_command = "fitshdr %s"%(file)
            p = subprocess.Popen(fitshdr_command.split(" "),stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            output, errors = p.communicate()

            for line in output.split("\n"):
                lineSplit = line.split(" ")
                lineSplit = [x for x in lineSplit if x != '']
                if len(lineSplit) == 0:
                    continue

                if lineSplit[0] == "MJD-OBS":
                    utcobs=float(lineSplit[2])

            type = file.split('.')[-3]
            if utcobs <= 56800:
                lens_type = "old"
            else:
                lens_type = "new"

            imstats_file = "%s/imstats_fish/imstats_fish_%s_%s.sh"%(params["codepath"],type,lens_type)
            fisheye_command = "%s %s"%(imstats_file,fileprefix)
            os.system(fisheye_command)

            if not os.path.isfile(outputfile):
                continue

            data_out = np.loadtxt(outputfile)

            if len(data_out) == 0:
                continue 

            if data_out.ndim == 1:
                print "only one entry in file... continuing"
                continue

            hist_vals = data_out[:,-1]

            plotDir = os.path.join(params["fisheyeplotpath"],file.split(".")[1])
            if not os.path.isdir(plotDir):
                os.mkdir(plotDir)

            plotName = os.path.join(plotDir,'hist_%s.png'%filter)
            hist, bins = np.histogram(hist_vals, bins=50)

            width = 0.7 * (bins[1] - bins[0])
            center = (bins[:-1] + bins[1:]) / 2
            plt.bar(center, hist, align='center', width=width)
            plt.show()
            plt.savefig(plotName,dpi=200)
            plt.close('all')

def photometry_tphot(params):

    outfile_nstars = os.path.join(params["dirpathname"],"%s.nstars"%(params["dirname"]))

    filters = params["filters"]
    colors = params["colors"]

    for filter in filters:
        outpath = os.path.join(params["dirpathname"],filter)

        rm_command = "rm %s/*.phot"%(outpath)
        #os.system(rm_command)
        rm_command = "rm %s/*.nstars"%(outpath)
        #os.system(rm_command)

        files = glob.glob(os.path.join(outpath,"*.%s.fits"%filter))
        for file in files:
            outfile = file.replace("fits","tph")

            runNumber = int(file.split('.')[-4])
            if runNumber > params["maxframes"]:
                continue 

            if not runNumber == 200:
                continue

            tphot_command = "imagebias=`gethead %s BIAS`; tphot %s -out %s -bias $imagebias"%(file,file,outfile)
            os.system(tphot_command)


def photometry(params):

    filters = params["filters"]
    colors = params["colors"]

    for filter in filters:
        expath = os.path.join(params["expathname"],filter)

        if not os.path.isdir(expath):
            os.mkdir(expath)

    for filter in filters:
        outpath = os.path.join(params["dirpathname"],filter)
        expath = os.path.join(params["expathname"],filter)

        if not os.path.isdir(expath):
            os.mkdir(expath)

        vals = []

        files = glob.glob(os.path.join(outpath,"*.%s.fits"%filter))
        for fitsfile in files:

            plotDir = os.path.join(params["fitsplotpath"],fitsfile.split(".")[2])
            if not os.path.isdir(plotDir):
                os.mkdir(plotDir)

            plotName = os.path.join(plotDir,'fits_%s.png'%filter)
            gc = aplpy.FITSFigure(fitsfile)
            gc.show_grayscale()
            gc.tick_labels.set_font(size='small')
            gc.save(plotName)
            gc.close()

            hdulist = astropy.io.fits.open(fitsfile)
            hdulist.info()
            scidata = hdulist[0].data
            im = scidata

            #indexes = np.where(im > 14000)[0]
            #print indexes
            #for index in indexes:
            #    im[index] = 0

            plotName = os.path.join(plotDir,'hist_%s.png'%filter)
            hist, bins = np.histogram(im[:], bins=50)

            width = 0.7 * (bins[1] - bins[0])
            center = (bins[:-1] + bins[1:]) / 2
            plt.bar(center, hist, align='center', width=width)
            plt.show()
            plt.savefig(plotName,dpi=200)
            plt.close('all')

            exfile = os.path.join(expath,fitsfile.split("/")[-1])
            ex_command = "sex -c allsky.sex %s -CATALOG_NAME %s"%(fitsfile,exfile)
            os.system(ex_command)

def makefits(params):

    filters = params["filters"]
    colors = params["colors"]

    for filter in filters:
        #outpath = os.path.join(params["dirpathname"],filter)
        outpath = os.path.join(params["fitspath"],filter)

        if not os.path.isdir(outpath):
            os.mkdir(outpath)
        rm_command = "rm %s/*.fits"%outpath
        #os.system(rm_command)

    print "Separating files into CR2 and color subdirectories"
    #outpath = os.path.join(params["dirpathname"],"CR2")
    outpath = os.path.join(params["fitspath"],"CR2")
    if not os.path.isdir(outpath):
         os.mkdir(outpath)

    framecounter = 1
    for i in xrange(params["maxframes"]):

        filename_short_gz = os.path.join(params["datapath"],"%s.%04d.short.cr2.gz"%(params["dirname"],framecounter))
        filename_short = os.path.join(outpath,"%s.%04d.short.cr2"%(params["dirname"],framecounter))
        framecounter = framecounter + 1
        filename_long_gz = os.path.join(params["datapath"],"%s.%04d.long.cr2.gz"%(params["dirname"],framecounter))
        filename_long = os.path.join(outpath,"%s.%04d.long.cr2"%(params["dirname"],framecounter))
        framecounter = framecounter + 1

        #if not framecounter == 201:
        #    continue

        if os.path.isfile(filename_short) and os.path.isfile(filename_long):
            continue

        if not (os.path.isfile(filename_short_gz) and os.path.isfile(filename_long_gz)):
            continue

        gzip_command = 'gzip -c -d %s > %s'%(filename_short_gz,filename_short)
        os.system(gzip_command)
        gzip_command = 'gzip -c -d %s > %s'%(filename_long_gz,filename_long)
        os.system(gzip_command)

        for filter,color in zip(filters,colors):

            print "Extracting %s FITS from CR2, raw2fits -%s ..."%(filter,color)
            #outpath = os.path.join(params["dirpathname"],filter)
            outpath = os.path.join(params["fitspath"],filter)

            plotdir = os.path.join(outpath,"plots")
            if not os.path.isdir(plotdir):
                os.mkdir(plotdir)

            cr2fits_command = "raw2fits -dir %s -%s %s"%(outpath,color,filename_short)
            os.system(cr2fits_command)
            cr2fits_command = "raw2fits -dir %s -%s %s"%(outpath,color,filename_long)
            os.system(cr2fits_command)

            filename_short_split = filename_short.split("/")
            fitsfile_short = os.path.join(outpath,filename_short_split[-1].replace("cr2","fits"))
            filename_long_split = filename_long.split("/")
            fitsfile_long = os.path.join(outpath,filename_long_split[-1].replace("cr2","fits"))

            mv_command = "mv %s %s"%(fitsfile_short,fitsfile_short.replace("fits","%s.fits"%filter))
            os.system(mv_command)
            fitsfile_short = fitsfile_short.replace("fits","%s.fits"%filter)
            mv_command = "mv %s %s"%(fitsfile_long,fitsfile_long.replace("fits","%s.fits"%filter))
            os.system(mv_command)
            fitsfile_long = fitsfile_long.replace("fits","%s.fits"%filter)

            fitsfile_short_median = fitsfile_short.replace("fits","median.fits")
            hdulist = astropy.io.fits.open(fitsfile_short)
            im = hdulist[0].data
            gethead_command = "gethead %s BIAS"%(fitsfile_short)
            p = subprocess.Popen(gethead_command.split(" "),stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            output, errors = p.communicate()
            bias = float(output)

            side = 10
            size=(side,side)
            im_filt = scipy.ndimage.filters.median_filter(im, size=size, mode='reflect')

            hdulist = astropy.io.fits.open(fitsfile_short)
            hdulist[0].data = im_filt
            hdulist.writeto(fitsfile_short_median,output_verify='warn',clobber=True)

            #if np.mod(i,50) == 0:
            if True:
                fitsSplit = fitsfile_short.split("/")
                plotName = os.path.join(plotdir,'%s.png'%fitsSplit[-1])

                gc = aplpy.FITSFigure(fitsfile_short)
                gc.show_grayscale()
                gc.tick_labels.set_font(size='small')
                gc.save(plotName)
                gc.close()

                fitsSplit = fitsfile_short_median.split("/")
                plotName = os.path.join(plotdir,'%s.png'%fitsSplit[-1])
                gc._data = im_filt
                gc.show_grayscale()
                gc.tick_labels.set_font(size='small')
                gc.save(plotName)
                gc.close()

            fitsfile_long_median = fitsfile_long.replace("fits","median.fits")
            hdulist = astropy.io.fits.open(fitsfile_long)
            im = hdulist[0].data
            gethead_command = "gethead %s BIAS"%(fitsfile_long)
            p = subprocess.Popen(gethead_command.split(" "),stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            output, errors = p.communicate()
            bias = float(output)

            side = 10
            size=(side,side)
            im_filt = scipy.ndimage.filters.median_filter(im, size=size, mode='reflect')

            hdulist = astropy.io.fits.open(fitsfile_long)
            hdulist[0].data = im_filt
            hdulist.writeto(fitsfile_long_median,output_verify='warn',clobber=True)

            #if np.mod(i,50) == 0:
            if True:
                fitsSplit = fitsfile_long.split("/")
                plotName = os.path.join(plotdir,'%s.png'%fitsSplit[-1])
                gc = aplpy.FITSFigure(fitsfile_long)
                gc.show_grayscale()
                gc.tick_labels.set_font(size='small')
                gc.save(plotName)
                gc.close()

                fitsSplit = fitsfile_long_median.split("/")
                plotName = os.path.join(plotdir,'%s.png'%fitsSplit[-1])
                gc._data = im_filt
                gc.show_grayscale()
                gc.tick_labels.set_font(size='small')
                gc.save(plotName)
                gc.close()

def makemovie(params):

    filters = params["filters"]
    colors = params["colors"]

    for filter,color in zip(filters,colors):

        #outpath = os.path.join(params["dirpathname"],filter)
        outpath = os.path.join(params["fitspath"],filter)
        plotdir = os.path.join(outpath,"plots")
        if not os.path.isdir(plotdir):
            os.mkdir(plotdir)

        outpath = os.path.join(params["dirpathname"],filter)
        if not os.path.isdir(outpath):
            os.mkdir(outpath)
        moviedir = os.path.join(outpath,"movie")
        if not os.path.isdir(moviedir):
            os.mkdir(moviedir)

        files = sorted(glob.glob(os.path.join(plotdir,'*.long.%s.fits.png'%filter)))
        n=1
        for file in files:
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

        files = sorted(glob.glob(os.path.join(plotdir,'*.long.%s.median.fits.png'%filter)))
        n=1
        for file in files:
            n = n + 1
            filename = os.path.join(moviedir,"fishy-%04d.png"%n)
            cp_command = "cp %s %s"%(file,filename)
            os.system(cp_command)

        moviefiles = os.path.join(moviedir,"fishy-%04d.png")
        filename = os.path.join(moviedir,"long_median.mpg")
        ffmpeg_command = 'ffmpeg -an -y -r 20 -i %s -b:v %s %s'%(moviefiles,'5000k',filename)
        os.system(ffmpeg_command)
        moviefiles = os.path.join(moviedir,"fishy-%04d.png")
        filename = os.path.join(moviedir,"long_median.gif")
        ffmpeg_command = 'ffmpeg -an -y -r 20 -i %s -b:v %s %s'%(moviefiles,'5000k',filename)
        os.system(ffmpeg_command)

        rm_command = "rm %s/*.png"%(moviedir)
        os.system(rm_command)

