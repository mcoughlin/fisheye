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
from matplotlib import cm

from matplotlib.mlab import griddata

import ephem
from utils import calcDist_haversine, mjd2djd, radec2pix, healbin
import healpy as hp

from grabESO import grabESO

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
    parser.add_option("--doCompareCloud",  action="store_true", default=False)
    parser.add_option("--doSkyCatalog",  action="store_true", default=False)
    parser.add_option("--doSkyCatalogInd",  action="store_true", default=False)
    parser.add_option("--doCombineSkyCatalog",  action="store_true", default=False)
    parser.add_option("--doFieldSkyCatalog",  action="store_true", default=False)
    parser.add_option("--doEphemSkyCatalog",  action="store_true", default=False)
    parser.add_option("--doZeropointSkyCatalog",  action="store_true", default=False)
    parser.add_option("--doCombineZeropointSkyCatalog",  action="store_true", default=False)
    parser.add_option("--doAnalyzeEphemSkyCatalog",  action="store_true", default=False)
    parser.add_option("--doMoonInfo",  action="store_true", default=False)

    parser.add_option("-r", "--ra", help="right ascension",type=float,default=0.4455) 
    parser.add_option("-d", "--dec", help="declination",type=float,default=-50.2887)

    opts, args = parser.parse_args()

    return opts

def analyzeephemskycatalog(params):

    filters = params["filters"]
    colors = params["colors"]

    for filter,color in zip(filters,colors):

        plotDir = os.path.join(params["analyzeephemplotpath"],filter)
        if not os.path.isdir(plotDir):
            os.mkdir(plotDir)

        # Read in all the photometry and package it up and save it
        # Set up LSST observatory
        filename = os.path.join(params["ephemplotpath"],"skyData_%s.npz"%filter)
        data_out = np.load(filename)
        skyData = data_out['skyData']

        # Let's look at some of the photometry data
        nside = 16
    
        #good=np.where(skyData['mjd'] == skyData['mjd'][0])[0]
        #plt.scatter(skyData['ra'][good], skyData['dec'][good], c=skyData['sky'][good])
    
        hpid = radec2pix(nside, np.radians(skyData['ra']),np.radians(skyData['dec']))

        sids = np.unique(hpid)

        count = 0
        for sid in sids:

            continue

            plt.figure()
            good=np.where((hpid == sid) & (skyData['sunAlt'] < np.radians(-20.)))[0]
            plt.scatter(np.degrees(skyData['moonAlt'][good]), skyData['moonPhase'][good], c=skyData['sky'][good], alpha=.1)
            plt.xlabel('Moon Altitude (deg)')
            plt.ylabel('Moon Phase')
            plt.title('Sun down, healpix id = %i'%sid)
            cb = plt.colorbar()
            cb.set_label('Sky brightness (mags/area)')
            #plotName = os.path.join(plotDir,'moonSB_%i.png'%sid)
            plotName = os.path.join(plotDir,'moonSB_%i.png'%count)
            plt.savefig(plotName)
    
            plt.figure()
            good=np.where((hpid == sid) &
                (skyData['alt'] > np.radians(10.)) &
                (skyData['sunAlt'] < np.radians(-20.)) &
                (skyData['moonAlt'] < np.radians(-20.)))[0]
    
            plt.plot(np.degrees(skyData['alt'][good]), skyData['sky'][good], 'ko', alpha=.1)
            plt.gca().invert_yaxis() # mags!
            plt.xlabel('Altitude (degrees)')
            plt.ylabel('Sky (mags/area)')
            plt.title('Moon and Sun down, healpix id = %i'%sid)
            #plotName = os.path.join(plotDir,'altSB_%i.png'%sid)
            plotName = os.path.join(plotDir,'altSB_%i.png'%count)
            plt.savefig(plotName)
    
            plt.figure()
            good=np.where((hpid == sid) &
              (skyData['sunAlt'] > np.radians(-20)) )[0]
            plt.scatter(np.degrees(skyData['sunAlt'][good]), skyData['sky'][good],
                c=np.degrees(skyData['moonAlt'][good]), alpha=.5)
            plt.gca().invert_yaxis() # mags!
            plt.xlabel('Sun Altitude (degrees)')
            plt.ylabel('Sky (mags/area)')
            plt.title('Sun up, healpix id = %i'%sid)
            cb = plt.colorbar()
            cb.set_label('Moon Alt (deg)')
            #plotName = os.path.join(plotDir,'sunSB_%i.png'%sid)
            plotName = os.path.join(plotDir,'sunSB_%i.png'%count)
            plt.savefig(plotName)

            count = count + 1

        moviedir = os.path.join(plotDir,"movie")
        if not os.path.isdir(moviedir):
            os.mkdir(moviedir)

        count = 0
        mjds = np.unique(skyData['mjd'])[0:1000]
        for mjd in mjds:

            #continue

            # Try plotting a single sky shot:
            good = np.where(skyData['mjd'] == mjd)[0]
 
            # So how do we model things?  Each HP gets a base spectrum, and then some added flux based on increased airmass, moon, and sun?

            dayMap = healbin(np.radians(skyData['ra'][good]), np.radians(skyData['dec'][good]), skyData['sky'][good], nside=nside)

            hp.mollview(dayMap, title='single Frame MJD=%f'%skyData['mjd'][good][0], unit='Surface Brightness (mag/area)')
            plotName = os.path.join(moviedir,'singleFrameMap-%04d.png'%count)
            plt.savefig(plotName)
           
            count = count + 1

        moviefiles = os.path.join(moviedir,"singleFrameMap-%04d.png")
        filename = os.path.join(moviedir,"singleFrameMap.mpg")
        ffmpeg_command = 'ffmpeg -an -y -r 20 -i %s -b:v %s %s'%(moviefiles,'5000k',filename)
        os.system(ffmpeg_command)
        filename = os.path.join(moviedir,"singleFrameMap.gif")
        ffmpeg_command = 'ffmpeg -an -y -r 20 -i %s -b:v %s %s'%(moviefiles,'5000k',filename)
        os.system(ffmpeg_command)

        rm_command = "rm %s/*.png"%(moviedir)
        #os.system(rm_command)

        # Let's make a minimum brightness map over the sky
        minSkyMap = np.ma.masked_all(hp.nside2npix(nside), dtype=float)
        minSkyMap.fill_value = hp.UNSEEN
        maxAltMap = minSkyMap.copy()

        altPad = 5. # Degrees
        perClip = 90. # Percentile clip to get rid of clouds and outliers

        order = np.argsort(hpid)
        skyData = skyData[order]
        hpid = hpid[order]
        ids = np.arange(hp.nside2npix(nside))
        left = np.searchsorted(hpid, ids)
        right = np.searchsorted(hpid, ids, side='right')
    
        for sid,le,ri in zip(ids,left,right):
            if le != ri:
                if not np.isnan(np.max(skyData['alt'][le:ri])):
                    maxAlt = np.max(skyData['alt'][le:ri])
                    good2 = np.where( skyData['alt'][le:ri] > maxAlt-altPad )
                    brightness = skyData['sky'][le:ri][good2]
                    minSkyMap[sid] = np.percentile(brightness, perClip)
                    maxAltMap[sid] = maxAlt
    
        hp.mollview(minSkyMap, unit='mag/area', title='Faintest Observed Sky Values')
        plotName = os.path.join(plotDir,'minSkyMap.png')
        plt.savefig(plotName)
    
        hp.mollview(np.degrees(maxAltMap), unit='degrees', title='Max Altitude')
        plotName = os.path.join(plotDir,'maxAltMap.png')
        plt.savefig(plotName)
    
def mooninfo(params):

    plotDir = params["moonplotpath"]

    lsstObs = ephem.Observer()
    lsstObs.lat = -0.527868529 #radians of '-30:14:40.7'
    lsstObs.lon = -1.2348102646986 #radians of '-70:44:57.9'
    lsstObs.elevation = 2662.75 #meters
    # Set the target alt and az
    targetAlt = 90.
    targetAz = 0.
    airmass = 1./np.cos(np.radians(90.-targetAlt) )
    esoobj = grabESO(restoreFile='')
    
    # Set to full moon
    esoobj.values['SKYMODEL.MOON.SUN.SEP'] = 180.
    # Middle of the night
    esoobj.values["SKYMODEL.TIME"] = 2
    esoobj.values["SKYMODEL.INCL.ZODIACAL"] = 'N'
    
    moonAlts = np.arange(-60,92,4)
    sbs=[] 
    for alt in moonAlts:
         esoobj.values['SKYMODEL.MOON.ALT'] = alt
         esoobj.values['SKYMODEL.MOON.TARGET.SEP'] =  90.-alt
         sbs.append(esoobj.query(esoobj.values))
    
    esoobj.saveValues()
    
    Rmags = []
    for surf in sbs:
        Rmags.append(surf['R'])
    
    plotName = os.path.join(plotDir,'moonAlt.png')
    plt.plot(moonAlts, Rmags, 'ko')
    plt.title('Full Moon')
    plt.xlabel('Moon Altitude (degrees)')
    plt.ylabel('R-band surface brightness (mag/sq arcsec)')
    plt.savefig(plotName)
    plt.close()

    esoobj = grabESO(restoreFile='')
    # Set to full moon
    esoobj.values['SKYMODEL.MOON.SUN.SEP'] = 180.
    # Middle of the night
    esoobj.values["SKYMODEL.TIME"] = 2
    esoobj.values["SKYMODEL.INCL.ZODIACAL"] = 'N'

    # Set to alt of 45 degrees
    alt = 45.
    esoobj.values['SKYMODEL.MOON.ALT'] = alt
    esoobj.values['SKYMODEL.MOON.TARGET.SEP'] =  90.-alt
    
    phases = np.arange(0,185,5)
    sbs=[] 
    for phase in phases:
         esoobj.values['SKYMODEL.MOON.SUN.SEP'] = phase
         sbs.append(esoobj.query(esoobj.values))
    
    esoobj.saveValues()
    
    Rmags = []
    for surf in sbs:
        Rmags.append(surf['R'])
    
    plotName = os.path.join(plotDir,'moonPhase.png')
    plt.plot(phases, Rmags, 'ko')
    plt.title('Moon Phases, altitude = %i degrees'%alt)
    plt.xlabel('Moon phase (degrees)')
    plt.ylabel('R-band surface brightness (mag/sq arcsec)')
    plt.savefig(plotName)
    plt.close()

def combinezeropointskycatalog(params):

    filters = params["filters"]
    colors = params["colors"]

    for filter,color in zip(filters,colors):

        plotDir = os.path.join(params["zeropointplotpath"],filter)
        if not os.path.isdir(plotDir):
            os.mkdir(plotDir)

        catalogFiles = glob.glob(os.path.join(plotDir,"*.txt"))
        catalogFiles = sorted(catalogFiles)

        catalogFile = os.path.join(params["zeropointplotpath"],'%s.txt'%filter)
        f = open(catalogFile,'w')
        for thiscatalogFile in catalogFiles:
            data_out = np.loadtxt(thiscatalogFile)
            f.write("%.5f %.5f %.5f\n"%(data_out[0],data_out[1],data_out[2]))            
        f.close()

def zeropointskycatalog(params):

    filters = params["filters"]
    colors = params["colors"]

    for filter,color in zip(filters,colors):

        # Read in all the photometry and package it up and save it
        # Set up LSST observatory

        filename = os.path.join(params["ephemplotpath"],"skyData_%s.npz"%filter)
        data_out = np.load(filename)
        skyData = data_out['skyData']

        #indexes = np.where( (skyData['ra'] == params['ra']) & (skyData['dec'] == params['dec']))[0]

        tol = 1e-3
        indexes1 = np.where( (skyData['ra'] > params['ra'] - tol) & (skyData['ra'] < params['ra'] + tol))[0]
        indexes2 = np.where( (skyData['dec'] > params['dec'] - tol) & (skyData['dec'] < params['dec'] + tol))[0]
        indexes = np.intersect1d(indexes1,indexes2)

        mjds = skyData['mjd'][indexes]
        ras = skyData['ra'][indexes]
        decs = skyData['dec'][indexes]
        skys = skyData['sky'][indexes]
        alts = skyData['alt'][indexes]
        sun2targets = skyData['sun2target'][indexes]
        moon2targets = skyData['moon2target'][indexes]
        moonPhases = skyData['moonPhase'][indexes]
        moonAlts = skyData['moonAlt'][indexes]
        sunAlts = skyData['sunAlt'][indexes]
        airmass = 1/np.cos(alts)

        mjds = np.round(mjds)
        mjds_unique = np.unique(mjds)

        plotDir = os.path.join(params["zeropointplotpath"],filter)
        if not os.path.isdir(plotDir):
            os.mkdir(plotDir)

        print len(airmass),len(skys)
     

        #z = np.polyfit(secant, data_all_cut[:,17]+data_all_cut[:,20], 1)
        try:
            z = np.polyfit(airmass,skys, 1)
            p = np.poly1d(z)
            val = p(1)
        except:
            val = np.interp(1,airmass,skys)

        catalogFile = os.path.join(plotDir,"%.5f_%.5f.txt"%(params['ra'],params['dec']))
        f = open(catalogFile,'w')
        f.write("%.5f %.5f %.5f\n"%(params['ra'],params['dec'],val))
        f.close()

        #f = open(catalogFile,'w')
        #for i in xrange(len(mjds_unique)):
        #    indexes = np.where(mjds_unique[i] == mjds)[0]
        #    these_alts = alts[indexes]
        #    these_skys = skys[indexes]
        #    these_airmass = airmass[indexes]

        #    this_airmass = np.interp(airmass_1, these_airmass, these_skys)

        #    f.write("%d %.5f\n"%(mjds_unique[i],this_airmass))
        #f.close()

def ephemskycatalog(params):

    filters = params["filters"]
    colors = params["colors"]

    for filter,color in zip(filters,colors):

        # Read in all the photometry and package it up and save it
        # Set up LSST observatory

        lsstObs = ephem.Observer()
        lsstObs.lat = -0.527868529 #radians of '-30:14:40.7'
        lsstObs.lon = -1.2348102646986 #radians of '-70:44:57.9'
        lsstObs.elevation = 2662.75 #meters
        # lists to store things
        mjd=[]
        ra=[]
        dec=[]
        sky=[]

        # Other things to calc
        alt = []
        moonPhase=[]
        moonAlt=[]
        sunAlt=[]
        moon2target=[]
        sun2target=[]
        
        sunRA = []
        sunDec = []
        moonRA = []
        moonDec= []
       
        catalogDir = os.path.join(params["catalogpath"],filter) 
        files = glob.glob(os.path.join(catalogDir,'*.txt'))
        
        names = ['MJD', 'RA','Dec','HA','xcalc', 'ycalc','BT','VT','V','(B-V)','(V-I)',
                 'm','x','y','HAcalc', 'Decalc','alt','rad','minst','dm','sky','dvig']
        types = [float]*len(names)

        for filename in files:
            try:
                temp_data = np.loadtxt(filename, dtype=zip(names,types))
            except:
                print 'filename %s died'%filename
            mjd.extend(temp_data['MJD'].tolist())
            ra.extend(temp_data['RA'].tolist())
            dec.extend(temp_data['Dec'].tolist())
            sky.extend(temp_data['sky'].tolist())
            alt.extend(temp_data['alt'].tolist())
        
        print 'Read in %i photometry points'%len(mjd)
        
        skyData = np.zeros(len(mjd), dtype=zip(['mjd','ra','dec','sky', 'alt',
                                                'sun2target', 'moon2target', 'moonPhase',
                                                'moonAlt', 'sunAlt'],[float]*10))
        
        tempData = np.zeros(len(mjd), dtype=zip(['sunRA','sunDec','moonRA','moonDec'],[float]*4))
        
        skyData['mjd'] = mjd
        skyData['ra'] = ra
        skyData['dec'] = dec
        skyData['sky'] = sky
        skyData['alt'] = np.radians(alt)

        skyData.sort(order='mjd')
        
        umjd = np.unique(skyData['mjd'])
        left = np.searchsorted(skyData['mjd'], umjd)
        right = np.searchsorted(skyData['mjd'], umjd, side='right')
        
        djd = mjd2djd(skyData['mjd'])
        udjd = np.unique(djd)
        
        print 'looping over %i unique mjd values'%udjd.size
        
        for dj,le,ri in zip(udjd, left,right):
            lsstObs.date = dj
            sun = ephem.Sun(lsstObs)
            moon = ephem.Moon(lsstObs)
            #spot = ephem.FixedBody(ra=r, dec=de)
            #spot.compute(lsstObs)
            #alt.append(spot.alt)
            tempData['sunRA'][le:ri] = sun.ra
            tempData['sunDec'][le:ri] = sun.dec
            tempData['moonRA'][le:ri] = moon.ra
            tempData['moonDec'][le:ri] = moon.dec
            skyData['moonPhase'][le:ri] = moon.phase
            skyData['moonAlt'][le:ri] = moon.alt
            skyData['sunAlt'][le:ri] = sun.alt
        
        skyData['sun2target'] = calcDist_haversine(tempData['sunRA'], tempData['sunDec'], skyData['ra'], skyData['dec'])
        skyData['moon2target'] = calcDist_haversine(tempData['moonRA'],tempData['moonDec'], skyData['ra'], skyData['dec'])

        filename = os.path.join(params["ephemplotpath"],"skyData_%s.npz"%filter)        
        np.savez(filename, skyData=skyData)

def fieldskycatalog(params):

    filters = params["filters"]
    colors = params["colors"]

    for filter,color in zip(filters,colors):

        filename = os.path.join(params["catalogpath"],"%s.txt"%filter)
        mags_data = np.loadtxt(filename)

        plotDir = os.path.join(params["catalogplotpath"],filter)
        if not os.path.isdir(plotDir):
            os.mkdir(plotDir)

        filename = os.path.join(params["catalogplotpath"],"%s.txt"%filter)
        data_out = np.loadtxt(filename)

        ras = data_out[:,0]
        decs = data_out[:,1]
        zz = data_out[:,2]

        indexes = np.where(zz >= 8)
        ras = ras[indexes]
        decs = decs[indexes]
        zz = zz[indexes]

        z_min = 8
        z_max = 9

        plotName = os.path.join(plotDir,'scatter_ra_dec.png')
        plt.scatter(ras,decs,s=50,c=zz, vmin=z_min, vmax=z_max)
        # plot data points.
        #plt.xlim([-180,180])
        #plt.ylim([-90,90])
        plt.xlim([0,60])
        plt.ylim([-80,40])
        plt.xlabel('Right Ascension')
        plt.ylabel('Declination')

        cb = plt.colorbar()
        cb.set_label('Sky brightness')
        #cb.set_label('(minst - m) - m0')
        plt.show()
        plt.savefig(plotName,dpi=200)
        plt.close('all')

        xi = np.linspace(-180,180,3601)
        yi = np.linspace(-90,90,1801)

        xi = np.linspace(-180,180,361)
        yi = np.linspace(-90,90,181)

        XI,YI = np.meshgrid(xi,yi)

        zi = griddata(ras,decs,zz,xi,yi,interp='nn')
        zi = np.array(zi)
        zi[np.isnan(zi)] = 0

        #z_min = np.min(zi)
        #z_max = np.max(zi)

        #contour_levels = np.arange(-2.0, 2.0, 0.1)
        plotName = os.path.join(plotDir,'pcolor_ra_dec.png')
        plt.pcolor(xi, yi, zi, vmin=z_min, vmax=z_max)
        # plot data points.
        plt.xlim([0,60])
        plt.ylim([-80,40])
        plt.xlabel('Right Ascension')
        plt.ylabel('Declination')

        cb = plt.colorbar()
        cb.set_label('Sky brightness')
        #cb.set_label('(minst - m) - m0')
        plt.show()
        plt.savefig(plotName,dpi=200)
        plt.close('all')

        data_out_field = np.loadtxt(params["fieldTable"])
        ras_field = data_out_field[:,2]
        decs_field = data_out_field[:,3]

        filename = os.path.join(params["catalogplotpath"],"%s_field.txt"%filter)
        fid = open(filename,'w')

        vals_field = []

        for i in xrange(len(data_out_field)):

            ra = ras_field[i]
            dec = decs_field[i]

            index1 = np.argmin(np.absolute(xi - ra))
            index2 = np.argmin(np.absolute(yi - dec))

            ra_dec_interp = zi[index2,index1]
            fid.write('%.5f %.5f %.5f\n'%(ra,dec,ra_dec_interp))

            vals_field.append(ra_dec_interp)

        fid.close()

        plotName = os.path.join(plotDir,'scatter_ra_dec_field.png')
        plt.scatter(ras_field,decs_field,s=50,c=vals_field, vmin=z_min, vmax=z_max)
        # plot data points.
        #plt.xlim([-180,180])
        #plt.ylim([-90,90])
        plt.xlim([0,60])
        plt.ylim([-80,40])
        plt.xlabel('Right Ascension')
        plt.ylabel('Declination')

        cb = plt.colorbar()
        cb.set_label('Sky brightness')
        #cb.set_label('(minst - m) - m0')
        plt.show()
        plt.savefig(plotName,dpi=200)
        plt.close('all')

def combineskycatalog(params):

    filters = params["filters"]
    colors = params["colors"]

    for filter,color in zip(filters,colors):

        filename = os.path.join(params["catalogpath"],"%s.txt"%filter)
        mags_data = np.loadtxt(filename)

        plotDir = os.path.join(params["catalogplotpath"],filter)
        if not os.path.isdir(plotDir):
            os.mkdir(plotDir)

        dirpath = os.path.join(params["fitsdisk"],"%s*"%params["folderName"])
        folders = glob.glob(dirpath)
        folders = sorted(folders)

        filename = os.path.join(params["catalogplotpath"],"%s.txt"%filter)
        fid = open(filename,'w')

        for i in xrange(len(mags_data[:,0])):

            ra = mags_data[i,0]
            dec = mags_data[i,1]

            catalogFile = os.path.join(params["catalogpath"],"%s/%.5f_%.5f.txt"%(filter,ra,dec))
            try:
                data_out = np.loadtxt(catalogFile)
            except:
                continue

            if len(data_out) == 0:
                continue

            tt = data_out[:,0]
            m = data_out[:,10+2]
            sky = data_out[:,-2]

            plotName = os.path.join(plotDir,"%.5f_%.5f.png"%(ra,dec))

            plt.figure()
            plt.plot(np.mod(tt,1),sky,'k*')
            plt.xlabel('Time [days]')
            plt.ylabel('Sky brightness')
            plt.xlim([0,0.4])
            plt.ylim([3,15])
            #cb = plt.colorbar()
            #cb.set_label('minst')
            #cb.set_label('minst - m')
            plt.show()
            plt.savefig(plotName,dpi=200)
            plt.close('all')

            median_val = np.median(sky)
            fid.write('%.5f %.5f %.5f\n'%(ra,dec,median_val))
        
        fid.close()

def skycatalogind(params):

    filters = params["filters"]
    colors = params["colors"]

    for filter,color in zip(filters,colors):

        filename = os.path.join(params["catalogpath"],"%s.txt"%filter)
        mags_data = np.loadtxt(filename)

        dirpath = os.path.join(params["fitsdisk"],"%s*"%params["folderName"])
        folders = glob.glob(dirpath)
        folders = sorted(folders)

        #ra = mags_data[0,0]
        #dec = mags_data[0,1]

        index1 = np.where(params["ra"] == mags_data[:,0])
        index2 = np.where(params["dec"] == mags_data[:,1])
        index = np.intersect1d(index1[0],index2[0])

        ra = params["ra"]
        dec = params["dec"]

        catalogFile = os.path.join(params["catalogpath"],"%s/%.5f_%.5f.txt"%(filter,ra,dec))
        f = open(catalogFile,'w')

        for folder in folders:

            outpath = os.path.join(folder,filter)

            files = glob.glob(os.path.join(outpath,"*.long.%s.fish"%filter))
            files = sorted(files)

            #for file in []:
            for file in files:

                folderName = file.split('.')[-4]
                runNumber = int(folderName)

                fitsfile = file.replace("fish","fits")
                fileprefix = file.replace(".fits","")
                outputfile = file.replace("fish","mch")

                if not os.path.isfile(outputfile):
                    continue

                try:
                    data_out = np.loadtxt(outputfile,comments='RA')
                except:
                    print "File %s missing data... continuing.\n"%outputfile
                    continue
                if len(data_out) == 0:
                    continue

                #if not runNumber == 400:
                #    continue

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

                data_all = data_out

                try:
                    a,b = data_all.shape
                except:
                    continue

                indexes = np.where(~np.isnan(data_all[:,17]))[0]
                data_all = data_all[indexes,:]
                magnitudes = data_all[:,10] - data_all[:,17]

                index1 = np.where(ra == data_all[:,0])[0]
                index2 = np.where(dec == data_all[:,1])[0]

                index = index1

                if len(index) == 0:
                    continue

                this_data = data_all[index,:][0]

                str = "%.5f"%utcobs
                for data in this_data:
                    str = "%s %.5e"%(str,data)
                print str

                f.write('%s\n'%(str))

        f.close()

def skycatalog(params):

    filters = params["filters"]
    colors = params["colors"]

    for filter,color in zip(filters,colors):

        filename = os.path.join(params["catalogpath"],"%s.txt"%filter)
        mags_data = np.loadtxt(filename)

        dirpath = os.path.join(params["fitsdisk"],"%s*"%params["folderName"])
        folders = glob.glob(dirpath)
        folders = sorted(folders)

        for i in xrange(len(mags_data[:,0])):

            ra = mags_data[i,0]
            dec = mags_data[i,1]

            catalogFile = os.path.join(params["catalogpath"],"%s/%.5f_%.5f.txt"%(filter,ra,dec))
            f = open(catalogFile,'w')

            for folder in folders:

                outpath = os.path.join(folder,filter)

                files = glob.glob(os.path.join(outpath,"*.long.%s.fish"%filter))
                files = sorted(files)

                #for file in []:
                for file in files:

                    folderName = file.split('.')[-4]
                    runNumber = int(folderName)

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

                    data_all = data_out

                    indexes = np.where(~np.isnan(data_all[:,17]))[0]
                    data_all = data_all[indexes,:]
                    magnitudes = data_all[:,10] - data_all[:,17]
          
                    index1 = np.where(ra == data_all[:,0])[0]
                    index2 = np.where(dec == data_all[:,1])[0]

                    index = index1

                    if len(index) == 0:
                        continue                

                    this_data = data_all[index,:][0]

                    str = "%.5f"%utcobs
                    for data in this_data:
                        str = "%s %.5e"%(str,data)
                    print str

                    f.write('%s\n'%(str))

            f.close()

             
def compareskybrightness(params):

    filters = params["filters"]
    colors = params["colors"]

    for filter,color in zip(filters,colors):

        plotDir = os.path.join(params["skybrightnessplotpath"],filter)
        if not os.path.isdir(plotDir):
            os.mkdir(plotDir)
        moviedir = os.path.join(plotDir,"movie")
        if not os.path.isdir(moviedir):
            os.mkdir(moviedir)

        dirpath = os.path.join(params["dirpath"],"%s09*"%params["folderName"])
        folders = glob.glob(dirpath)
        folders = sorted(folders)

        data_all = []
        data_date = []

        range = np.linspace(0, 1, len(folders))
        colorsall = cm.rainbow(range)

        n = 0
        plotName = os.path.join(plotDir,'skybrightness_daily.png')
        for folder in folders:

            magsDir = os.path.join(folder,"skybrightnessplots/%s"%filter)
            filename = os.path.join(magsDir,"skybrightness.txt")
            if not os.path.isfile(filename):
                continue

            data_out = np.loadtxt(filename)
            if len(data_out) == 0:
                continue

            a,b = data_out.shape

            if len(data_all) == 0:
                data_all = data_out
            else:
                data_all = np.vstack((data_all,data_out))

            folderName = folder.split("/")[-1]
            data_date.append(folderName)

            tt = data_out[:,0]-data_out[0,0]
            tt = tt * 24.0
            color = colorsall[n]
            plt.semilogy(tt,data_out[:,1],'--',color=color)

            n = n + 1

        plt.xlabel('Hours since twilight')
        plt.ylabel('Sky brightness')
        plt.xlim([0,12])
        plt.ylim([1e-1,1e3])
        #cb = plt.colorbar()
        #cb.set_label('minst')
        #cb.set_label('minst - m')
        plt.show()
        plt.savefig(plotName,dpi=200)
        plt.close('all')

        plotName = os.path.join(plotDir,'skybrightness.png')
        plt.plot(data_all[:,0],data_all[:,1],'*')
        #plt.xlim([0,360])
        #plt.ylim([11,15])
        plt.xlabel('')
        plt.ylabel('m - minst')
        #cb = plt.colorbar()
        #cb.set_label('minst')
        #cb.set_label('minst - m')
        plt.show()
        plt.savefig(plotName,dpi=200)
        plt.close('all')

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

def comparecloud(params):

    filters = params["filters"]
    colors = params["colors"]

    for filter,color in zip(filters,colors):

        plotDir = os.path.join(params["cloudsplotpath"],filter)
        if not os.path.isdir(plotDir):
            os.mkdir(plotDir)
        moviedir = os.path.join(plotDir,"movie")
        if not os.path.isdir(moviedir):
            os.mkdir(moviedir)

        dirpath = os.path.join(params["dirpath"],"%s09*"%params["folderName"])
        folders = glob.glob(dirpath)
        folders = sorted(folders)

        data_all = []
        data_date = []

        range = np.linspace(0, 1, len(folders)-10)
        colorsall = cm.rainbow(range)

        n = 0
        plotName = os.path.join(plotDir,'clouds_daily.png')
        for folder in folders:

            magsDir = os.path.join(folder,"cloudsplots/%s/combine"%filter)
            filename = os.path.join(magsDir,"perc.txt")
            if not os.path.isfile(filename):
                continue

            data_out = np.loadtxt(filename)
            if len(data_out) == 0:
                continue

            a,b = data_out.shape

            if len(data_all) == 0:
                data_all = data_out
            else:
                data_all = np.vstack((data_all,data_out))

            folderName = folder.split("/")[-1]
            data_date.append(folderName)

            tt = data_out[:,0]-data_out[0,0]
            tt = tt * 24.0
            color = colorsall[n]
            plt.plot(tt,data_out[:,1],'--',color=color)

            n = n + 1

        plt.xlabel('Hours since twilight')
        plt.ylabel('percentage greater than 1 magnitude')
        plt.xlim([0,12])
        plt.ylim([0,100])
        #cb = plt.colorbar()
        #cb.set_label('minst')
        #cb.set_label('minst - m')
        plt.show()
        plt.savefig(plotName,dpi=200)
        plt.close('all')

        plotName = os.path.join(plotDir,'clouds.png')
        plt.plot(data_all[:,0],data_all[:,1],'*')
        #plt.xlim([0,360])
        #plt.ylim([11,15])
        plt.xlabel('')
        plt.ylabel('m - minst')
        #cb = plt.colorbar()
        #cb.set_label('minst')
        #cb.set_label('minst - m')
        plt.show()
        plt.savefig(plotName,dpi=200)
        plt.close('all')

        n=1
        for folder in folders:

            file = os.path.join(folder,"cloudsplots/%s/combine/tf.png"%filter)
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

        print plotDir

def initialize(params):

    # set up directory names and data paths
    params["datadisk"]='/lsst/all-sky'
    params["dirpath"]='/lsst/home/coughlin/allsky/data'
    params["fitsdisk"]='/lsst/home/coughlin/allsky/data/FITS'
    params["catalogpath"] = "/lsst/home/coughlin/git-repo/fisheye/catalog"
    params["fieldTable"] = "/home/coughlin/git-repo/fisheye/fields/Field_table.txt"
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

    params["catalogplotpath"]=os.path.join(params["dirpathname"],"catalogplots")
    if not os.path.isdir(params["catalogplotpath"]):
        os.mkdir(params["catalogplotpath"])

    params["fieldplotpath"]=os.path.join(params["dirpathname"],"fieldplots")
    if not os.path.isdir(params["fieldplotpath"]):
        os.mkdir(params["fieldplotpath"])

    params["ephemplotpath"]=os.path.join(params["dirpathname"],"ephemplots")
    if not os.path.isdir(params["ephemplotpath"]):
        os.mkdir(params["ephemplotpath"])

    params["zeropointplotpath"]=os.path.join(params["dirpathname"],"zeropointplots")
    if not os.path.isdir(params["zeropointplotpath"]):
        os.mkdir(params["zeropointplotpath"])

    params["analyzeephemplotpath"]=os.path.join(params["dirpathname"],"analyzeephemplots")
    if not os.path.isdir(params["analyzeephemplotpath"]):
        os.mkdir(params["analyzeephemplotpath"])

    params["moonplotpath"]=os.path.join(params["dirpathname"],"moonplots")
    if not os.path.isdir(params["moonplotpath"]):
        os.mkdir(params["moonplotpath"])

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
    ra = opts.ra
    dec = opts.dec
    params = {}

    params["outputFolder"] = outputFolder
    params["folderName"] = folderName
    params["ra"] = ra
    params["dec"] = dec

    params = initialize(params)

    # Compare fisheye fits across nights
    params["doCompareCloud"] = opts.doCompareCloud
    # Compare sky brightness across nights
    params["doCompareSkyBrightness"] = opts.doCompareSkyBrightness
    # Generate sky catalog
    params["doSkyCatalog"] = opts.doSkyCatalog
    # Generate sky catalog individual
    params["doSkyCatalogInd"] = opts.doSkyCatalogInd
    # Combine sky catalog
    params["doCombineSkyCatalog"] = opts.doCombineSkyCatalog
    # LSST field / sky catalog
    params["doFieldSkyCatalog"] = opts.doFieldSkyCatalog
    # LSST ephemeris field / sky catalog
    params["doEphemSkyCatalog"] = opts.doEphemSkyCatalog
    # Zeropoint LSST ephemeris field / sky catalog
    params["doZeropointSkyCatalog"] = opts.doZeropointSkyCatalog
    # Combine zeropoint LSST ephemeris field / sky catalog
    params["doCombineZeropointSkyCatalog"] = opts.doCombineZeropointSkyCatalog
    # Analyze LSST ephemeris field / sky catalog
    params["doAnalyzeEphemSkyCatalog"] = opts.doAnalyzeEphemSkyCatalog 
    # Moon info plots
    params["doMoonInfo"] = opts.doMoonInfo

    if params["doCompareSkyBrightness"]:
        compareskybrightness(params)
    if params["doCompareCloud"]:
        comparecloud(params)
    if params["doSkyCatalog"]:
        skycatalog(params)
    if params["doSkyCatalogInd"]:
        skycatalogind(params)
    if params["doCombineSkyCatalog"]:
        combineskycatalog(params)
    if params["doFieldSkyCatalog"]:
        fieldskycatalog(params)
    if params["doEphemSkyCatalog"]:
        ephemskycatalog(params)
    if params["doZeropointSkyCatalog"]:
        zeropointskycatalog(params)
    if params["doCombineZeropointSkyCatalog"]:
        combinezeropointskycatalog(params)
    if params["doAnalyzeEphemSkyCatalog"]:
        analyzeephemskycatalog(params)
    if params["doMoonInfo"]:
        mooninfo(params)
