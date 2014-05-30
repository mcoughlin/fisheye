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

    parser.add_option("-m", "--maxframes", help="max frames",default=500,type=int)
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

    opts, args = parser.parse_args()

    return opts

def clouds(params):

    filters = params["filters"]
    colors = params["colors"]

    for filter in filters:
        outpath = os.path.join(params["dirpathname"],filter)

        baseplotDir = os.path.join(params["cloudsplotpath"],filter)
        if not os.path.isdir(baseplotDir):
            os.mkdir(baseplotDir)

        magsDir = os.path.join(params["fisheyeplotpath"],filter)
        filename = os.path.join(magsDir,"mags.txt")
        mags_data = np.loadtxt(filename)

        files = glob.glob(os.path.join(outpath,"*.long.%s.fish"%filter))
        files = sorted(files)

        zi_ave = []

        #for file in []:
        for file in files:

            folderName = file.split('.')[-4]
            runNumber = int(folderName)
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
            plt.scatter(data_all[:,3],data_all[:,4],s=10,c=magnitudes,vmin=12,vmax=15)
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
            for i in xrange(len(data_all[:,0])):

                ra = data_all[i,0]
                ha = LST - ra
                index = np.where(ra == mags_data[:,0])[0]
                if len(mags_data[index,2]) == 0:
                    mags.append(magnitudes[i])
                else:
                    mags.append(mags_data[index,2][0])


            data_all_new = []
            #mags = []
            #for i in xrange(len(data_all[:,0])):
            for i in xrange(len(mags_data[:,0])):

                continue

                ra = mags_data[i,0]
                dec = mags_data[i,1]
                ha = LST - ra
                index = np.where(ra == data_all[:,0])[0]

                if len(index) > 1:
                    index = [index[0]]

                mags.append(mags_data[i,2])

                if len(mags_data[index,2]) == 0:
 
                    fisheye_command = "echo '%.5f %.5f' | fisheye -ha %.5f -dec %.5f -az %.5f -scale %.5f -quad %.5f -cube %.5f -cx %.5f -cy %.5f -lat %.5f -sky2xy"%(ha,dec,HA,DEC,AZ,SCALE,QUAD,CUBE,CX,CY,LAT)            
                    output = os.popen(fisheye_command).readlines()[0].strip()
                    outputSplit = output.split(" ")
                    outputSplit = [x for x in outputSplit if x != '']

                    this_data = data_all[0,:]
                    this_data[0] = ra
                    this_data[1] = dec
                    this_data[2] = ha
                    this_data[3] = float(outputSplit[0])
                    this_data[4] = float(outputSplit[1])
                    this_data[10] = mags_data[i,2]
                    this_data[17] = -2

                    if len(data_all_new) == 0:
                        data_all_new = np.array(this_data)
                    else:
                        data_all_new = np.vstack((data_all_new,this_data))
                else:
                    if len(data_all_new) == 0:
                        data_all_new = np.array(data_all[index,:])
                    else:
                        data_all_new = np.vstack((data_all_new,data_all[index,:]))

            #data_all = data_all_new
            magnitudes = data_all[:,10] - data_all[:,17]
            mags = np.array(mags).T
            diffs = magnitudes - mags

            if np.sum(diffs) == 0.0:
                continue

            if len(diffs) < 3:
                continue

            plotName = os.path.join(plotDir,'diff_x_y.png')
            plt.scatter(data_all[:,3],data_all[:,4],s=10,c=diffs,vmin=-1,vmax=1)
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

            xi = np.linspace(0,2897,1000)
            yi = np.linspace(0,1935,1000)
            zi = griddata(data_all[:,3],data_all[:,4],diffs,xi,yi,interp='nn')
            zi = np.array(zi)

            zi_sum_this = zi
            zi_sum_this[np.isnan(zi)] = 0
            zi_sum_this[np.abs(zi)>0] = 1

            if len(zi_ave) == 0:
                zi_ave = zi
                #zi_sum = zi_sum_this
            else:
                zi_ave = zi_ave + zi
                #zi_sum = zi_sum + zi_sum_this

            contour_levels = np.arange(-3, 3, 0.2)
            plotName = os.path.join(plotDir,'contour_x_y.png')
            plt.scatter(data_all[:,3],data_all[:,4],s=10,c=diffs,zorder=10,vmax=3, vmin=-3)
            cs = plt.contour(xi,yi,zi,contour_levels,linewidths=0.5,colors='k')
            plt.contourf(xi,yi,zi,contour_levels,cmap=plt.cm.rainbow)
            cb = plt.colorbar() # draw colorbar
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

            #p = cs.collections[0].get_paths()[0]
            #v = p.vertices
            #x = v[:,0]
            #y = v[:,1]

        zi_ave_norm = np.array(zi_ave) / float(len(files))

        plotDir = os.path.join(baseplotDir,"ave")
        if not os.path.isdir(plotDir):
            os.mkdir(plotDir)

        contour_levels = np.arange(-2, 2, 0.2)
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

def combinefisheye(params):

    filters = params["filters"]
    colors = params["colors"]

    for filter in filters:
        outpath = os.path.join(params["dirpathname"],filter)

        data_all = []

        files = glob.glob(os.path.join(outpath,"*.long.%s.fish"%filter))
        files = sorted(files)

        runNumbers = []
        nstars = []

        for file in files:

            runNumber = int(file.split('.')[-4])
            if runNumber > params["maxframes"]:
                continue

            fileprefix = file.replace(".fits","")
            outputfile = file.replace("fish","mch")

            if not os.path.isfile(outputfile):
                continue
            data_out = np.loadtxt(outputfile,comments='RA')

            runNumbers.append(runNumber)
            nstars.append(len(data_out))

            if len(data_out) == 0:
                continue            

            if len(data_all) == 0:
                data_all = data_out
            else:
                data_all = np.vstack((data_all,data_out))

            #if len(data_out) < 2000:
            #    print file
            #    print len(data_out)

        plotDir = os.path.join(params["fisheyeplotpath"],filter)
        if not os.path.isdir(plotDir):
            os.mkdir(plotDir)

        if len(data_all) == 0:
            continue

        plotName = os.path.join(plotDir,'nstars.png')
        plt.plot(runNumbers,nstars,'*')
        plt.xlabel('run number')
        plt.ylabel('number of stars')
        plt.show()
        plt.savefig(plotName,dpi=200)
        plt.close('all')

        indexes = np.where(~np.isnan(data_all[:,17]))[0]
        data_all = data_all[indexes,:]
        magnitudes = data_all[:,10] - data_all[:,17]

        indexes = np.where(data_all[:,18] < 0.1)[0]
        data_all = data_all[indexes,:]
        magnitudes = data_all[:,10] - data_all[:,17]

        zenith = (90 - data_all[:,15]) * 2 * np.pi / 360.0
        secant = (1.002432 * np.cos(zenith)**2 + 0.148386 * np.cos(zenith) + 0.0096467) / (np.cos(zenith)**3 + 0.149864 * np.cos(zenith)**2 + 0.0102963 * np.cos(zenith) + 0.00303978)

        plotName = os.path.join(plotDir,'mag_x_y.png')
        plt.scatter(data_all[:,3],data_all[:,4],s=10,c=magnitudes,vmin=12,vmax=15)
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

        data_all = data_all[data_all[:,0].argsort(),:]
        ar,return_index,return_inverse = np.unique(data_all[:,0],return_index=True, return_inverse=True)
        indexes = np.where(return_inverse == mode(return_inverse)[0])[0]
        data_all_cut = data_all[indexes,:]
        magnitudes_cut = magnitudes[indexes]

        plotName = os.path.join(plotDir,'mag_x_y_max.png')
        #plt.scatter(data_all_cut[:,3],data_all_cut[:,4],s=10,c=data_all_cut[:,17]+data_all_cut[:,20])
        plt.scatter(data_all_cut[:,3],data_all_cut[:,4],s=10,c=magnitudes_cut)
        plt.xlabel('x')
        plt.ylabel('y')
        cb = plt.colorbar()
        cb.set_label('minst - m')
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

        zenith = (90 - data_all_cut[:,15]) * 2 * np.pi / 360.0
        secant = (1.002432 * np.cos(zenith)**2 + 0.148386 * np.cos(zenith) + 0.0096467) / (np.cos(zenith)**3 + 0.149864 * np.cos(zenith)**2 + 0.0102963 * np.cos(zenith) + 0.00303978) 

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

            zenith = (90 - data_all_cut[:,15]) * 2 * np.pi / 360.0
            secant = (1.002432 * np.cos(zenith)**2 + 0.148386 * np.cos(zenith) + 0.0096467) / (np.cos(zenith)**3 + 0.149864 * np.cos(zenith)**2 + 0.0102963 * np.cos(zenith) + 0.00303978)

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

            #plt.plot(secant,magnitudes_cut,'*')
            #plt.plot(xp, p(xp),'b-')

            #print magnitudes_cut
            #print p(magnitudes_cut)
            rmserror = np.mean(np.absolute(magnitudes_cut - p(secant)))
            #if rmserror > 1.5:
            #    print data_all_cut[0,0],data_all_cut[0,1],p
            #    print magnitudes_cut

            p = np.median(magnitudes_cut)

            rmserrors.append(rmserror)
            plt.plot(data_all_cut[0,0],rmserror,'*')

            fid.write("%.5f %.5f %.10e\n"%(data_all_cut[0,0],data_all_cut[0,1],p))

        plt.xlabel('secant')
        plt.ylabel('minst')
        #plt.xlim([1,3])
        plt.show()
        plt.savefig(plotName,dpi=200)
        plt.close('all')

        #print stop

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

            fileprefix = file.replace(".fits","")

            fisheye_command = "./imstats_fish.sh %s"%fileprefix
            os.system(fisheye_command)

            outputfile = file.replace("fits","fish")

            if not os.path.isfile(outputfile):
                continue

            data_out = np.loadtxt(outputfile)

            if len(data_out) == 0:
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

def initialize(params):

    # set up directory names and data paths
    params["datadisk"]='/lsst/all-sky'
    params["fitsdisk"]='/lsst/all-sky/FITS'    
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
    params["outputFolder"] = outputFolder
    params["maxframes"] = maxframes
    params["user"] = user

    # Create fits files from .gz files
    params["doMakeFits"] = opts.doMakeFits
    # Do analysis where we look for clouds
    params["doClouds"] = opts.doClouds
    # Run fisheye code, which converts from pixel x,y to RA and Dec
    params["doFisheye"] = opts.doFisheye
    # Track stars across the night
    params["doCombineFisheye"] = opts.doCombineFisheye
    # Create movie of fits images
    params["doMakeMovie"] = opts.doMakeMovie
    # Run tphot / source extractor on the images
    params["doPhotometry"] = opts.doPhotometry

    params = initialize(params)

    if params["doMakeFits"]:
        makefits(params)
    if params["doMakeMovie"]:
        makemovie(params)
    if params["doPhotometry"]:
        photometry(params)
        photometry_tphot(params)
    if params["doFisheye"]:
        fisheye(params)
    if params["doCombineFisheye"]:
        combinefisheye(params)
    if params["doClouds"]:
        clouds(params)

