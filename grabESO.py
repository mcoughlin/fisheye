import numpy as np
from collections import OrderedDict
import sys, os
import urllib
import urllib2
import warnings

# make a python function that can grab the sky brightness values from the web form.


class grabESO(object):

    def __init__(self, restoreFile='callArchive.npy', values=None):
        """

        """
        self.restoreFile = restoreFile
        #url = 'https://www.eso.org/observing/etc/bin/gen/form?INS.MODE=swspectr+INS.NAME=SKYCALC'
        self.url='https://www.eso.org/observing/etc/bin/simu/skycalc'
        self.keys = ['INS.NAME',"INS.MODE",'SKYMODEL.TARGET.ALT','SKYMODEL.TARGET.AIRMASS',"SKYMODEL.SEASON",
                "SKYMODEL.TIME","SKYMODEL.PWV.MODE","SKYMODEL.MSOLFLUX","SKYMODEL.INCL.MOON",
                "SKYMODEL.MOON.SUN.SEP","SKYMODEL.MOON.TARGET.SEP","SKYMODEL.MOON.ALT",
                "SKYMODEL.MOON.EARTH.DIST","SKYMODEL.INCL.STARLIGHT","SKYMODEL.INCL.ZODIACAL",
                "SKYMODEL.ECL.LON","SKYMODEL.ECL.LAT","SKYMODEL.INCL.MOLEC.EMIS.LOWER.ATM",
                "SKYMODEL.INCL.MOLEC.EMIS.UPPER.ATM","SKYMODEL.INCL.AIRGLOW","SKYMODEL.INCL.THERMAL",
                "SKYMODEL.WAVELENGTH.MIN","SKYMODEL.WAVELENGTH.MAX","SKYMODEL.WAVELENGTH.GRID.MODE",
                "SKYMODEL.WAVELENGTH.RESOLUTION","SKYMODEL.LSF.KERNEL.TYPE","SKYCALC.RAD.PLOT.FLAG",
                "SKYCALC.TRANS.PLOT.FLAG","SKYCALC.MAG.FLAG","SKYCALC.LSF.PLOT.FLAG"]
        if values is None:
            self.values =  {'INS.NAME':"SKYCALC",
                      "INS.MODE":"swspectr",
                      'SKYMODEL.TARGET.ALT': 90.00,
                      'SKYMODEL.TARGET.AIRMASS': 1.00,
                      "SKYMODEL.SEASON": 0,
                      "SKYMODEL.TIME": 0,
                      "SKYMODEL.PWV.MODE": 2.5,
                      "SKYMODEL.MSOLFLUX": 130.00,
                      "SKYMODEL.INCL.MOON": "Y",
                      "SKYMODEL.MOON.SUN.SEP": 90., #0-180
                      "SKYMODEL.MOON.TARGET.SEP": 45.00, #0-180
                      "SKYMODEL.MOON.ALT":45.00, #-90-90
                      "SKYMODEL.MOON.EARTH.DIST": 1., #[0.945, 1.055], mean=1
                      "SKYMODEL.INCL.STARLIGHT": "Y",
                      "SKYMODEL.INCL.ZODIACAL": "Y",
                      "SKYMODEL.ECL.LON": 135., #[-180, 180]
                      "SKYMODEL.ECL.LAT": 90., #[-180, 180]
                      "SKYMODEL.INCL.MOLEC.EMIS.LOWER.ATM": "Y",
                      "SKYMODEL.INCL.MOLEC.EMIS.UPPER.ATM": "Y",
                      "SKYMODEL.INCL.AIRGLOW": "Y",
                      "SKYMODEL.INCL.THERMAL": "N",
                      "SKYMODEL.WAVELENGTH.MIN": 1000.00,
                      "SKYMODEL.WAVELENGTH.MAX": 2000.00,
                      "SKYMODEL.WAVELENGTH.GRID.MODE": "fixed_spectral_resolution",
                      "SKYMODEL.WAVELENGTH.RESOLUTION": 20000,
                      "SKYMODEL.LSF.KERNEL.TYPE": "none",
                      "SKYCALC.RAD.PLOT.FLAG": 0,
                      "SKYCALC.TRANS.PLOT.FLAG": 0,
                      "SKYCALC.MAG.FLAG": 1,
                      "SKYCALC.LSF.PLOT.FLAG":0}

        # If the restore file exists, restore it
        if os.path.isfile(restoreFile):
            self.oldCalls = np.load(restoreFile)[()]
        else:
            self.oldCalls={}

    def vals2key(self, values):
        thisKey = ''
        for key in self.keys:
            thisKey = thisKey + str(values[key])
        return thisKey

    def checkOld(self, values):
        """
        Check if this query has already been made
        """
        key = self.vals2key(self.values)
        if key in self.oldCalls.keys():
            result = self.oldCalls[key]
        else:
            result = None
        return result

    def query(self, values):
        """
        values: Dict with all the parameters for the web page
        """

        prevResult = self.checkOld(values)
        if prevResult is not None:
            print 'Restoring previous query'
            return prevResult
        else:
            print 'Making new query'
            data = urllib.urlencode(self.values)
            req = urllib2.Request(self.url, data)
            response = urllib2.urlopen(req)
            the_page = response.read().split('\n')

            # Now need to read through the output and pull out the surface brightnesses
            sbs = {}
            filters = ['U','B','V','R','I','Z','Y','J','H','K','L','M','N','Q']
            filters = [ filt+':' for filt in filters]
            for line in the_page:
                line = line.lstrip().split(' ')
                line = [x for x in line if x != '']
                if len(line) > 0:
                    if line[0] in filters:
                        sbs[line[0][0] ] = float(line[1])
            if sbs == {}:
                warnings.warn('Did not get a result')
            newKey = self.vals2key(values)
            self.oldCalls[newKey] = sbs
            return sbs

    def saveValues(self):
        """
        Save any new results.
        """
        np.save(self.restoreFile, self.oldCalls)


