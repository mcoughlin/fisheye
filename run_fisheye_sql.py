import numpy as np
import os,glob, sys
import ephem
#from lsst.sims.maf.utils.telescopeInfo import TelescopeInfo
#from lsst.sims.skybrightness.utils import mjd2djd
from telescopeInfo import TelescopeInfo
from skybrightness_utils import mjd2djd


def create_sql(filename,catalognightDir,photodiodefiles,skybrightnessfiles):

    f = open(filename,'w')

    f.write("""
CREATE TABLE stars(
   ID INT PRIMARY KEY     NOT NULL,
   ra             REAL,
   dec            REAL,
   catalogmag     REAL
);

CREATE TABLE dates(
    ID INT PRIMARY KEY  NOT NULL,
    mjd      REAL,
    sunAlt   REAL,
    moonAlt  REAL,
    moonPhase REAL
);

CREATE TABLE obs(
  ID  INT PRIMARY KEY  NOT NULL,
  starID     INT,
  dateID     INT,
  x          REAL,
  y          REAL,
  alt        REAL,
  starMag    REAL,
  starMag_inst REAL,
  starMag_err REAL,
  sky         REAL,
  filter     TEXT
);

CREATE TABLE photdiode(
 mjd    REAL,
 R      REAL,
 Y      REAL,
 Z      REAL
);

CREATE TABLE skybrightness(
 mjd    REAL,
 M      REAL,
 R      REAL,
 G      REAL,
 B      REAL
);\n\n
""")

    f.write(
"""
.separator ","
.import %s/obsTable.dat obs
.import %s/mjdTable.dat dates
.import %s/starTable.dat stars
    
CREATE INDEX obsStarID on obs (starID);
CREATE INDEX obsDateID on obs (dateID);

.separator " "
"""%(catalognightDir,catalognightDir,catalognightDir))


    for photodiodefile in photodiodefiles:
        f.write(".import %s photdiode\n"%photodiodefile)

    f.write("\n\nCREATE INDEX diodeDateID on photdiode (mjd);\n")

    for skybrightnessfile in skybrightnessfiles:
        f.write(".import %s skybrightness\n"%skybrightnessfile)

    f.write("\n\nCREATE INDEX skyDateID on skybrightness (mjd);\n")

    f.close()

catalognightDir = '/home/coughlin/git-repo/fisheye/catalognight' 

FITSfolder = "/scratch/coughlin/allsky/data"
folders = glob.glob(os.path.join(FITSfolder,"ut*"))
folders = sorted(folders)
photodiodefiles = []
for folder in folders:
    filename = os.path.join(folder,'photodiodeplots/photodiode.txt')
    if os.path.isfile(filename):
        photodiodefiles.append(filename)
photodiodefiles = sorted(photodiodefiles)

skybrightnessfiles = []
for folder in folders:
    filename = os.path.join(folder,'skybrightnessplots/skybrightness.txt')
    if os.path.isfile(filename):
        skybrightnessfiles.append(filename)

filename = os.path.join(catalognightDir,'createDatabase.sql')
create_sql(filename,catalognightDir,photodiodefiles,skybrightnessfiles)

print stop

# Set up LSST telescope
telescope = TelescopeInfo('LSST')
Observatory = ephem.Observer()
Observatory.lat = telescope.lat
Observatory.lon = telescope.lon
Observatory.elevation = telescope.elev

sun = ephem.Sun()
moon = ephem.Moon()

# Need to read through all the files to make the starID's
bfiles = glob.glob('%s/B/*.txt'%catalognightDir)
rfiles = glob.glob('%s/R/*.txt'%catalognightDir)
gfiles = glob.glob('%s/G/*.txt'%catalognightDir)
mfiles = glob.glob('%s/M/*.txt'%catalognightDir)

allFiles = bfiles + rfiles + gfiles + mfiles

ra=np.zeros(0, dtype=float)
dec=np.zeros(0, dtype=float)
radec=np.zeros(0, dtype=float)
names = ['mjd', 'ra','dec']
types = [float]*3
mjds = np.zeros(0, dtype=float)

for filename in allFiles:
    data = np.loadtxt(filename, dtype=zip(names,types))
    tempRaDec = np.round(data['dec']*1000)*10+data['ra']
    tempRaDec, uind = np.unique(tempRaDec, return_index=True)

    radec = np.append(radec, tempRaDec)
    ra = np.append(ra, data['ra'][uind])
    dec = np.append(dec, data['dec'][uind])

    radec, uind = np.unique(radec, return_index=True)
    ra = ra[uind]
    dec = dec[uind]
    mjds = np.append(mjds, data['mjd'])
    mjds = np.unique(mjds)


# Generate starID table
starids = np.arange(ra.size)
f = open('%s/starTable.dat'%catalognightDir, 'w')
for i,rai in enumerate(ra):
    print >>f, '%i,%f,%f,0' % (i,ra[i],dec[i])
f.close()

# Generate mjd table
f = open('%s/mjdTable.dat'%catalognightDir, 'w')
mjdID = np.arange(mjds.size)
for mjdid,mjd in zip(mjdID,mjds):
    Observatory.date = mjd2djd(mjd)
    sun.compute(Observatory)
    moon.compute(Observatory)
    print >>f, '%i,%f,%f,%f,%f' % (mjdid, mjd, sun.alt, moon.alt, moon.phase)
f.close()


# now to loop through and write the obs table and the mjd table
names = ['mjd', 'ra','dec', 'm', 'x', 'y', 'alt', 'minst', 'dm', 'sky', 'band']
types = [float]*(len(names)-1)
types.append('|S1')

obsidMax = 0

f = open('%s/obsTable.dat'%catalognightDir, 'w')
maxJ = float(len(allFiles))

for j,filename in enumerate(allFiles):
    # Maybe read in a dummy column, set it to the filter and then stack all of these so they can quickly be sorted?
    data = np.genfromtxt(filename, dtype = zip(names,types),
                         usecols=(0,1,2,11,12,13,16,18,19,20,21))
    if data.size > 0:
        filenameSplit = filename.split("/")
        thisfilter = filenameSplit[-2]
        data['band'] = thisfilter

        data.sort(order=['mjd'])
        # Look up starIDs and mjdIDs
        starIDs = np.zeros(data.size, dtype=int)
        mjdIDs = np.zeros(data.size, dtype=int)
        obsIDs = np.arange(data.size) + obsidMax
        obsidMax = obsIDs.max()+1

        left = np.searchsorted(mjds, data['mjd']  )
        right = np.searchsorted(mjds, data['mjd'], side='right')

        for i,ack in enumerate(left):
            mjdIDs[i] = mjdID[left[i]:right[i]]

        tempRaDec =  np.round(data['dec']*1000)*10+data['ra']
        ord = np.argsort(tempRaDec)
        tempRaDec = tempRaDec[ord]
        data = data[ord]
        mjdIDs = mjdIDs[ord]


        left = np.searchsorted(radec, tempRaDec)
        right = np.searchsorted(radec, tempRaDec, side='right')
        for i,ack in enumerate(left):
            starIDs[i] = starids[left[i]:right[i]]

        for i,starid in enumerate(starIDs):
            print >>f, '%i,%i,%i,%f,%f,%f,%f,%f,%f,%f,%s' % (obsIDs[i], starIDs[i], mjdIDs[i], 
                                                    data['x'][i], data['y'][i], data['alt'][i],
                                                    data['m'][i], data['minst'][i], data['dm'][i], 
                                                    data['sky'][i],data['band'][i])
        progress = j/maxJ*100
        text = "\rprogress = %.1f%%"%progress
        sys.stdout.write(text)
        sys.stdout.flush()
f.close()









