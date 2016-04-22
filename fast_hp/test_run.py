import glob
import numpy as np
import healpy as hp
from lsst.sims.utils import _healbin, _hpid2RaDec
from lsst.sims.skybrightness import stupidFast_RaDec2AltAz, stupidFast_altAz2RaDec
from lsst.sims.utils import Site
from astropy.io import fits
import os
import subprocess
from shutil import copyfile
import matplotlib.pylab as plt



if __name__ == '__main__':
	

	altaz = np.load('coordims.npz')
	alt = altaz['alt'].copy()
	az = altaz['az'].copy()
	altaz.close()

	# Crop down to only the pixels with alt > 0
	above_ground = np.where(alt > 0.)

	site = Site('LSST')

	output_dir = 'output'

	nside = 32
	hpindex = np.arange(hp.nside2npix(nside))
	hpra, hpdec = _hpid2RaDec(nside, hpindex)
	filters = {'R':0, 'G':1,'B':2}
	night_dirs = ['ut012716']#['ut012916']
	# Don't bother with frames that are above this median counts
	countLimit = 13000.

	for night_dir in night_dirs:
		zipped_files = glob.glob(os.path.join(night_dir, '*.long*.gz'))
		if not os.path.exists(os.path.join(output_dir, night_dir)):
			os.makedirs(os.path.join(output_dir, night_dir))
		outfile = open(os.path.join(output_dir, night_dir, 'healmaps.dat'), 'w')
		# I hate to leave off the header, but it makes it easier to slurp into sqlite
		# print >>outfile, '# hpid (nside=%i), R, G, B, airmass, mjd' % nside
		for filename in zipped_files:
			# copy and unzip the file
			tempfile = os.path.split(filename)[-1]
			copyfile(filename, tempfile)
			subprocess.call('gunzip -f %s' % tempfile, shell=True)
			crfile = tempfile[:-3]
			if os.path.isfile(crfile):
				subprocess.call('raw2fits -rgb %s' % crfile, shell=True)
				fitsfile = crfile[:-3] + 'fits'
				# read in the fits, assuming it's rgb order?
				hdulist = fits.open(fitsfile)
				mjd = hdulist[0].header['MJD-OBS']
				ra, dec = stupidFast_altAz2RaDec(alt[above_ground], az[above_ground], site.latitude_rad, site.longitude_rad,  mjd)
				hpAlt, hpAz = stupidFast_RaDec2AltAz(hpra, hpdec, site.latitude_rad, site.longitude_rad, mjd)
				airmass = 1./np.cos(np.pi/2. - hpAlt)
				hpAbove = np.where(hpAlt >= 0.)
				hpmaps = {}
				for filtername in filters.keys():
					image = hdulist[0].data[filters[filtername],:,:] - hdulist[0].header['bias']
					# Crop to only images above ground.
					#image = image.ravel()[above_ground]
					# Put a check that this isn't a very saturated frame
					if np.median(image)+hdulist[0].header['bias'] < countLimit:
						hpmap = _healbin(ra, dec, image[above_ground], nside=nside, reduceFunc=np.median)
					else:
						hpmap = np.empty(hp.nside2npix(nside))
						hpmap.fill(hp.UNSEEN)
					hpmaps[filtername] = hpmap
				maskcheck = np.unique(np.array([np.median(val[hpAbove]) for val in hpmaps.values()]))
				if np.max(maskcheck) != hp.UNSEEN:
					good = np.where(hpAlt >= 0.)[0]
					for indx in good:
						# Skip if everythin in the pixel is masked
						if not ((hpmaps['R'][indx] == hp.UNSEEN) & (hpmaps['G'][indx] == hp.UNSEEN) & (hpmaps['B'][indx] == hp.UNSEEN)): 
							print >>outfile, '%i, %f, %f, %f, %f, %f' % (hpindex[indx], hpmaps['R'][indx], hpmaps['G'][indx],
						                                             	 hpmaps['B'][indx],airmass[indx], mjd)
				# delete the temp files that were created
				hdulist.close()
				os.unlink(crfile)
				os.unlink(fitsfile)
			else:
				os.unlink(tempfile)
			

		outfile.close()


