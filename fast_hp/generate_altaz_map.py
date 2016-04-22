import numpy as np
import ephem
from lsst.sims.skybrightness import stupidFast_RaDec2AltAz
from lsst.sims.utils import Site
from astropy.io import fits




def make_altaz():
    """
    Load up the alt - az coordinates from a single frame
    """

    # Coordinates are hour angle and declination.
    coords = np.genfromtxt('coord_im/ut012716.0100.long.M.xxyy',
                           dtype=[int,int,float,float])

    hdulist = fits.open('coord_im/ut012716.0100.long.M.fits')

    image = hdulist[0].data
    xim = np.zeros(image.shape, dtype=int)-666
    yim = np.zeros(image.shape, dtype=int)-666
    haim = np.zeros(image.shape, dtype=float)-666
    decim = np.zeros(image.shape, dtype=float)-666

    for coord in coords:
        xim[coord['f1'], coord['f0']] = coord['f0']
        yim[coord['f1'], coord['f0']] = coord['f1']
        haim[coord['f1'], coord['f0']] = coord['f2']
        decim[coord['f1'], coord['f0']] = coord['f3']


    mjd = hdulist[0].header['MJD-OBS']
    obs = ephem.Observer()
    site = Site('LSST')
    obs.lon = site.longitude_rad
    obs.lat = site.latitude_rad
    obs.elevation = site.height
    doff = ephem.Date(0)-ephem.Date('1858/11/17')
    obs.date = mjd - doff
    lst = obs.sidereal_time()

    ra = lst - np.radians(haim) 
    while ra.min() < 0:
        ra += 2.*np.pi
    ra = ra % (2.*np.pi)

    alt, az = stupidFast_RaDec2AltAz(ra, np.radians(decim), site.latitude_rad, site.longitude_rad, mjd)
    mask = np.where(xim == -666)
    alt[mask] = -666
    az[mask] = -666


    np.savez('coordims.npz', alt=alt, az=az)  #xim=xim,yim=yim,haim=haim,decim=decim, ra=ra, alt=alt, az=az)


if __name__ == '__main__':

    make_altaz()
    data = np.load('coordims.npz')

    