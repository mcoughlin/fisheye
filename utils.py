
import numpy as np
import healpy as hp
import ephem

def mjd2djd(inDate):
    """
    Convert Modified Julian Date to Dublin Julian Date (what pyephem uses).
    """
    if not hasattr(mjd2djd, 'doff'):
        mjd2djd.doff = ephem.Date(0)-ephem.Date('1858/11/17')
    djd = inDate-mjd2djd.doff
    return djd

def radec2pix(nside, ra, dec):
    """Convert ra,dec to the nearest heaplixel id."""
    lat = np.pi/2. - dec
    hpid = hp.ang2pix(nside, lat, ra )
    return hpid

def calcDist_haversine(RA1, Dec1, RA2, Dec2):
    """Calculates distance on a sphere using the haversine formula.
    Give this function RA/Dec values in radians. Returns angular distance(s), in radians.
    Note that since this is all numpy, you could input arrays of RA/Decs."""
    # This formula can have rounding errors for antipodal cases.
    D = (np.sin((Dec2-Dec1)/2))**2 + np.cos(Dec1)*np.cos(Dec2)*(np.sin((RA2-RA1)/2))**2
    D = np.sqrt(D)
    D = 2.0 * np.arcsin(D)
    return D

# Being bad and copying code over from selfcal
def healbin(ra, dec, values, nside=128, reduceFunc=np.mean, dtype='float'):
    """Take arrays of ra, dec, and value and bin into healpixels  """

    lat = np.pi/2. - dec
    hpids = hp.ang2pix(nside, lat, ra)

    order = np.argsort(hpids)
    hpids = hpids[order]
    values = values[order]
    pixids = np.arange(hp.nside2npix(nside))

    left = np.searchsorted(hpids, pixids)
    right = np.searchsorted(hpids, pixids, side='right')

    mapVals = np.zeros(pixids.size, dtype=dtype)+hp.UNSEEN

    for idx in pixids:
        if left[idx] != right[idx]:
            mapVals[idx] = reduceFunc(values[left[idx]:right[idx]] )

    return mapVals

def spectral_histogram(specgram,bins=None,lowBin=None,highBin=None,nbins=None):
    """@calculate spectral histogram from spectrogram

    @param specgram
        spectrogram structure
    @param bins
        spectral bins
    @param lowBin
        low bin
    @param highBin
        high bin
    @param nbins
        number of spectral bins

    """

    # Define bins for the spectral variation histogram
    if lowBin == None:
        lowBin = np.log10(np.min(specgram)/2)
    if highBin == None:
        highBin = np.log10(np.max(specgram)*2)
    if nbins == None:
        nbins = 500
    if bins == None:
        bins = np.logspace(lowBin,highBin,num=nbins)

    # Ensure we work with numpy array data
    data = np.array(specgram)

    spectral_variation_norm = []
    rows, columns = data.shape

    # Loop over frequencies
    for i in xrange(columns):
        # calculate histogram for this frequency bin
        this_spectral_variation, bin_edges = np.histogram(data[:,i],bins)
        this_spectral_variation = np.array(this_spectral_variation)
        # Calculate weights for bins (to normalize)
        weight = (100/float(sum(this_spectral_variation))) + np.zeros(this_spectral_variation.shape)
        # stack output array
        if spectral_variation_norm == []:
            spectral_variation_norm = this_spectral_variation * weight
        else:
            spectral_variation_norm = np.vstack([spectral_variation_norm,this_spectral_variation * weight])
    spectral_variation_norm = np.transpose(spectral_variation_norm)

    return bins,spectral_variation_norm

def spectral_percentiles(specvar,bins,percentile):
    """@calculate spectral percentiles from spectral variation histogram

    @param specvar
        spectral variation histogram
    @param bins
        spectral bins
    @param percentile
        percentile of the bins to compute
    """

    # Ensure we work with numpy array data
    data = np.array(specvar)

    percentiles = []
    rows, columns = specvar.shape

    # Loop over frequencies
    for i in xrange(columns):
        # Calculate cumulative sum for array
        cumsumvals = np.cumsum(data[:,i])

        # Find value nearest requested percentile
        abs_cumsumvals_minus_percentile = abs(cumsumvals - percentile)
        minindex = abs_cumsumvals_minus_percentile.argmin()
        val = bins[minindex]

        percentiles.append(val)

    percentiles = np.array(percentiles)
    return percentiles


