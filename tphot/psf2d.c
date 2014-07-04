/* psf2d.c -- fit a 2D Waussian profile */
/* 130108 extend pixmedave to trailed apertures */
/* 121219 fixed a bug in traileval */
/* 121013 added trailpsf().  FIXME: trailpsf flux is junk. */
/* 120929 v1.0 John Tonry adapted from psf2dim.f */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "psf2d.h"

typedef struct {
   int nx;		/* Number of columns in DATA */
   int ny;		/* Number of rows in DATA */
   int NX;		/* Column dimension of DATA */
   int init;		/* Initialize params?
                           0, start with whats in PAR already
                           1, initialize all but x0,y0 = PAR(1,2)
                           2, also find peak to start x0,y0 = PAR(1,2) */
   double x0;		/* X offset of subarray in entire array */
   double y0;		/* Y offset of subarray in entire array */
   double eadu;		/* E/ADU for calculating sigma */
   double sat;		/* saturation value */
   double extra;	/* Added value so that variance = E/ADU*(DATA+EXTRA) */
   double ignore;	/* Pixel value to ignore in fitting */
} PSF2D_AUX;

#define PSFCORE 60	/* Ring count where sum changes to median */
// #define ROBUST_MEDIAN_PHOTOMETRY	/* Photometry from median*area */

#define RESIDCUT 2.0
#define MAXITER 20
#define TRAILITER 20
#define MAXRAD 64
#define MAXPT 512
#define GAUSSRAD 2
#define FWHMGAUSS (2*sqrt(2*log(2.0)))

// #define DEBUG		/* Debug output? */
// #define PRINTRAIL	/* Dump table of data for each trail? */

/* fitlm() fit function prototype */
typedef int LMFUNC(int kpix, void *dataptr, void *usr1, void *usr2, 
		   double *dataval, double *wgt, 
		   int npar, double *par, double *model, double *deriv);

/* Compute the value of a Waussian fit at (x,y)  [leaving out the sky] */
double wauss(double x, double y, PSF2D_PARAM *wpar, double *z2);

/* Compute the value of a trailed Waussian at (x,y)  (leaving out the sky) */
double wtrail(double x, double y, PSF2D_PARAM *wpar, double *z2);

/* Find from a rough position (IX,IY), the highest point PEAK at (MX,MY), the
 * rough FWHM, and a rough BACKground of the image in the wings of the star */
int crudemax(int ix, int iy, int nx, int ny, float *a, 
	     int *mx, int *my, double *peak, double *fwhm, double *back);

/* From a rough position (IX,IY), estimate the center of a trail
 * (MX,MY), its PEAK value, a rough LENG, and a rough BACKground of
 * the image in the wings of the trail */
int trailctr(int RMAX, int ix, int iy, int nx, int ny, float *a, 
	     int *mx, int *my, double *peak, double *leng, double *back);

/* Collect average and median of all annuli out to a radius rmax */
int pixmedave(int mx, int my, double phi, double len,
	      int nx, int ny, float *a, int rmax, 
	      int *np, int *ntot, double *ave, double *med, 
	      double *rms, double *buf);

/* Get a decent estimate of SKY and RMS from the PROFILE and PRMS, 
 * and then get a decent estimate of the FWHM */
int estsky(int np, double *profile, double *prms, 
	   double *fwhm, double *sky, double *err, double *rms);

/* Fit the profile as AMPL * SMOOTH + SKY, and return the error as ERR */
int fitsky(int n, double *profile, double *smooth, 
	   double *ampl, double *amperr, double *sky, double *err);

/* Collect average and sum of the flux out to a radius R */
int pixsum(int mx, int my, int nx, int ny, float *a,
	   int r, int *np, double *ave, double *total, double *peak);

/* Average and sum of the flux in a lozenge of radius rmaj x rap at phi */
int pixlozenge(int mx, int my, int nx, int ny, float *a,
	       int rmaj, int rap, double phi, 
	       int *np, double *ave, double *total, double *peak);

/* Calculate sigma's and position angles from sx2, sxy, sy2 */
int domajmin(PSF2D_PARAM *par);

/* Calculate sx2, sxy, sy2 from sigma's and position angles */
int dosxy(PSF2D_PARAM *par);

/* Evaluate waussian fit and derivatives at a point */
LMFUNC weval;

/* Evaluate a trailed Waussian fit and derivatives at a point */
LMFUNC traileval;

/* Fit a source with a Waussian = Wingy Gaussian */
int waussfit(float *data, int npar, PSF2D_PARAM *par, 
	     PSF2D_AUX *aux, int *niter, double *cov, double *chisq);

/* Fit a Waussian using only PEAK and SKY; he's linear, Jim. */
int wausstwo(float *data, PSF2D_PARAM *par, PSF2D_AUX *aux, 
	     double *cov, double *chisq);

/* Fit a source with a trailed Waussian */
int trailfit(float *data, int npar, PSF2D_PARAM *par, 
	     PSF2D_AUX *aux, int *niter, double *cov, double *chisq);

/* Quicksort routine for doubles: 1980 (John Tonry) */
int qsort8(int n, double *x);

/* fitlm control parameters */
#define FITLM_VERBOSE 0x01	/* Verbose? */
#define FITLM_DERIV   0x02	/* Analytic first derivatives? */
#define FITLM_LINEAR  0x04	/* Linear function? */
#define FITLM_WGT     0x08	/* Data point weights provided? */

/* Fit a model to data, minimizing chi^2 using Levenberg-Marquardt */
int fitlm(int ndata, void *dataptr, void *usr1, void *usr2, LMFUNC func,
	  int npar, double *par, int *usepar, 
	  double *chi, double *cov, int *ctrl);

/* Give user access to QUITFRAC LAMINIT, etc.  Not for the novice! */
int fitlm_param(int maxiter, int laminit, double quit, double parvary);


/*
 * Data in array a, dimensions nx ny
 *
 *   eadu = 1.0;	// Electrons per ADU to estimate noise
 *   sat = 100000;	// Saturation level
 *   aprad = 20;	// Aperture radius for photometry
 *   skyrad = 60;	// Sky radius for photometry
 *   nwpar = 7;		// How many Waussian params to fit?
 *                         if nwpar = 4, wpar->major, minor, and phi
 *                         must be populated!
 *
 * Candidate star position in wpar->x0 wpar->y0
 *
 * If wpar->beta4 < 0 use Gaussian instead of Waussian (requires code change!)
 *
 * Return error code 0: no error, 
 *                  -1: off image, 
 *                   1: fit didn't converge,
 *                N*10: N zero'ed pixels within FWHM
 *
 */
int wpsf(int nx, int ny, float *a, double eadu, double sat, double badata,
          int aprad, int skyrad, int nwpar, 
          PSF2D_PARAM *wpar, PSF2D_FLUX *flux)
{
   double pmed[MAXRAD+1], pave[MAXRAD+1], prms[MAXRAD+1];
   int npix[MAXRAD+1], ntot[MAXRAD+1];
   double smooth[MAXRAD+1], profile[MAXRAD+1];
   double buf[768*MAXRAD]; // 768 = MAXPT (from pixmedave) * sqrt(2) * (1+fudge)

   PSF2D_AUX waux;

   int ix, iy, err, mx, my, mxr;
   double peak, crudefw, crudeback;

   int i, j, m1, n, iwxs, iwys, ntotal, nfit;
   int niter, ierr, i0, i1, j0, j1, npt, nresid;

   double fwhm, esky, eskyerr, ave, rms, ampl, amperr;
   double sky, dsky, apflux, dflux, WRAD, WCHI, extrasky;
   double chisq, x, y, diff, z2, biggie, absave, resid, fw;
   double cov[9*9] = {0.0};

   profile[0] = 0;

   ix = wpar->x0;
   iy = wpar->y0;

/* Abort this star if it's off the image */
   if(ix < 0 || ix > nx-1 || iy < 0 || iy > ny-1) return(-1);

/****************************************************************/
/* Find the highest pixel and get some crude information */
   err = crudemax(ix,iy, nx,ny,a, &mx, &my, &peak, &crudefw, &crudeback);

#ifdef DEBUG
   printf("crudemax: err= %d mx= %d my= %d peak= %.1f crudefw= %.2f crudeback= %.1f\n",
	  err, mx, my, peak, crudefw, crudeback);
#endif

/* Get medians and averages at all radii... */
//      mxr = max(16, min(NINT(15*crudefw),MAXRAD))
   mxr = skyrad;
   if(mxr <= 0 || mxr > MAXRAD) mxr = MAXRAD;
   err = pixmedave(mx,my, 0.0,0.0, nx,ny, a, mxr, npix,ntot, pave,pmed,prms, buf);

/* Get an estimated values for SKY and FWHM */
   err = estsky(mxr+1, pmed,prms, &fwhm, &esky, &eskyerr, &rms);

#ifdef DEBUG
   printf("estsky: err= %d fwhm= %.2f esky= %.1f eskyerr= %.1f rms= %.2f\n", 
	  err, fwhm, esky, eskyerr, rms);
#endif

/* Get a fitted value for SKY and FWHM */
   for(i=0; i<=mxr; i++) smooth[i] = pow((3.*fwhm)/MAX(1,i), 3.0);
   m1 = (mxr+1) / 2;
   nfit = mxr+1 - m1;

   err = fitsky(nfit, pmed+m1, smooth+m1, &ampl, &amperr, &sky, &dsky);

#ifdef DEBUG
   printf("fitsky: err= %d amp= %.1f amperr= %.1f sky= %.1f dsky=%.1f\n", 
	  err, ampl, amperr, sky, dsky);
#endif

#ifdef ROBUST_MEDIAN_PHOTOMETRY
/* Use sum and median to suppress adjacent stars */
   for(i=0, ntotal=0, apflux=0.0; i<=mxr && i<aprad; i++) {
      ntotal += ntot[i];
      apflux += ntot[i] * ((ntot[i] < PSFCORE) ? pave[i] : pmed[i]);
   }
#else
/* Add up the total flux, don't worry about adjacent stars */
   err = pixsum(mx,my,nx,ny,a,aprad, &ntotal,&ave,&apflux,&peak);
#endif

   apflux = apflux - ntotal*sky;
   dflux = apflux/eadu + pow(ntotal*dsky, 2.0) + rms*rms*ntotal;
   dflux = (dflux >= 0) ? sqrt(dflux) : -sqrt(-dflux);
   flux->sky = sky;
   flux->dsky = dsky;
   flux->flux = apflux;
   flux->dflux = dflux;
   flux->skyrms = rms;

#ifdef DEBUG
   printf("pixsum: ntot= %d sky= %.1f dsky= %.1f apflux= %.1f dflux= %.1f RMS= %.2f ave= %.1f peak= %.1f\n",
	  ntotal, sky, dsky, apflux, dflux, rms, ave, peak);
#endif

/************************************************************************/
/* Now do a Waussian fit to the data using all we know for initial values */
/* How many FWHM do we fit? */
   WRAD = 2.5;
/* How many FWHM do we use for evaluating chi/N? */
   WCHI = 2.0;
/* Guess at extrasky to meet noise seen by estsky */
   if(rms != 0 && esky != 0) {
      extrasky = eadu*rms*rms - esky;
   } else {
/* Disable weighting if something's really wrong with the "sky" level and rms */
      extrasky = -1e10;
   }
// WRITE(6,*) EADU, RMS, ESKY, EXTRASKY

/* What do we want for a center guess? */
   if(nwpar >= 7) {
/* Use the positive peak... */
      n = MAX(3, MIN(31, MAX(GAUSSRAD, NINT(WRAD*fwhm))));
      iwxs = MAX(0,mx-n);
      iwys = MAX(0,my-n);
   } else {
      err = dosxy(wpar);
      n = MAX(3, MIN(31, MAX(GAUSSRAD, 
			     NINT(WRAD*sqrt(ABS(wpar->major*wpar->minor))))));
      iwxs = MAX(0, NINT(wpar->x0)-n);
      iwys = MAX(0, NINT(wpar->y0)-n);
   }

/* WPAR: X0, Y0, P, SKY, SX2, SXY, SY2, B4, B6, FWX, FWY, PHI */
   if(nwpar >= 7) {
      wpar->x0 = mx;
      wpar->y0 = my;
   }
/* Waussian fit (WPAR(8) = BETA4 = 1, WPAR(9) = BETA6 = 0.5) */
/* NOTE: wpar->beta4 < 0 to request Gaussian fit instead */
   wpar->beta4 = 1.0;
   wpar->beta6 = 0.5;
/* WAUX: NXPATCH, NYPATCH, NX, XOFF, YOFF, EADU, EXTRASKY, IGNORE_VALUE, INIT */
   waux.nx = MIN(n,nx-1-mx) + n + 1;
   waux.ny = MIN(n,ny-1-my) + n + 1;
   waux.NX = nx;
   waux.x0 = iwxs;
   waux.y0 = iwys;
   waux.eadu = eadu;
   waux.sat = sat;
   waux.extra = extrasky;
/* Bad data value, by convention 0.0 for JT */
   waux.ignore = badata;
   waux.init = 1;
/* Don't dare leap to all 9 params at once */
   if(nwpar > 7) {
      niter = MAXITER;
      err = waussfit(a+iwxs+iwys*nx, 7, wpar, &waux, &niter, cov, &chisq);
   }
   if(nwpar > 2) {
      niter = MAXITER;
      err = waussfit(a+iwxs+iwys*nx, nwpar, wpar, &waux, &niter, cov, &chisq);
   } else {
      err = wausstwo(a+iwxs+iwys*nx, wpar, &waux, cov, &chisq);
   }
   wpar->chin = chisq;

/* Fill in the uncertainties */
   wpar->dx0 = wpar->dy0 = wpar->dpeak = wpar->dsky = wpar->dsx2 =
      wpar->dsxy = wpar->dsy2 = wpar->dbeta4 = wpar->dbeta6 = 0.0;
   if(chisq > 0) {
      if(nwpar <= 2) {
	 wpar->dpeak  = cov[0+0*nwpar];
	 wpar->dsky   = cov[1+1*nwpar];
      } else {
	 if(nwpar >= 4) {
	    wpar->dx0    = cov[0+0*nwpar];
	    wpar->dy0    = cov[1+1*nwpar];
	    wpar->dpeak  = cov[2+2*nwpar];
	    wpar->dsky   = cov[3+3*nwpar];
	 }
	 if(nwpar >= 7) {
	    wpar->dsx2   = cov[4+4*nwpar];
	    wpar->dsxy   = cov[5+5*nwpar];
	    wpar->dsy2   = cov[6+6*nwpar];
	 }
	 if(nwpar >= 9) {
	    wpar->dbeta4 = cov[7+7*nwpar];
	    wpar->dbeta6 = cov[8+8*nwpar];
	 }
      }
   }

   if(niter == MAXITER || chisq == 0) return(niter);

   if(chisq < 0) return(NINT(chisq));

   ierr = 0;
/* Improved estimate of chi^2/N over WCHI(=2) * FWHM, using observed RMS */
   if(rms > 0 && wpar->major > 0 && wpar->minor > 0) {
      fw = sqrt(wpar->major*wpar->minor);
      i0 = MAX(0,    NINT(wpar->x0 - WCHI*fw));
      i1 = MIN(nx-1, NINT(wpar->x0 + WCHI*fw));
      j0 = MAX(0,    NINT(wpar->y0 - WCHI*fw));
      j1 = MIN(ny-1, NINT(wpar->y0 + WCHI*fw));
      chisq = 0;
      npt = 0;
      resid = 0;
      nresid = 0;
      biggie = 0;
      absave = 0;

      for(j=j0; j<=j1; j++) {
	 y = j + 0.5;
	 for(i=i0; i<=i1; i++) {
/* Count bad pixels within r<FWHM, then skip */
	    if(a[i+nx*j] == waux.ignore) {
	       if(pow(i-wpar->x0, 2.0)+pow(j-wpar->y0, 2.0) <= fw*fw) 
		  ierr = ierr + 10;
	       continue;
	    }
	    x = i + 0.5;
	    diff = a[i+nx*j] - (wauss(x, y, wpar, &z2) + wpar->sky);
	    chisq += pow(diff/rms, 2.0);
	    npt++;
/* Get residual statistics within (-0.5*r^2/sig^2) < RESIDCUT(=2) */
	    if(z2 < RESIDCUT) {
	       resid += diff*diff;
	       nresid++;
	       if(ABS(diff) > ABS(biggie)) biggie = diff;
	       absave += ABS(diff);
	    }
	 }
      }
/* Replace wpar->chin with this better chi/N */
      wpar->chin = chisq / MAX(1, npt-nwpar);
      wpar->rms = sqrt(resid / MAX(1, nresid));
      wpar->absresid = absave / MAX(1, nresid);
      wpar->maxresid = biggie;
   }
   return(ierr);
}

/* Compute the value of a Waussian fit at (x,y)  [leaving out the sky] */
double wauss(double x, double y, PSF2D_PARAM *wpar, double *z2)
{
   *z2 = wpar->sx2*(x-wpar->x0)*(x-wpar->x0) + 
         wpar->sxy*(x-wpar->x0)*(y-wpar->y0) +
         wpar->sy2*(y-wpar->y0)*(y-wpar->y0);
   if(wpar->beta4 < 0) {		// Gaussian, not Waussian
      return(wpar->peak * exp(-*z2));
   } else {				// Waussian, not Gaussian
      return(wpar->peak/(1 + (*z2)*(1 + (*z2)*(wpar->beta4/2 + 
					       (*z2)*wpar->beta6/6))));
   }
}

/* Find from a rough position (IX,IY), the highest point PEAK at (MX,MY), the
 * rough FWHM, and a rough BACKground of the image in the wings of the star */
int crudemax(int ix, int iy, int nx, int ny, float *a, 
	 int *mx, int *my, double *peak, double *fwhm, double *back)
{
   double pi=4*atan(1.0);
   double peak2;
   int i, j, k, l, n, im=0, jm=0, idx, idy, maxstep, nstep;

/* Find the local maximum */
   maxstep = 30;
   *mx = ix;
   *my = iy;
   *peak = a[ix+iy*nx];
   for(nstep=0; nstep<maxstep; nstep++) {
      peak2 = -1e10;
      for(j=*my-1; j<=*my+1; j++) {
	 if(j > ny-1 || j < 0) continue;
	 for(i=*mx-1; i<=*mx+1; i++) {
	    if(i > nx-1 || i < 0) continue;
	    if( a[i+j*nx] > peak2) {
	       peak2 = a[i+j*nx];
	       im = i;
	       jm = j;
	    }
	 }
      }
      if(peak2 <= *peak) break;
      *mx = im;
      *my = jm;
      *peak = peak2;
   }

/* Run down the profile in the +/-x and y directions to where it turns up */
   *back = 0;
   n = 0;
   for(k=0; k<4; k++) {
      idx = NINT(cos(k*pi/2));
      idy = NINT(sin(k*pi/2));
      for(l=1; l<maxstep; l++) {
	 i = *mx + idx*l;
	 j = *my + idy*l;
	 if(i >= nx-1 || i <= 0 || j >= ny-1 || j <= 0) break;
	 if(a[i+j*nx] > a[i-idx+(j-idy)*nx]) {
	    n++;
	    *back += a[i+j*nx];
	    break;
	 }
      }
   }
   *back /= MAX(n,1);

/* Find the FWHM using this value for BACK */
   *fwhm = 0;
   n = 0;
   for(k=0; k<4; k++) {
      idx = NINT(cos(k*pi/2));
      idy = NINT(sin(k*pi/2));
      for(l=1; l<maxstep; l++) {
	 i = *mx + idx*l;
	 j = *my + idy*l;
	 if(i >= nx || i < 0 || j >= ny || j < 0) break;
	 if(a[i+j*nx] < *peak-(*peak-*back)/2) {
	    n++;
	    *fwhm += l;
	    break;
	 }
      }
   }
   *fwhm = *fwhm * 2.0 / MAX(1,n) - 1;
   return(0);
}

/* Collect average and median of all annuli out to a radius rmax */
/*  center point mx, my;  streak half-length len, angle phi [rad],
 *  data array a[ny][nx]
 *  collect arrays out to a radius rmax:
 *     np[] = number of pixels in annulus falling within the array
 *   ntot[] = total number of pixels in annulus
 *    ave[] = average value in annulus
 *    med[] = median within annulus, trimmed of bumps
 *    rms[] = median-derived RMS within annulus, trimmed of bumps
 *    buf[] = buffer of annuli, size MAXPT*768 (= sqrt(2) * (1+fudge))
 */

int pixmedave(int mx, int my, double phi, double len,
	      int nx, int ny, float *a, int rmax, 
	      int *np, int *ntot, double *ave, double *med, 
	      double *rms, double *buf)
{
   int i, j, k, jp, ix1, ix2, err, nsig, n, n1, n2, nuke;
   double r, pixel, amed, arms, trigger, reset, c, s;
   for(i=0; i<=rmax; i++) {
      np[i] = ntot[i] = 0;
      ave[i] = 0;
   }

/* Accumulate pixels from the bottom to top on the right and
 * then top to bottom on the left to get the pixels in azimuthal order */
   c = cos(phi);
   s = sin(phi);

   for(jp=0; jp<=2*(2*rmax+1); jp++) {
      if(jp < 2*rmax+1) {
	 j = jp - rmax;
	 ix1 = 0;
	 ix2 = rmax;
      } else {
	 j = 3*rmax + 1 - jp;
	 ix1 = -rmax;
	 ix2 = -1;
      }
      for(i=ix1; i<=ix2; i++) {
/* Lozenge shape around trail */
	 if(i*c+j*s > len) {
	    r = sqrt((i-c*len)*(i-c*len)+(j-s*len)*(j-s*len));
	 } else if(i*c+j*s < -len) {
	    r = sqrt((i+c*len)*(i+c*len)+(j+s*len)*(j+s*len));
	 } else {
	    r = ABS(-i*s+j*c);
	 }
//	 r = sqrt( (double)(i*i+j*j));
	 k = NINT(r);
	 if(k > rmax) continue;
	 ntot[k] += 1;
	 if(j+my >= ny-1 || j+my < 0) continue;
	 if(i+mx >= nx-1 || i+mx < 0) continue;
	 pixel = a[i+mx + (j+my)*nx];
	 buf[np[k] + k*MAXPT] = pixel;
	 np[k] += 1;
	 if(np[k] > MAXPT-1) {
//               write(6,*) 'pixmedave: np exceeded MAXPT!'
	    np[k] = MAXPT-1;
	 }
	 ave[k] += pixel;
//	 if(k==0) printf("k=0 %d %.1f %.1f\n", np[k], pixel, buf[0]);
      }
   }

/* Now compute medians.  If there are sufficient points, correct the median
 * for skewness by assessing the rms, and counting all points around
 * points of greater than 3-sigma which themselves exceed 1-sigma, and
 * throwing out those points from the sorted data.
 */
   for(k=0; k<=rmax; k++) {
      ave[k] /= MAX(np[k], 1);
      for(i=0; i<np[k]; i++) buf[i] = buf[i+k*MAXPT];
      err = qsort8(np[k], buf);
      amed = 0.5*(buf[(np[k]-1)/2] + buf[np[k]/2]);
      nsig = MAX(0, NINT(0.1587*np[k]));
      arms = amed - buf[nsig];
//      printf("%d %d %.1f %.1f %.1f\n", k, np[k], buf[0], amed, arms);

/* Don't do any funny stuff with radii less than 10 pixels... (N ~ 60) */
      if(np[k] <= PSFCORE) {
	 med[k] = amed;
	 rms[k] = arms;
      } else {
	 trigger = amed + 3*arms;
//            reset = amed + arms
	 reset = amed;

/* Find a spot which is lower than the median */
	 n1 = 1;
	 while(buf[n1+k*MAXPT] > reset) n1++;

/* Now advance around the circle, counting all pixels higher than RESET
 * which are adjacent to a pixel higher than TRIGGER */
	 nuke = 0;
/* But fed to a mod we want to start at n1+1 - 1 */
	 n2 = n1;
	 while(n2-n1 < np[k]) {
	    if(buf[(n2%np[k])+k*MAXPT] > trigger) {

/* Back up to find out where the pixels higher than RESET began */
	       i = n2 - 1;
	       while(buf[(i%np[k])+k*MAXPT] > reset) {
		  nuke++;
		  i--;
	       }

/* Advance to find out where the pixels higher than RESET end */
	       while(buf[(n2%np[k])+k*MAXPT] > reset) {
		  nuke++;
		  n2++;
	       }
	    } else {
	       n2++;
	    }
	 }

/* Now recompute the median and rms with NUKE pixels removed off the top */
	 n = np[k] - nuke;
	 med[k] = 0.5*(buf[(n-1)/2] + buf[n/2]);
	 rms[k] = med[k] - buf[NINT(0.1587*n)];
      }
   }
   return(0);
}

/* Collect average and sum of the flux out to a radius IR */
int pixsum(int mx, int my, int nx, int ny, float *a,
	   int ir, int *np, double *ave, double *total, double *peak)
{
   int i, j, ix, iy, ir2;

   *total = *peak = 0.0;
   *np = 0;

   for(j=-ir; j<=ir; j++) {
      iy = j + my;
      if(iy < 0 || iy >= ny) continue;
      for(i=-ir; i<=ir; i++) {
	 ix = i + mx;
	 if(ix < 0 || ix >= nx) continue;
	 ir2 = i*i + j*j;
	 if(ir2 > ir*ir) continue;
	 *np += 1;
	 *total += a[ix+iy*nx];
	 *peak = MAX(*peak, a[ix+iy*nx]);
      }
   }
   *ave = *total / MAX(1, *np);
   return(0);
}

/* Collect average and sum of the flux out to a radius IR */
int pixlozenge(int mx, int my, int nx, int ny, float *a,
	       int rmaj, int rap, double phi, 
	       int *np, double *ave, double *total, double *peak)
{
   int i, j, ix, iy;
   double c, s, r;

   *total = *peak = 0.0;
   *np = 0;

   c = cos(phi);
   s = sin(phi);

   for(j=-rmaj; j<=rmaj; j++) {
      iy = j + my;
      if(iy < 0 || iy >= ny) continue;
      for(i=-rmaj; i<=rmaj; i++) {
	 ix = i + mx;
	 if(ix < 0 || ix >= nx) continue;
	 if(i*c+j*s > rmaj) {
	    r = sqrt((i-c*rmaj)*(i-c*rmaj)+(j-s*rmaj)*(j-s*rmaj));
	 } else if(i*c+j*s < -rmaj) {
	    r = sqrt((i+c*rmaj)*(i+c*rmaj)+(j+s*rmaj)*(j+s*rmaj));
	 } else {
	    r = ABS(-i*s+j*c);
	 }
	 if(NINT(r) > rap) continue;
	 *np += 1;
	 *total += a[ix+iy*nx];
	 *peak = MAX(*peak, a[ix+iy*nx]);
      }
   }
   *ave = *total / MAX(1, *np);
   return(0);
}


/* Get a decent estimate of SKY and RMS from the PROFILE and PRMS, 
 * and then get a decent estimate of the FWHM */
int estsky(int np, double *profile, double *prms, 
	   double *fwhm, double *sky, double *err, double *rms)
{
   int i, n, minus, ir1=1, ir2=0, ierr;
   double buf[21], half, budge;

/* First get SKY from a median of the last points in PROFILE */
   n = MIN(np/2, 21);
   for(i=0; i<n; i++) buf[i] = profile[np-1-i];
   ierr = qsort8(n, buf);
   *sky = 0.5*(buf[(n-1)/2] + buf[n/2]);
   *err = (*sky - buf[MAX(0,(n+1)/6)]) / sqrt((double)n);

/* Now get FWHM from the estimated SKY and the central intensity */
   minus = (profile[0] >= *sky) ? +1 : -1;
   half = 0.5*(profile[0] + *sky);
   for(i=0; i<np; i++) {
//      printf("%d %.1f %.1f %d\n", i, profile[i], half, minus);
      if(minus*profile[i] <= minus*half) {
	 ir2 = i;
	 break;
      }
      if(minus*profile[i] > minus*half) ir1 = i;
   }


   budge = 0;
   if(profile[ir1] != profile[ir2]) {
      budge = (ir2-ir1)*(profile[ir1]-half) / (profile[ir1]-profile[ir2]);
//      printf("%d %d %.3f\n", ir1, ir2, budge);
      *fwhm = 2*(ir1 + budge);
   } else {
      *fwhm = 2*ir2;
   }

   if(*fwhm < 0) *fwhm = 0.5;

/* Finally get the RMS as the median of the last several rms points */
   for(i=0; i<n; i++) buf[i] = prms[np-1-i];
   ierr = qsort8(n, buf);
   *rms = 0.5*(buf[(n-1)/2] + buf[n/2]);

   return(0);
}

/* Fit the profile as AMPL * SMOOTH + SKY, and return the error as ERR */
int fitsky(int n, double *profile, double *smooth, 
	   double *ampl, double *amperr, double *sky, double *err)
{
   int j;
   double v[2], am[3], sum, sum2, det, tmp, rms;

   v[0] = v[1] = am[0] = am[1] = am[2] = 0.0;

/* Accumulate sums */
   for(j=0; j<n; j++) {
//      printf("%d %.1f %.3e\n", j, profile[j], smooth[j]);
      v[0] += profile[j];
      v[1] += profile[j] * smooth[j];
      am[0] += 1.0;
      am[1] += smooth[j] * smooth[j];
      am[2] += smooth[j];
   }

   det = am[0]*am[1] - am[2]*am[2];
   if(det == 0.0) {
      *sky = *ampl = *amperr = *err = 0.0;
      return(-1);
   }

   tmp = am[1];
   am[1] = am[0] / det;
   am[0] = tmp / det;
   am[2] = -am[2] / det;

   *sky =  am[0]*v[0] + am[2]*v[1];
   *ampl = am[2]*v[0] + am[1]*v[1];

   if(*ampl < 0) *ampl = 0.0;

   sum = sum2 = 0.0;
   for(j=0; j<n; j++) {
      sum += profile[j] - (*sky + *ampl*smooth[j]);
      sum2 += pow(profile[j] - (*sky + *ampl*smooth[j]), 2.0);
   }

   sum /= n;
   rms = sqrt(ABS(sum2/n-sum*sum));

   if(*ampl == 0.0) *sky = *sky + sum;

   *err = sqrt(am[0]) * rms;
   *amperr = sqrt(am[1]) * rms;
//   cov = am[2] / sqrt(am[0]*am[1]);

   return(0);
}


#define MAXSRT 257

/* From a rough position (IX,IY), estimate the center of a trail
 * (MX,MY), its PEAK value, a rough LENG, and a rough BACKground of
 * the image in the wings of the trail
 * RMAX = half size of search box for sky
 */
int trailctr(int RMAX, int ix, int iy, int nx, int ny, float *a, 
	 int *mx, int *my, double *peak, double *leng, double *back)
{
   double SRTBUF[MAXSRT];
   double thresh;
   int i, j, k, l, n, imin, imax, jmin, jmax, i0, i1;

/* Get the local background */
   l = (2*RMAX+1)*(2*RMAX+1) / MAXSRT;	// Stride through search box
   n = 0;
   for(k=0; k<MAXSRT; k++) {
      i = ((l*k) % (2*RMAX+1)) + ix-RMAX;
      j = ((l*k) / (2*RMAX+1)) + iy-RMAX;
      if(j>ny-1 || j<0 || i>nx-1 || i<0) continue;
      SRTBUF[n++] = a[i+j*nx];
   }
   i = qsort8(n, SRTBUF);
   *back = SRTBUF[(2*n)/5];		// A bit below median

   thresh = 0.5*(a[ix+iy*nx] + *back);	// Threshold for being part of trail

/* Sniff around for the center of the trail and its length */
   imin = imax = ix;
   i0 = i1 = ix;			// Range above thresh in previous row
   for(jmax=iy; jmax<MIN(ny,iy+RMAX); jmax++) {	// Step up from center
// Is there an abutting pixel that is greater than thresh?
      for(i0=MAX(0,i0-1); i0<=MIN(nx-1,i1+1) && a[i0+jmax*nx]<thresh; i0++);
      if(i0>=nx || a[i0+jmax*nx]<thresh) break;	// No pixel above thresh
      for(        ; i0>=0 && a[i0+jmax*nx]>=thresh; i0--); // Extend left
      for(i1=i0+1 ; i1<nx && a[i1+jmax*nx]>=thresh; i1++); // Extend right
      imin = MIN(imin, i0+1);
      imax = MAX(imax, i1-1);
   }
   i0 = i1 = ix;			// Range above thresh in previous row
   for(jmin=iy; jmin>=MAX(0,iy-RMAX); jmin--) {	// Step down from center
// Is there an abutting pixel that is greater than thresh?
      for(i0=MAX(0,i0-1); i0<=MIN(nx-1,i1+1) && a[i0+jmin*nx]<thresh; i0++);
      if(i0>=nx || a[i0+jmin*nx]<thresh) break;	// No pixel above thresh
      for(        ; i0>=0 && a[i0+jmin*nx]>=thresh; i0--); // Extend left
      for(i1=i0+1 ; i1<nx && a[i1+jmin*nx]>=thresh; i1++); // Extend right
      imin = MIN(imin, i0+1);
      imax = MAX(imax, i1-1);
   }

   *mx = (imin+imax) / 2;
   *my = (jmin+jmax) / 2;
   *leng = (imax-imin)*(imax-imin) + (jmax-jmin)*(jmax-jmin);
   *leng = sqrt(*leng);
   *peak = a[ix+iy*nx];

#ifdef DEBUG   
   printf("\ncrud: i,j= %d %d  %d %d thresh= %.1f x,y= %d %d l=%.1f p=%.1f b=%.1f\n", 
	  imin, imax, jmin, jmax, thresh, *mx, *my, *leng, *peak, *back);
#endif

   return(0);
}

/* Calculate sigma's and position angles from sx2, sxy, sy2 */
int domajmin(PSF2D_PARAM *par)
{
   double pi=4*atan(1.0), gausshalf=2.3548, tmp;
   if(par->sx2 == par->sy2) {
      par->phi = pi/4;
      par->major = par->minor = par->sx2 + par->sy2;	/* sx2 = 0.5/sig/sig */
   } else {
      par->phi = 0.5*atan(par->sxy / (par->sx2 - par->sy2));
      par->major = par->sx2 + par->sy2 + (par->sx2-par->sy2) / cos(2*par->phi);
      par->minor = par->sx2 + par->sy2 - (par->sx2-par->sy2) / cos(2*par->phi);
      if(par->major > par->minor) {	/* sic */
	 tmp = par->major;
	 par->major = par->minor;
	 par->minor = tmp;
	 par->phi = par->phi - pi/2;
      }
   }
   if(par->phi < 0) par->phi = par->phi + pi;

/* Take careful square roots */
   if(par->major > 0) {
      par->major = gausshalf / sqrt(par->major);
   } else if(par->major == 0) {
      par->major = 0;
   } else {
      par->major = -gausshalf / sqrt(-par->major);
   }

   if(par->minor > 0) {
      par->minor = gausshalf / sqrt(par->minor);
   } else if(par->minor == 0) {
      par->minor = 0;
   } else {
      par->minor = -gausshalf / sqrt(-par->minor);
   }
   return(0);
}

/* Calculate sx2, sxy, sy2 from sigma's and position angles */
int dosxy(PSF2D_PARAM *par)
{
   double gausshalf=2.3548, majs2, mins2, sum, diff;

   if(par->major <= 0 || par->minor <= 0) return(-1);

   if(par->major == par->minor) {
      par->sxy = 0.0;
      par->sx2 = par->sy2 = 0.5 * pow(gausshalf/par->major, 2.0);
      return(0);
   }

   majs2 = pow(gausshalf / par->major, 2.0);
   mins2 = pow(gausshalf / par->minor, 2.0);

   sum = 0.5 * (majs2 + mins2);
   diff = 0.5 * (majs2 - mins2) * cos(2*par->phi);

   par->sx2 = 0.5 * (sum + diff);
   par->sy2 = 0.5 * (sum - diff);
   par->sxy = diff * tan(2*par->phi);

   return(0);
}

/* Evaluate waussian fit and derivatives at a point */
int weval(int kpix, void *dataptr, void *auxvoid, void *usr2, 
	  double *ydat, double *wgt, 
	  int npar, double *parvoid, double *yfit, double *dyda)
{
   double half=0.5, sixth=1.0/6.0;
   double x, y, z2, fn, dfdz2;
   int i, j, err;
   float *data=(float *)dataptr;
   PSF2D_AUX *aux=(PSF2D_AUX *)auxvoid;
   PSF2D_PARAM *par=(PSF2D_PARAM *)parvoid;

   if(aux->nx == 0) return(-1);

   i = kpix % aux->nx;
   j = kpix / aux->nx;

   x = i + 0.5;
   y = j + 0.5;

   z2 = par->sx2*(x-par->x0)*(x-par->x0) + 
        par->sxy*(x-par->x0)*(y-par->y0) + 
        par->sy2*(y-par->y0)*(y-par->y0);
   fn = dfdz2 = 0.0;
   if(z2 < 85) {
      if(par->beta4 < 0 && npar <= 7) {		// Gaussian, not Waussian
	 fn = exp(-z2);
	 dfdz2 = 1/fn;
      } else {					// Waussian, not Gaussian
	 fn = 1/(1 + z2*(1 + z2*(half*par->beta4 + z2*sixth*par->beta6)));
	 dfdz2 = 1 + z2*(par->beta4 + z2*half*par->beta6);
      }
   }

   *yfit = par->sky + par->peak*fn;

   dyda[0] = par->peak*fn*fn*dfdz2*(2*par->sx2*(x-par->x0) + 
				    par->sxy*(y-par->y0));
   dyda[1] = par->peak*fn*fn*dfdz2*(par->sxy*(x-par->x0) + 
				    2*par->sy2*(y-par->y0));
   dyda[2] = fn;
   dyda[3] = 1.0;
   if(npar > 4) {
      dyda[4] = -par->peak*fn*fn*dfdz2*(x-par->x0)*(x-par->x0);
      dyda[5] = -par->peak*fn*fn*dfdz2*(x-par->x0)*(y-par->y0);
      dyda[6] = -par->peak*fn*fn*dfdz2*(y-par->y0)*(y-par->y0);
   }
   if(npar > 7) {
      dyda[7] = -par->peak*fn*fn*half*z2*z2;
      dyda[8] = -par->peak*fn*fn*sixth*z2*z2*z2;
   }

   err = 0;
   *ydat = data[i+aux->NX*j];

   if(aux->extra > -9e9) {
      *wgt = aux->eadu/(*ydat+aux->extra);
      if(*wgt < 0) *wgt = 1;
   } else {
      *wgt = 1;
   }

   if(*ydat == aux->ignore || *ydat > aux->sat) {
      *ydat = *yfit;
      *wgt = 0;
      err = 1;
   }

   return(err);
}

/* int waussfit(float *data, int npar, PSF2D_PARAM *wpar, 
 *	     PSF2D_AUX *waux, int *niter, double *cov, double *chisq)
 */
// *
// *     Input:	data	Pixel array
// *     NPAR      How many parameters to fit: 4, 7, 9
// *               Note: NPAR=9 wants parameters to be first fitted with NPAR=7!
// *     NITER     Max number of iterations requested
// * Parameter structure wpar (input)
//   double x0;		/* initial guess for x0 of fit (TSK convention!) */
//   double y0;		/* initial guess for y0 of fit */
//   double beta4;	/* r^4 coeff of Waussian */
//   double beta6;	/* r^6 coeff of Waussian */
// * Auxilary structure waux (input)
//   int nx;		/* Number of columns in DATA */
//   int ny;		/* Number of rows in DATA */
//   int NX;		/* Column dimension of DATA */
//   int init;		/* Initialize params?
//                           0, start with whats in PAR already
//                           1, initialize all but x0,y0 = PAR(1,2)
//                           2, also find peak to start x0,y0 = PAR(1,2) */
//   double x0;		/* X offset of subarray in entire array */
//   double y0;		/* Y offset of subarray in entire array */
//   double eadu;	/* E/ADU for calculating sigma */
//   double extra;	/* Added value so that variance = E/ADU*(DATA+EXTRA) */
//   double ignore;	/* Pixel value to ignore in fitting */
// * Parameter structure wpar (output)
//   double x0;		/* initial guess for x0 of fit (TSK convention!) */
//   double y0;		/* initial guess for y0 of fit */
//   double peak;	/* peak of Waussian */
//   double sky;	/* sky of Waussian fit */
//   double sx2;	/* sx2 */
//   double sxy;	/* sxy */
//   double sy2;	/* sy2 */
//   double beta4;	/* r^4 coeff of Waussian (untouched if npar <= 7 */
//   double beta6;	/* r^6 coeff of Waussian (untouched if npar <= 7 */
//   double major;	/* major axis fwhm of Waussian fit */
//   double minor;	/* minor axis fwhm of Waussian fit */
//   double phi;	/* angle of major axis (CCW from x axis) [rad] */
// *      NITER 	Number of iterations carried out
/* Fit a source with a Wingy Gaussian */
int waussfit(float *data, int npar, PSF2D_PARAM *wpar, 
	     PSF2D_AUX *waux, int *niter, double *cov, double *chisq)
{
   int i, j, nx, ny, nxdim, init, imax, jmax;

   double half, width, d0, d1;
   int err, usepar[9], ctrl;

   nx = NINT(waux->nx);
   ny = NINT(waux->ny);
   nxdim = NINT(waux->NX);
   init = NINT(waux->init);

   for(i=0; i<9; i++) usepar[i] = i < npar;

/* Set up initial values for parameters... */
   if(init > 1) {
/* Find the peak for INIT .GE. 2 */
      imax = jmax = 0;
      for(j=0; j<ny; j++) {
	 for(i=0; i<nx; i++) {
	    if(data[i+j*nxdim] > data[imax+jmax*nxdim]) {
	       imax = i;
	       jmax = j;
	    }
	 }
      }
      wpar->x0 = imax - 0.5 + waux->x0;
      wpar->y0 = jmax - 0.5 + waux->y0;
   }

/* Internal to waussfit the position is relative to the subarray... */
   wpar->x0 -= waux->x0;
   wpar->y0 -= waux->y0;

//   printf("%d %d %d %d %.1f %.1f %.1f %.1f\n", 
//   nx, ny, nxdim, init, wpar->x0, wpar->y0, waux->x0, waux->y0);


/* Continue initialization for INIT .GE. 1 */
   if(init > 0) {
      imax = wpar->x0;
      jmax = wpar->y0;

      d0 = MIN(data[0], data[(ny-1)*nxdim]);
      d1 = MIN(data[nx-1], data[nx-1+(ny-1)*nxdim]);
      wpar->sky = MIN(d0, d1);
      wpar->peak = data[imax+jmax*nxdim] - wpar->sky;
      half = 0.5*(wpar->peak-wpar->sky) + wpar->sky;
//      printf("%d %d %.1f %.1f %.1f %.1f %.1f %.1f\n", 
//	     imax, jmax, d0, d1, wpar->sky, wpar->peak, half, data[imax+jmax*nxdim]);
      if(npar > 4) {
	 width = 1;
	 for(i=imax; i<nx; i++) {
	    if(data[i+jmax*nxdim] < half) {
	       width = i - imax;
	       break;
	    }
	 }
	 if(i == nx) width = nx / 2;
	 wpar->sx2 = 1 / pow(width/1.2, 2.0);
	 wpar->sxy = 0.001;
	 wpar->sy2 = wpar->sx2;
      }
      if(npar > 7) {
	 wpar->beta4 = 1;
	 wpar->beta6 = 1;
      }
   }

   ctrl = FITLM_DERIV + FITLM_WGT;
//   ctrl = FITLM_DERIV + FITLM_WGT + FITLM_VERBOSE;

   err = fitlm(nx*ny, (void *)data, (void *)waux, NULL, weval,
	       npar, (double *)wpar, usepar, chisq, cov, &ctrl);

   *niter = ctrl / 256;		/* Iteration count */

#ifdef DEBUG
   printf("waussfit: %2d %2d %2d x0= %6.2f y0= %6.2f peak= %7.1f sky= %7.1f P= %6.4f %6.4f %6.4f chi= %6.2f nx= %3d\n", err, ctrl, *niter,
	     wpar->x0, wpar->y0, wpar->peak, wpar->sky,
	     wpar->sx2, wpar->sxy, wpar->sy2, *chisq, nx);
#endif

/* Patch up and fill in the parameters */
   wpar->x0 += waux->x0;
   wpar->y0 += waux->y0;

/* Calculate sigma's and position angles from sx2, sxy, sy2 */
   err = domajmin(wpar);

   return(0);
}


/* Fit a Waussian using only PEAK and SKY; he's linear, Jim. */
int wausstwo(float *data, PSF2D_PARAM *par, PSF2D_AUX *aux, 
	     double *cov, double *chisq)
{
   double v[2], det;
   double x, y, z2, fit, dat, wgt=1.0;
   int i, j;

/* Internal to waussfit the position is relative to the subarray... */
   par->x0 -= aux->x0;
   par->y0 -= aux->y0;

   v[0] = v[1] = cov[0] = cov[1] = cov[2] = cov[3] = 0.0;
/* Accumulate sums */
   for(j=0; j<aux->ny; j++) {
      for(i=0; i<aux->nx; i++) {
	 x = i + 0.5;
	 y = j + 0.5;
	 dat = data[i+j*aux->NX];

	 if(dat == aux->ignore || dat > aux->sat) continue;
	 if(aux->extra > -9e9) wgt = aux->eadu/(dat+aux->extra);
	 if(aux->extra <= -9e9 || wgt < 0 ) wgt = 1.0;

	 z2 = par->sx2*(x-par->x0)*(x-par->x0) +
	      par->sxy*(x-par->x0)*(y-par->y0) +
	      par->sy2*(y-par->y0)*(y-par->y0);
	 if(par->beta4 < 0) {	// Gaussian
	    fit = exp(-z2);
	 } else {		// Waussian
	    fit = 1 / (1+z2*(1+z2*(par->beta4+z2*par->beta6)));
	 }
	 v[0] += wgt*dat;
	 v[1] += wgt*dat*fit;
	 cov[0] += wgt;
	 cov[1] += wgt*fit;
	 cov[3] += wgt*fit*fit;
	 cov[2] += wgt*dat*dat;		// Accumulate for chi^2 only
//	 printf("%3d %3d %9.1e %6.1f %9.4f %8.1f %8.1f %8.1f %8.1f %8.1f\n", 
//		i, j, wgt, dat, fit, v[0], v[1], cov[0], cov[1], cov[3]);
      }
   }

   det = cov[0]*cov[3] - cov[1]*cov[1];
   if(det == 0.0) {
      par->peak = par->sky = 0.0;
      *chisq = -2.0;
      return(-1);
   }
   par->sky =  ( cov[3]*v[0] - cov[1]*v[1]) / det;
   par->peak = (-cov[1]*v[0] + cov[0]*v[1]) / det;
   *chisq = cov[2] + par->sky*par->sky*cov[0] + par->peak*par->peak*cov[3]
      - 2*par->sky*v[0] - 2*par->peak*v[1] + 2*par->sky*par->peak*cov[1];

   cov[0] = sqrt(cov[0] / det);		// backwards to conform to param order
   cov[3] = sqrt(cov[3] / det);		// peak then sky
   cov[1] = cov[2] = (-cov[1] / det) / (cov[0]*cov[3]);

/* Fix up the position */
   par->x0 += aux->x0;
   par->y0 += aux->y0;

   domajmin(par);	// In case we have only sx2, sxy, sy2...
//   printf("%6.1f %6.1f %6.1f %6.1f %8.3f\n", 
//	  par->sky, par->peak, cov[0], cov[3], *chisq);
   return(0);
}

/* Evaluate a trailed Waussian fit and derivatives at a point */
int traileval(int kpix, void *dataptr, void *auxvoid, void *usr2, 
	      double *ydat, double *wgt, 
	      int npar, double *parvoid, double *yfit, double *dyda)
{
   double half=0.5, sixth=1.0/6.0;
   double x, y, z2, xp, yp, c, s, x0, y0, phi, tau, len, w, fn, dfdz2, dgdz2;
   int i, j, err;
   float *data=(float *)dataptr;
   PSF2D_AUX *aux=(PSF2D_AUX *)auxvoid;
   PSF2D_PARAM *par=(PSF2D_PARAM *)parvoid;

   if(aux->nx == 0) return(-1);

   i = kpix % aux->nx;
   j = kpix / aux->nx;
   x = i + 0.5;
   y = j + 0.5;

   err = 0;
   *ydat = data[i+aux->NX*j];

   if(aux->extra > -9e9) *wgt = aux->eadu/(*ydat+aux->extra);
   if(aux->extra <= -9e9 || *wgt < 0 ) *wgt = 1.0;

   if(*ydat == aux->ignore || *ydat > aux->sat) {
      *wgt = 0.0;
      err = 1;
   }

   x0 = par->x0;	// 0: streak center
   y0 = par->y0;	// 1: streak center
   w = par->sx2;	// 4: w = 1/sig^2
   phi = par->sxy;	// 5: streak rotated CCW by phi from x axis
   tau = par->sy2;	// 6: full length = 2*exp(tau)
   len = exp(tau);	//    half length

/* Offset and rotate by phi */
   c = cos(phi);
   s = sin(phi);
   xp =  (x-x0) * c + (y-y0) * s;
   yp = -(x-x0) * s + (y-y0) * c;


/* Argument of Waussian */
   z2 = 0.5*w*yp*yp;
   if(xp > len)  z2 += 0.5*w*(xp-len)*(xp-len);
   if(xp < -len) z2 += 0.5*w*(xp+len)*(xp+len);

/* Waussian and 1/(df/dz2) */
   fn = dgdz2 = 0.0;
   if(z2 < 85) {
      if(par->beta4 < 0 && npar <= 7) {		// Gaussian, not Waussian
	 fn = exp(-z2);
	 dgdz2 = 1/fn;
      } else {					// Waussian, not Gaussian
	 fn = 1/(1 + z2*(1 + z2*(half*par->beta4 + z2*sixth*par->beta6)));
	 dgdz2 = 1 + z2*(par->beta4 + z2*half*par->beta6);
      }
   }
   dfdz2 = -par->peak * fn*fn * dgdz2;		// deriv wrt z2

/* Function and yp partial derivatives */
   *yfit = par->sky + par->peak*fn;

   if(dyda == NULL) return(err);
   dyda[0] =  dfdz2 * w*yp*s;		// deriv wrt x0
   dyda[1] = -dfdz2 * w*yp*c;		// deriv wrt y0
   dyda[2] = fn;			// deriv wrt peak
   dyda[3] = 1.0;			// deriv wrt sky
   dyda[4] =  dfdz2 * 0.5*yp*yp;	// deriv wrt w=1/sig^2
   dyda[5] = -dfdz2 * w*yp*xp;		// deriv wrt phi
   dyda[6] = 0.0;			// deriv wrt tau=ln(halflen)

/* add in xp partial derivatives */
   if(xp > len) {
      dyda[0] += -dfdz2 * w*(xp-len)*c;
      dyda[1] += -dfdz2 * w*(xp-len)*s;
      dyda[4] +=  dfdz2 * 0.5*(xp-len)*(xp-len);
      dyda[5] +=  dfdz2 * w*(xp-len)*yp;
      dyda[6] += -dfdz2 * w*(xp-len)*len;
   } else if(xp < -len) {
      dyda[0] += -dfdz2 * w*(xp+len)*c;
      dyda[1] += -dfdz2 * w*(xp+len)*s;
      dyda[4] +=  dfdz2 * 0.5*(xp+len)*(xp+len);
      dyda[5] +=  dfdz2 * w*(xp+len)*yp;
      dyda[6] +=  dfdz2 * w*(xp+len)*len;
   }

   return(err);
}

/* Compute the value of a trailed Waussian at (x,y)  (leaving out the sky) */
double wtrail(double x, double y, PSF2D_PARAM *wpar, double *z2)
{
   double xp, yp, c, s, w, len;

/* Offset and rotate by phi */
   c = cos(wpar->sxy);
   s = sin(wpar->sxy);
   xp =  (x-wpar->x0) * c + (y-wpar->y0) * s;
   yp = -(x-wpar->x0) * s + (y-wpar->y0) * c;
   w = wpar->sx2;
   len = exp(wpar->sy2);

/* Argument of Waussian */
   *z2 = 0.5*w*yp*yp;
   if(xp > len)  *z2 += 0.5*w*(xp-len)*(xp-len);
   if(xp < -len) *z2 += 0.5*w*(xp+len)*(xp+len);
   if(wpar->beta4 < 0) {		// Gaussian, not Waussian
      return(wpar->peak * exp(-*z2));
   } else {				// Waussian, not Gaussian
      return(wpar->peak / (1+(*z2)*(1+(*z2)*
				    (wpar->beta4/2+(*z2)*wpar->beta6/6))));
   }
}

/* Fit a trailed PSF */
int trailpsf(int nx, int ny, float *a, double eadu, double sat, double badata,
          int aprad, int skyrad, int nwpar, PSF2D_PARAM *wpar, PSF2D_FLUX *flux)
{
   double pmed[MAXRAD+1], pave[MAXRAD+1], prms[MAXRAD+1];
   int npix[MAXRAD+1], ntot[MAXRAD+1];
   double smooth[MAXRAD+1], profile[MAXRAD+1];
   double buf[768*MAXRAD]; // 768 = MAXPT (from pixmedave) * sqrt(2) * (1+fudge)
   PSF2D_AUX waux;

   int ix, iy, err, mx, my, mxr;
   double peak, crudefw, crudeback;

   int i, j, m1, n, iwxs, iwys, ntotal, nfit;
   int niter, ierr, i0, i1, j0, j1, npt, nresid;

   double fwhm, esky, eskyerr, ave, rms, ampl, amperr;
   double sky, dsky, apflux, dflux, WRAD, WCHI, extrasky;
   double chisq, x, y, diff, z2, biggie, absave, resid, fw;
   double cov[9*9] = {0.0};

   profile[0] = 0;

   ix = NINT(wpar->x0);
   iy = NINT(wpar->y0);

/* Abort this star if it's off the image */
   if(ix < 0 || ix > nx-1 || iy < 0 || iy > ny-1) return(-1);

/*********************************************************/
/* Find the highest pixel and get some crude information */
   err = trailctr(32, ix,iy, nx,ny,a, &mx, &my, &peak, &crudefw, &crudeback);

#ifdef DEBUG
   printf("crudemax: err= %d ixy= %d %d mxy= %d %d peak= %.1f crudefw= %.2f crudeback= %.1f\n",
	  err, ix, iy, mx, my, peak, crudefw, crudeback);
#endif

/* Get medians and averages at all radii... */
//      mxr = max(16, min(NINT(15*crudefw),MAXRAD))
   mxr = skyrad;
   if(mxr <= 0 || mxr > MAXRAD) mxr = MAXRAD;
   err = pixmedave(mx,my, 0.0, 0.0, nx,ny, a, mxr, npix,ntot, pave,pmed,prms, buf);

/* Get an estimated values for SKY and FWHM */
   err = estsky(mxr+1, pmed,prms, &fwhm, &esky, &eskyerr, &rms);

#ifdef DEBUG
   printf("estsky: err= %d fwhm= %.2f esky= %.1f eskyerr= %.1f rms= %.2f\n", 
	  err, fwhm, esky, eskyerr, rms);
#endif

/* Get a fitted value for SKY and FWHM */
   for(i=0; i<=mxr; i++) smooth[i] = pow((3.*fwhm)/MAX(1,i), 3.0);
   m1 = (mxr+1) / 2;
   nfit = mxr+1 - m1;

   err = fitsky(nfit, pmed+m1, smooth+m1, &ampl, &amperr, &sky, &dsky);

#ifdef DEBUG
   printf("fitsky: err= %d amp= %.1f amperr= %.1f sky= %.1f dsky=%.1f\n", 
	  err, ampl, amperr, sky, dsky);
#endif

#ifdef ROBUST_MEDIAN_PHOTOMETRY
/* Use sum and median to suppress adjacent stars */
   for(i=0, ntotal=0, apflux=0.0; i<=mxr && i<aprad; i++) {
      ntotal += ntot[i];
      apflux += ntot[i] * ((ntot[i] < PSFCORE) ? pave[i] : pmed[i]);
   }
#else
/* Add up the total flux, don't worry about adjacent stars */
   err = pixlozenge(mx,my,nx,ny,a, aprad,aprad,0.0, 
		    &ntotal,&ave,&apflux,&peak);
//   err = pixsum(mx,my,nx,ny,a,aprad, &ntotal,&ave,&apflux,&peak);
#endif

   apflux = apflux - ntotal*sky;
   dflux = apflux/eadu + pow(ntotal*dsky, 2.0) + rms*rms*ntotal;
   dflux = (dflux >= 0) ? sqrt(dflux) : -sqrt(-dflux);
   flux->sky = sky;
   flux->dsky = dsky;
   flux->flux = apflux;
   flux->dflux = dflux;
   flux->skyrms = rms;

#ifdef DEBUG
   printf("pixsum: ntot= %d sky= %.1f dsky= %.1f apflux= %.1f dflux= %.1f RMS= %.2f ave= %.1f peak= %.1f\n",
	  ntotal, sky, dsky, apflux, dflux, rms, ave, peak);
#endif

/************************************************************************/
/* Now do a Waussian fit to the data using all we know for initial values */
   if(nwpar < 7) {
      wpar->x0 = mx + 0.5;
      wpar->y0 = my + 0.5;

// FIXME: this really is a great way to wander off into never-never land!
   } else {
/* How many FWHM do we fit? */
      WRAD = 2.5;
/* How many FWHM do we use for evaluating chi/N? */
      WCHI = 2.0;
/* Highest pixel mx,my, peak, crudeback, sky => average position */
      x = y = z2 = 0.0;
      n = 0;
      for(j=my-WRAD*crudefw; j<=my+WRAD*crudefw; j++) {
	 if(j < 0 || j > ny-1) continue;
	 for(i=mx-WRAD*crudefw; i<=mx+WRAD*crudefw; i++) {
	    if(i < 0 || i > nx-1) continue;
	    if((a[i+j*nx]-crudeback) > 0.3*(peak-crudeback)) {
	       x += (a[i+j*nx]-crudeback) * (i+0.5);
	       y += (a[i+j*nx]-crudeback) * (j+0.5);
	       z2 += a[i+j*nx]-crudeback;
	       n++;
	    }
	 }
      }
/* Use average position as a center guess */
      if(n > 0 && z2 > 0) {
	 wpar->x0 = x / z2;
	 wpar->y0 = y / z2;
      } else {
	 wpar->x0 = mx + 0.5;
	 wpar->y0 = my + 0.5;
      }
   }

/* Guess at extrasky to meet noise seen by estsky */
   if(rms != 0 && esky != 0) {
      extrasky = eadu*rms*rms - esky;
   } else {
/* Disable weighting if something's really wrong with the "sky" level and rms */
      extrasky = -1e10;
   }
// WRITE(6,*) EADU, RMS, ESKY, EXTRASKY

/* How big an area do we fit?  Box of 2n+1  */
   n = MAX(3, MIN(51, MAX(GAUSSRAD, NINT(WRAD*fwhm))));
   if(nwpar < 7) n = MAX(n, wpar->major);

   iwxs = MAX(0,mx-n);
   iwys = MAX(0,my-n);

   wpar->beta4 = 1.0;
   wpar->beta6 = 0.5;
/* WAUX: NXPATCH, NYPATCH, NX, XOFF, YOFF, EADU, EXTRASKY, IGNORE_VALUE, INIT */
   waux.nx = MIN(n,nx-1-mx) + n + 1;
   waux.ny = MIN(n,ny-1-my) + n + 1;
   waux.NX = nx;
   waux.x0 = iwxs;
   waux.y0 = iwys;
   waux.eadu = eadu;
   waux.sat = sat;
   waux.extra = extrasky;
/* Bad data value, by convention 0.0 for JT */
   waux.ignore = badata;
   err = trailfit(a+iwxs+iwys*nx, nwpar, wpar, &waux, &niter, cov, &chisq);
   wpar->chin = chisq;

/* Fill in the uncertainties */
   wpar->dx0 = wpar->dy0 = wpar->dpeak = wpar->dsky = wpar->dsx2 =
      wpar->dsxy = wpar->dsy2 = wpar->dbeta4 = wpar->dbeta6 = 0.0;
   if(niter < MAXITER && chisq > 0) {
      wpar->dx0    = cov[0+0*nwpar];
      wpar->dy0    = cov[1+1*nwpar];
      wpar->dpeak  = cov[2+2*nwpar];
      wpar->dsky   = cov[3+3*nwpar];
      wpar->dsx2   = cov[4+4*nwpar];	// Actually 1/sig^2
      wpar->dsxy   = cov[5+5*nwpar];	// Actually phi
      wpar->dsy2   = cov[6+6*nwpar];	// Actually ln(halflen)
   }

//   if(niter == TRAILITER || chisq == 0) return(niter);
   if(chisq == 0) return(niter);

   if(chisq < 0) return(NINT(chisq));

   ierr = 0;
/* Improved estimate of chi^2/N over WCHI(=2) * FWHM, using observed RMS */
   if(rms > 0 && wpar->major > 0 && wpar->minor > 0) {
      fw = sqrt(wpar->major*wpar->minor);
      i0 = MAX(0,    NINT(wpar->x0 - WCHI*fw));
      i1 = MIN(nx-1, NINT(wpar->x0 + WCHI*fw));
      j0 = MAX(0,    NINT(wpar->y0 - WCHI*fw));
      j1 = MIN(ny-1, NINT(wpar->y0 + WCHI*fw));
      chisq = 0;
      npt = 0;
      resid = 0;
      nresid = 0;
      biggie = 0;
      absave = 0;

      for(j=j0; j<=j1; j++) {
	 y = j + 0.5;
	 for(i=i0; i<=i1; i++) {
/* Count bad pixels within r<FWHM, then skip */
	    if(a[i+nx*j] == waux.ignore) {
	       if(pow(i-wpar->x0, 2.0)+pow(j-wpar->y0, 2.0) <= fw*fw) 
		  ierr = ierr + 10;
	       continue;
	    }
	    x = i + 0.5;
	    diff = a[i+nx*j] - (wtrail(x, y, wpar, &z2) + wpar->sky);
	    chisq += pow(diff/rms, 2.0);
	    npt++;
/* Get residual statistics within (-0.5*r^2/sig^2) < RESIDCUT(=2) */
	    if(z2 < RESIDCUT) {
	       resid += diff*diff;
	       nresid++;
	       if(ABS(diff) > ABS(biggie)) biggie = diff;
	       absave += ABS(diff);
	    }
	 }
      }
/* Replace wpar->chin with this better chi/N */
      wpar->chin = chisq / MAX(1, npt-nwpar);
      wpar->rms = sqrt(resid / MAX(1, nresid));
      wpar->absresid = absave / MAX(1, nresid);
      wpar->maxresid = biggie;
   }

/* Now given a trail shape, redo medians and averages at all radii... */
   err = pixmedave(mx,my, wpar->phi, 0.5*wpar->major, 
		   nx,ny, a, mxr, npix,ntot, pave,pmed,prms, buf);

/* Get a new estimated values for SKY and FWHM */
   err = estsky(mxr+1, pmed,prms, &fwhm, &esky, &eskyerr, &rms);

#ifdef DEBUG
   printf("estsky: err= %d fwhm= %.2f esky= %.1f eskyerr= %.1f rms= %.2f\n", 
	  err, fwhm, esky, eskyerr, rms);
#endif

/* Update fitted value for SKY and FWHM */
   for(i=0; i<=mxr; i++) smooth[i] = pow((3.*fwhm)/MAX(1,i), 3.0);
   err = fitsky(nfit, pmed+m1, smooth+m1, &ampl, &amperr, &sky, &dsky);

#ifdef DEBUG
   printf("fitsky: err= %d amp= %.1f amperr= %.1f sky= %.1f dsky=%.1f\n", 
	  err, ampl, amperr, sky, dsky);
#endif

#ifdef ROBUST_MEDIAN_PHOTOMETRY
/* Use sum and median to suppress adjacent stars */
   for(i=0, ntotal=0, apflux=0.0; i<=mxr && i<aprad; i++) {
      ntotal += ntot[i];
      apflux += ntot[i] * ((ntot[i] < PSFCORE) ? pave[i] : pmed[i]);
   }
#else
/* Add up the total flux, don't worry about adjacent stars */
   err = pixlozenge(mx,my,nx,ny,a, NINT(exp(wpar->sy2)),aprad,wpar->phi,
		    &ntotal,&ave,&apflux,&peak);
//   err = pixsum(mx,my,nx,ny,a,aprad, &ntotal,&ave,&apflux,&peak);
#endif

   apflux = apflux - ntotal*sky;
   dflux = apflux/eadu + pow(ntotal*dsky, 2.0) + rms*rms*ntotal;
   dflux = (dflux >= 0) ? sqrt(dflux) : -sqrt(-dflux);
   flux->sky = sky;
   flux->dsky = dsky;
   flux->flux = apflux;
   flux->dflux = dflux;
   flux->skyrms = rms;

#ifdef DEBUG
   printf("pixsum: ntot= %d sky= %.1f dsky= %.1f apflux= %.1f dflux= %.1f RMS= %.2f ave= %.1f peak= %.1f\n",
	  ntotal, sky, dsky, apflux, dflux, rms, ave, peak);
#endif

   return(ierr);
}

/* Fit a source with a trailed Waussian */
int trailfit(float *data, int npar, PSF2D_PARAM *wpar, 
	     PSF2D_AUX *waux, int *niter, double *cov, double *chisq)
{
   int i, j, n, nx, ny, nxdim, imax, jmax;

   double d0, d1, x2, xy, y2;
   int err, usepar[9], ctrl;

   nx = NINT(waux->nx);
   ny = NINT(waux->ny);
   nxdim = NINT(waux->NX);

   for(i=0; i<9; i++) usepar[i] = i < npar;

/* Internal to trailfit the position is relative to the subarray... */
   wpar->x0 -= waux->x0;
   wpar->y0 -= waux->y0;

#ifdef DEBUG
   printf("n= %d %d %d ctr= %.1f %.1f corner= %.1f %.1f\n", 
	  nx, ny, nxdim, wpar->x0, wpar->y0, waux->x0, waux->y0);
#endif
   imax = wpar->x0;
   jmax = wpar->y0;

/* Initialize parameters... */
   d0 = MIN(data[0], data[(ny-1)*nxdim]);
   d1 = MIN(data[nx-1], data[nx-1+(ny-1)*nxdim]);
   wpar->sky = MIN(d0, d1);
   wpar->peak = data[imax+jmax*nxdim] - wpar->sky;
   x2 = xy = y2 = d1 = 0.0;
   n = 0;
   if(npar > 4) {
/* Collect ~central moments */
      for(j=0; j<ny; j++) {
	 for(i=0; i<nx; i++) {
	    d0 = data[i+j*nxdim] - wpar->sky;
	    if(d0 >= 0.3*wpar->peak) {
	       d1 += d0;
	       x2 += d0 * (i+0.5-wpar->x0) * (i+0.5-wpar->x0);
	       xy += d0 * (i+0.5-wpar->x0) * (j+0.5-wpar->y0);
	       y2 += d0 * (j+0.5-wpar->y0) * (j+0.5-wpar->y0);
	       n++;
	    }
	 }
      }
      if(n>0 && d1>0 && x2>0 && y2>0) {
	 x2 /= d1;
	 xy /= d1;
	 y2 /= d1;
	 wpar->sxy = 0.5*atan2(2*xy, x2-y2);	// angle of trail
/* Covariance matrix (x2 xy / xy y2) eigenvalues */
	 d0 = 0.5 * (x2+y2 + sqrt(4*xy*xy+(x2-y2)*(x2-y2)));
	 d1 = 0.5 * (x2+y2 - sqrt(4*xy*xy+(x2-y2)*(x2-y2)));
	 wpar->sx2 = 1/d1;			// 1/sig^2
	 wpar->sy2 = log(sqrt(d0));		// log(trail half length)

      } else {		// Uh oh, fake it
	 wpar->sx2 = wpar->sy2 = 1.0;
	 wpar->sxy = 0.1;
      }
   } else {
/* 4 parameter fit?  Set wpar "shape" parameters from major, minor, phi */
      wpar->sxy = wpar->phi;
      wpar->sx2 = pow(2.3548/wpar->minor, 2.0);
      wpar->sy2 = sqrt(wpar->major*wpar->major+wpar->minor*wpar->minor) - 
	               wpar->minor;
      if(wpar->sy2 > 0) wpar->sy2 = log(0.5*wpar->sy2);
   }

#ifdef DEBUG
   printf("i,j,n= %d %d %d  dat= %.1f sky= %.1f peak= %.1f  x2,xy,y2= %.1f %.1f %.1f  FWHM,trail,phi= %.2f %.2f %.2f %.1f %.1f\n", 
	  imax, jmax, n, data[imax+jmax*nxdim], wpar->sky, wpar->peak,
	  x2, xy, y2, sqrt(1/wpar->sx2), exp(wpar->sy2), wpar->sxy*57.296, waux->x0, waux->y0);
#endif

#ifdef PRINTRAIL
   for(i=-1; i<nx; i++) printf(" %4d", NINT(i+waux->x0));
   printf("\n");
   for(j=ny-1; j>=0; j--) {
      printf(" %4d", NINT(j+waux->y0));
      for(i=0; i<nx; i++) printf(" %4.0f", data[i+j*nxdim]-wpar->sky);
      printf("\n");
   }
#endif

   ctrl = FITLM_DERIV + FITLM_WGT;
/* Use just numerical derivatives */
//   ctrl = FITLM_WGT;
/* Turn on iteration verbosity? */
//   ctrl += FITLM_VERBOSE;

   err = fitlm_param(TRAILITER/*iter*/, -2/*lam*/, 1e-6/*quit*/, 1e-5/*vary*/);

   err = fitlm(nx*ny, (void *)data, (void *)waux, NULL, traileval,
	       npar, (double *)wpar, usepar, chisq, cov, &ctrl);

   *niter = ctrl / 256;		/* Iteration count */

#ifdef DEBUG
   printf("waussfit: %2d %2d x0= %6.2f y0= %6.2f peak= %7.1f sky= %7.1f P= %6.4f %6.4f %6.4f chi= %6.2f nx= %3d\n", err, ctrl,
	     wpar->x0, wpar->y0, wpar->peak, wpar->sky,
	     wpar->sx2, wpar->sxy, wpar->sy2, 
	  *chisq, nx);
#endif

/* Patch up and fill in the parameters */
   wpar->x0 += waux->x0;
   wpar->y0 += waux->y0;

/* Calculate "minor" = PSF FWHM, "major" = total trail length, and angle */
   wpar->minor = 2.3548 / sqrt(ABS(wpar->sx2));
   wpar->major = pow(2*exp(wpar->sy2)+wpar->minor, 2.0) - pow(wpar->minor, 2.0);
   if(wpar->major > 0) wpar->major = sqrt(wpar->major);
   wpar->phi = wpar->sxy;

   return(0);
}
