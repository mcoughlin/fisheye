/* Read an image, find objects, write table */
/* Usage:
 *
 *        tphot image_file [options]
 *
 */
/* v1.0 120929 John Tonry */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define ABS(a) (((a) > 0) ? (a) : -(a))
#define NINT(x) (x<0?(int)((x)-0.5):(int)((x)+0.5))

/* Failure codes */
#define FAIL_XBORDER  0x0001	/* Failed because outside x border */
#define FAIL_YBORDER  0x0002	/* Failed because outside y border */
#define FAIL_XMOVE    0x0004	/* Failed because x moved too far */
#define FAIL_YMOVE    0x0008	/* Failed because y moved too far */
#define FAIL_MAJBIG   0x0010    /* Failed because major > fwmax */
#define FAIL_MINBIG   0x0020    /* Failed because minor > fwmax */
#define FAIL_MAJSMALL 0x0040    /* Failed because major < fwmin */
#define FAIL_MINSMALL 0x0080    /* Failed because minor < fwmin */
#define FAIL_PEAK     0x0100    /* Failed because peak < netmin */
#define FAIL_CHIN     0x0200    /* Failed because chin > chinmax */
#define FAIL_FLUXNEG  0x0400    /* Failed because flux < 0 */
#define FAIL_FERRNEG  0x0800    /* Failed because dflux < 0 */
#define FAIL_SNR      0x1000    /* Failed because flux/dflux < fitsig */
#define FAIL_OKFIT    0x2000    /* Failed because fit was not OK */

typedef struct {
   double x0;		/* x0 of fit (TSK convention!) */
   double y0;		/* y0 of fit */
   double peak;		/* peak of Waussian */
   double sky;		/* sky of Waussian fit */
   double sx2;		/* sx2 */
   double sxy;		/* sxy */
   double sy2;		/* sy2 */
   double beta4;	/* beta4 */
   double beta6;	/* beta6 */
   double major;	/* major axis fwhm of Waussian fit */
   double minor;	/* minor axis fwhm of Waussian fit */
   double phi;		/* angle of major axis (CCW from x axis) [rad] */
   double chin;		/* chi^2/Ndof of fit */
   double rms;		/* rms residual within a radius of 2 sigma */
   double absresid;	/* average abs residual within a radius of 2 sigma */
   double maxresid;	/* max residual within a radius of 2 sigma */
   double dx0;		/* uncertainty in x0 of fit (TSK convention!) */
   double dy0;		/* uncertainty in y0 of fit */
   double dpeak;	/* uncertainty in peak of Waussian */
   double dsky;		/* uncertainty in sky of Waussian fit */
   double dsx2;		/* uncertainty in sx2 */
   double dsxy;		/* uncertainty in sxy */
   double dsy2;		/* uncertainty in sy2 */
   double dbeta4;	/* uncertainty in beta4 */
   double dbeta6;	/* uncertainty in beta6 */
} PSF2D_PARAM;

typedef struct {
   double flux;		/* flux */
   double dflux;	/* flux uncertainty */
   double sky;		/* sky level */
   double dsky;		/* sky uncertainty */
   double skyrms;	/* rms in sky */
} PSF2D_FLUX;

int syntax(char *prog);

/* Print reasons for failure */
int failcode(int failure);

/* Compute median and median RMS for an array */
int imgstat(int x0, int y0, int x1, int y1, int NX, int NY, 
	    float *data, double *med, double *rms, int *nok);

/* Evaluate the 4 point interpolation of sky at image index k */
double skinterp(int k, int nx, int ny, int nxsky, int nysky, double *skymed,
		double *skyrms, double *rms);

/* Quicksort program for doubles: 1980 (John Tonry) */
int qsort8(int n, double *x);

/* Fit a waussian PSF */
int wpsf(int nx, int ny, float *a, double eadu, double sat, double badata,
	 int aprad, int skyrad, int nwpar, 
	 PSF2D_PARAM *wpar, PSF2D_FLUX *flux);

/* Fit a trailed waussian */
int trailpsf(int nx, int ny, float *a, double eadu, double sat, double badata,
	     int aprad, int skyrad, int nwpar, 
	     PSF2D_PARAM *wpar, PSF2D_FLUX *flux);

/* Compute the value of a Waussian fit at (x,y) (leaving out sky) */
double wauss(double x, double y, PSF2D_PARAM *wpar, double *z2);

/* Compute the value of a trailed waussian fit at (x,y) (leaving out sky) */
double wtrail(double x, double y, PSF2D_PARAM *wpar, double *z2);

int TEST=0;

#define MAXSRT (4093)
#define RESIDRAD (2.5)
#define MAXBIN (10)

/* Global control parameters */
static double srchcut, srchmax, srchrad, srchsig, eadu, bias, sat;
static double fwmax, fwmin, netmin, chinmax, maxmove, fitsig, badata;
static double major, minor, phi;
static int nwpar, aprad, skyrad, okfit, x0border, x1border, y0border, y1border;
static int rdntf, trail, subtract, residsky;
static char *fout, *residfile, *objinput;
static int binfactor, nxsky, nysky;

int main(int argc, char *argv[])
{
   int i, j, npix, nx, NX, ny, sx, sy, failure;
   int ii, jj, err, mefits, bitpix, nfind, nfit, nequal, idxeq, rsub;
   char *imname, *head, *flatfile, *headflat;
   char line[10240];
   int nxflat, nyflat;
   float *data, *resid=NULL, *flat;
   double pi=4*atan(1.0), z2;
   int nskyok, neadu;
   double *skymed, *skyrms, eaduin, sky, rms, flatnorm;
   FILE *fp, *infp=NULL;
   PSF2D_PARAM wpar;
   PSF2D_FLUX wflux;

/* Default parameters */
   bitpix = -32;	// Want floating point data
   binfactor = 1;	// Bin the input data?

   x0border = 2;	// Image border to avoid
   x1border = y0border = y1border = x0border;
   srchsig = 5;		// Minimum value = sky+sig*rms to trigger
   srchcut = 100;	// Minimum value above 0 to trigger
   srchmax = 100000;	// Maximum peak to trigger
   srchrad = 5.0;	// Local max required over radius r to trigger
   fwmax = 100.0;	// Maximum FW to accept
   fwmin = 0.5;		// Minimum FW to accept
   netmin = 10;		// Minimum peak-sky to keep
   maxmove = 3;		// Fit center must coincide with peak this closely
   okfit = 1;		// Does the fit have to be OK to keep?
   fitsig = -1.0;	// Minimum SNR in the flux determination
   chinmax = 1000;	// Maximum chi/N to keep
   bias = 0;		// bias = Image_background - true_sky
   sat = 100000;	// Saturated pixel level
   nxsky = 0;		// Subdivide image in x for sky into N chunks
   nysky = 0;		// Subdivide image in y for sky into N chunks

   fout = "-";		// Output file name
   rdntf = 0;		// Read NITF instead of FITS?
   trail = 0;		// Fit a trailed PSF?

   eaduin = -1.0;	// Electrons per ADU to estimate noise [-1 => auto]
   aprad = 15;		// Aperture radius for photometry
   skyrad = 40;		// Sky radius for photometry

   nwpar = 7;		// How many Waussian params to fit?
   major = -1.0;	// External specification of major axis?
   minor = -1.0;	// External specification of minor axis?
   phi = -100.0;	// External specification of angle of PSF?

   badata = 0.0;	// Bad data value

   residfile = NULL;	// Residual image
   subtract = 0;	// Subtract fit from image as we go?
   residsky = 0;	// Subtract sky estimate from residual image?
   objinput = NULL;	// File with x,y values to fit

   flatfile = NULL;	// Flatfield file
   flatnorm = 0.0;	// Flatfield normalization factor

/* Bomb out if no arguments, now that defaults are set */
   if(argc < 2) {
      syntax(argv[0]);
      exit(0);
   }

   imname = argv[1];	// First argument is the image file

/* Parse the arguments */
   for(i=2; i<argc; i++) {

      if(strcmp(argv[i], "-out") == 0) {
	 fout = argv[++i];

      } else if(strcmp(argv[i], "-obj") == 0) {
	 objinput = argv[++i];

      } else if(strcmp(argv[i], "-resid") == 0) {
	 residfile = argv[++i];

      } else if(strcmp(argv[i], "-subtract") == 0) {
	 subtract = 1;

      } else if(strcmp(argv[i], "-residsky") == 0) {
	 residsky = 1;

      } else if(strcmp(argv[i], "-nitf") == 0) {
	 rdntf = 1;

      } else if(strcmp(argv[i], "-flat") == 0) {
	 flatfile = argv[++i];

      } else if(strcmp(argv[i], "-flatnorm") == 0) {
	 sscanf(argv[++i], "%lf", &flatnorm);

      } else if(strcmp(argv[i], "-trail") == 0) {
	 trail = 1;

      } else if(strcmp(argv[i], "-verb") == 0) {
	 TEST = 1;

      } else if(strcmp(argv[i], "-VERB") == 0) {
	 TEST = 2;

      } else if(strcmp(argv[i], "-bin") == 0) {
	 sscanf(argv[++i], "%d", &binfactor);
	 if(binfactor > MAXBIN || binfactor < 1) {
	    fprintf(stderr, "Error: binfactor %d > MAXBIN %d\n", binfactor, MAXBIN);
	    exit(1);
	 }

      } else if(strcmp(argv[i], "-border") == 0) {
	 sscanf(argv[++i], "%d", &x0border);
	 x1border = y0border = y1border = x0border;

      } else if(strcmp(argv[i], "-xlft") == 0) {
	 sscanf(argv[++i], "%d", &x0border);

      } else if(strcmp(argv[i], "-xrgt") == 0) {
	 sscanf(argv[++i], "%d", &x1border);

      } else if(strcmp(argv[i], "-ybot") == 0) {
	 sscanf(argv[++i], "%d", &y0border);

      } else if(strcmp(argv[i], "-ytop") == 0) {
	 sscanf(argv[++i], "%d", &y1border);

      } else if(strcmp(argv[i], "-nxsky") == 0) {
	 sscanf(argv[++i], "%d", &nxsky);

      } else if(strcmp(argv[i], "-nysky") == 0) {
	 sscanf(argv[++i], "%d", &nysky);

      } else if(strcmp(argv[i], "-rad") == 0) {
	 sscanf(argv[++i], "%lf", &srchrad);

      } else if(strcmp(argv[i], "-badata") == 0) {
	 sscanf(argv[++i], "%lf", &badata);

      } else if(strcmp(argv[i], "-sig") == 0) {
	 sscanf(argv[++i], "%lf", &srchsig);

      } else if(strcmp(argv[i], "-min") == 0) {
	 sscanf(argv[++i], "%lf", &srchcut);

      } else if(strcmp(argv[i], "-max") == 0) {
	 sscanf(argv[++i], "%lf", &srchmax);

      } else if(strcmp(argv[i], "-major") == 0) {
	 sscanf(argv[++i], "%lf", &major);

      } else if(strcmp(argv[i], "-minor") == 0) {
	 sscanf(argv[++i], "%lf", &minor);

      } else if(strcmp(argv[i], "-phi") == 0) {
	 sscanf(argv[++i], "%lf", &phi);
	 phi *= pi/180;	// convert to rad

      } else if(strcmp(argv[i], "-npar") == 0) {
	 sscanf(argv[++i], "%d", &nwpar);

      } else if(strcmp(argv[i], "-fwmax") == 0) {
	 sscanf(argv[++i], "%lf", &fwmax);

      } else if(strcmp(argv[i], "-fwmin") == 0) {
	 sscanf(argv[++i], "%lf", &fwmin);

      } else if(strcmp(argv[i], "-net") == 0) {
	 sscanf(argv[++i], "%lf", &netmin);

      } else if(strcmp(argv[i], "-move") == 0) {
	 sscanf(argv[++i], "%lf", &maxmove);

      } else if(strcmp(argv[i], "-okfit") == 0) {
	 sscanf(argv[++i], "%d", &okfit);

      } else if(strcmp(argv[i], "-chin") == 0) {
	 sscanf(argv[++i], "%lf", &chinmax);

      } else if(strcmp(argv[i], "-snr") == 0) {
	 sscanf(argv[++i], "%lf", &fitsig);

      } else if(strcmp(argv[i], "-eadu") == 0) {
	 if(strcmp(argv[++i], "auto") == 0) {
	    eaduin = -1;
	 } else {
	    sscanf(argv[i], "%lf", &eaduin);
	 }

      } else if(strcmp(argv[i], "-bias") == 0) {
	 sscanf(argv[++i], "%lf", &bias);

      } else if(strcmp(argv[i], "-sat") == 0) {
	 sscanf(argv[++i], "%lf", &sat);

      } else if(strcmp(argv[i], "-aprad") == 0) {
	 sscanf(argv[++i], "%d", &aprad);

      } else if(strcmp(argv[i], "-skyrad") == 0) {
	 sscanf(argv[++i], "%d", &skyrad);

      } else {
	 fprintf(stderr, "Unrecognized argument `%s'\n\n", argv[i]);
	 syntax(argv[0]);
	 exit(1);
      }   
   }

/* Are we being requested to do a 4 param fit? */
   if(major > 0 && minor > 0 && phi > -pi && nwpar >= 4) nwpar = 4;

/* Verify that we have shape and external positions for 2 param fit */
   if(nwpar <= 2 && 
      (major < 0 || minor < 0 || phi < -pi || objinput == NULL)) {
      fprintf(stderr, "Findobj: must specify input file, major,minor,phi for npar=2\n");
      exit(1);
   }

/* Sky is either constant or subdivided in both dimensions */
   if(nxsky == 0 || nysky == 0) nxsky = nysky = 0;

/* Open the output file */
   if(strcmp(fout, "-") == 0) {
      fp = stdout;
   } else if( (fp = fopen(fout, "w")) == NULL) {
      fprintf(stderr, "Cannot open output file '%s'\n", fout);
      syntax(argv[0]);
      exit(1);
   }

/* Open an input file? */
   if(objinput != NULL) {
      if(strcmp(objinput, "-") == 0) {
	 infp = stdin;
      } else if( (infp = fopen(objinput, "r")) == NULL) {
	 fprintf(stderr, "Cannot open input file '%s'\n", objinput);
	 syntax(argv[0]);
	 exit(1);
      }
   }

/* Does the image look like FITS? */
   if(!rdntf) {
      err = testfits_(imname, &mefits, &bitpix, strlen(imname)+1);
      rdntf = mefits != 1;
      bitpix = -32;
   }

/* Read the image */
   if(!rdntf) {		// Read FITS
      data = NULL;
      rfitsreal(&head, &nx, &ny, &data, imname);

   } else {		// Read NITF
      err = nitfread(imname, bitpix, &head, &nx, &ny, &data, TEST);
      if(err) {
	 fprintf(stderr, "Error from nitfread %d\n", err);
	 syntax(argv[0]);
	 exit(0);
      }
   }

/* Correct for the bias level */
   for(i=0; i<nx*ny; i++) {
      if(data[i] < 0 || data[i] >= 0) {	/* strip NaNs! */
	 if(data[i] != badata) data[i] -= bias;
      } else {
	 data[i] = badata;
      }
   }

/* Is there a flatfield? */
   if(flatfile != NULL) {
      err = rfitsreal(&headflat, &nxflat, &nyflat, &flat, flatfile);
      if(err) {
	 fprintf(stderr, "Error reading flat file %d\n", err);
	 exit(1);
      }
      if(nxflat != nx || nyflat != ny) {
	 fprintf(stderr, "Error: flat dims %d %d do not match image %d %d\n",
		 nxflat, nyflat, nx, ny);
	 exit(1);
      }
/* Normalize to max if no explicit value given */
      if(flatnorm == 0.0) {
	 for(i=0; i<nx*ny; i++) if(flat[i] > flatnorm) flatnorm = flat[i];
      }
/* Flatten the data */
      for(i=0; i<nx*ny; i++) data[i] /= flat[i] / flatnorm;
   }

/* Bin the image? */
   if(binfactor > 1) {
      for(j=0; j<ny/binfactor; j++) {
	 for(i=0; i<nx/binfactor; i++) {
	    z2 = 0.0;
	    for(jj=0; jj<binfactor; jj++) {
	       for(ii=0; ii<binfactor; ii++) {
		  z2 += data[ii+i*binfactor+(jj+j*binfactor)*nx];
	       }
	    }
	    data[i+j*(nx/binfactor)] = z2;
	 }
      }
      nx /= binfactor;
      ny /= binfactor;
      j = chfitshead(&i, head, "NAXIS1  ", "INTEGER", nx, 0.0);
      j = chfitshead(&i, head, "NAXIS2  ", "INTEGER", ny, 0.0);
   }

   NX = nx;
   if(ifitshead(head, "CNPIX1  ", &sx)) sx = 0;
   if(ifitshead(head, "CNPIX2  ", &sy)) sy = 0;

   if(TEST > 0) {
      fprintf(stderr, "Read %s\n  size %d x %d, offset (%d, %d),", 
	     imname, nx, ny, sx, sy);
      fprintf(stderr, "  data %6.1f %6.1f %6.1f\n", data[nx/2-1+ny/2*NX], 
	     data[nx/2+ny/2*NX], data[nx/2+1+ny/2*NX]);
   }

/* Get estimate of sky level and sky rms */
   skymed = (double *)calloc((nxsky+1)*(nysky+1), sizeof(double));
   skyrms = (double *)calloc((nxsky+1)*(nysky+1), sizeof(double));
   neadu = 0;
   eadu = 1.0;
   for(j=0; j<=nysky; j++) {
      for(i=0; i<=nxsky; i++) {
	 if(nxsky == 0 && nysky == 0) {		/* Constant sky */
	    imgstat(0, 0, nx-1, ny-1, NX, ny, data,
		    &skymed[0], &skyrms[0], &nskyok);
	 } else {				/* Chunked sky */
	    imgstat((int)((i-0.5)*(nx/nxsky)), (int)((j-0.5)*(ny/nysky)), 
		    (int)((i+0.5)*(nx/nxsky)), (int)((j+0.5)*(ny/nysky)),
		    NX, ny, data, 
		    &skymed[i+j*(nxsky+1)], &skyrms[i+j*(nxsky+1)], &nskyok);
	 }

/* Accumulate e/ADU as an average of mean/variance */
	 if(eaduin < 0) {
	    if(skyrms[i+j*(nxsky+1)]>0 && skymed[i+j*(nxsky+1)]>0) {
	       if(neadu == 0) eadu = 0.0;
	       eadu += skymed[i+j*(nxsky+1)] / 
		      (skyrms[i+j*(nxsky+1)]*skyrms[i+j*(nxsky+1)]);
	       neadu++;
	    }
	 } else {
	    eadu = eaduin;
	    neadu = 1;
	 }
	 if(TEST > 0) {
	    fprintf(stderr, "Sky= %7.1f  RMS= %7.1f  %4d  bias= %7.1f  e/ADU= %7.3f  sat= %7.1f\n", 
		    skymed[i+j*(nxsky+1)], skyrms[i+j*(nxsky+1)], 
		    nskyok, bias, eadu/neadu, sat);
	 }

      }
   }
   eadu /= neadu;

/* Are we going to write a residual image? */
   if(residfile != NULL) {
      resid = (float *)calloc(NX*ny, sizeof(float));
      for(ii=0; ii<NX*ny; ii++) {
	 resid[ii] = data[ii];
	 if(residsky) {
	    sky = skinterp(ii, nx, ny, nxsky, nysky, skymed, skyrms, &rms);
	    resid[ii] -= sky;
	 }
      }
   }

   if(!trail) {
      fprintf(fp, "#   x        y     peakval  skyval  peakfit   dpeak  skyfit     flux     dflux    major  minor    phi  err chi/N\n");
   } else {
      fprintf(fp, "#   x        y     peakval  skyval  peakfit   dpeak  skyfit     flux     dflux    trail   FWHM    phi  err chi/N\n");
   }

/* Find the objects */
   nfind = nfit = npix = 0;
   while(1) {

      if(infp != NULL) {	// Read an object from the input file
	 if(fgets(line, 10240, infp) == NULL) break;
	 if(line[0] == '#') continue;
	 if(nwpar >= 2) {
	    if(sscanf(line, "%lf %lf", &wpar.x0, &wpar.y0) != 2) {
	       fprintf(stderr, "Findobj: cannot read x,y from `%s'\n", line);
	       exit(1);
	    }
	 } else if(nwpar == 0) {	// No fit: strictly difference image
	    if(sscanf(line, "%lf %lf %lf", &wpar.x0, &wpar.y0, &wpar.peak) != 3) {
	       fprintf(stderr, "Findobj: cannot read x,y,peak from `%s'\n", line);
	       exit(1);
	    }
/* Fill in the rest of wpar since we're not going to fit anything */
	    wpar.beta4 = 1.0;
	    wpar.beta6 = 0.5;
	    if(!trail) {
	       err = dosxy(wpar);
	    } else {
	       wpar.sxy = phi;
	       wpar.sx2 = pow(2.3548/minor, 2.0);
	       wpar.sy2 = sqrt(major*major+minor*minor) - minor;
	       if(wpar.sy2 > 0) wpar.sy2 = log(0.5*wpar.sy2);
	    }
	 }
	 i = wpar.x0 + 0.5;
	 j = wpar.y0 + 0.5;

      } else {			// Step through the image
	 if(npix >= nx*ny) break;
	 j = npix / nx;
	 i = npix % nx;
	 npix++;
/* Skip pixels too close to the borders */
	 if(j < y0border || j > ny-1-y1border) continue;
	 if(i < x0border || i > nx-1-x1border) continue;

/* Pass the cut?  Weirdness with double negatives is in case of NaN */
	 sky = skinterp(npix, nx, ny, nxsky, nysky, skymed, skyrms, &rms);
	 if(!(data[i+j*NX] > sky+srchcut && data[i+j*NX] < sky+srchmax)) continue;
	 if(!(data[i+j*NX] > sky + srchsig*rms)) continue;

/* Maximum among closest neighbors? */
	 if(!(data[i+j*NX] >  data[i+1+j*NX]))   continue;
	 if(!(data[i+j*NX] >= data[i-1+j*NX]))   continue;
	 if(!(data[i+j*NX] >  data[i+(j+1)*NX])) continue;
	 if(!(data[i+j*NX] >= data[i+(j-1)*NX])) continue;

/* Maximum among all 8 neighbors? */
	 if(srchrad > 1) {
	    if(!(data[i+j*NX] >  data[i+1+(j+1)*NX])) continue;
	    if(!(data[i+j*NX] >= data[i-1+(j+1)*NX])) continue;
	    if(!(data[i+j*NX] >  data[i+1+(j-1)*NX])) continue;
	    if(!(data[i+j*NX] >= data[i-1+(j-1)*NX])) continue;
	 }
	 
/* Bigger search area?  One tie allowed. */
	 if(srchrad > sqrt(2.0)) {
	    nequal = 0;
	    idxeq = 0;
	    for(jj=-NINT(srchrad); jj<=NINT(srchrad); jj++) {
	       if(j+jj < 0 || j+jj >= ny) continue;
	       for(ii=-NINT(srchrad); ii<=NINT(srchrad); ii++) {
		  if(i+ii < 0 || i+ii >= nx) continue;
		  if(!(data[i+ii+(j+jj)*NX] != badata)) continue;
		  if(!(data[i+j*NX] >= data[i+ii+(j+jj)*NX])) {
		     nequal = 99;
		     break;
		  }
		  if(!(data[i+j*NX] != data[i+ii+(j+jj)*NX])) {
		     nequal++;
		     idxeq = ii<0 || (ii==0 && jj<0);
		  }
	       }
	       if(nequal > 2) break;
	    }

/* Local maximum?  Tie on the RHS? */
	    if(nequal > 2 || idxeq) continue;
	 }

/* OK, got one... analyze it */
	 wpar.x0 = i;
	 wpar.y0 = j;
      }
      wpar.major = major;
      wpar.minor = minor;
      wpar.phi = phi;

      if(TEST > 1) {
	 printf("\nNew object at x= %d y= %d  peak= %.1f\n", 
		i, j, data[i+j*NX]);
      }

      if(nwpar > 0) {	// Fit at all?
	 if(!trail) {
	    err = wpsf(nx, ny, data, eadu, sat, badata, aprad, skyrad,
		       nwpar, &wpar, &wflux);
	 } else {
	    err = trailpsf(nx, ny, data, eadu, sat, badata, aprad, skyrad,
			   nwpar, &wpar, &wflux);
	 }
	 nfit++;

/* Does this pass the remainder of the tests? */
	 failure = 0;
	 if(wpar.x0<x0border || wpar.x0>nx-1-x1border) failure |= FAIL_XBORDER;
	 if(wpar.y0<y0border || wpar.y0>ny-1-y1border) failure |= FAIL_YBORDER;
	 if(ABS(wpar.x0-(i+0.5)) > maxmove) failure |= FAIL_XMOVE;
	 if(ABS(wpar.y0-(j+0.5)) > maxmove) failure |= FAIL_YMOVE;
	 if(wpar.major > fwmax) failure |= FAIL_MAJBIG;
	 if(wpar.minor > fwmax) failure |= FAIL_MINBIG;
	 if(wpar.major < fwmin && !trail) failure |= FAIL_MAJSMALL;
	 if(wpar.minor < fwmin) failure |= FAIL_MINSMALL;
	 if(wpar.peak < netmin) failure |= FAIL_PEAK;
	 if(wpar.chin > chinmax) failure |= FAIL_CHIN;
	 if(wpar.chin == 0.0) failure |= FAIL_CHIN;
	 if(fitsig > 0) {
	    if(wflux.flux <= 0) failure |= FAIL_FLUXNEG;
	    if(wflux.dflux <= 0) failure |= FAIL_FERRNEG;
	    if(wflux.dflux > 0 && wflux.flux/wflux.dflux < fitsig) failure |= FAIL_SNR;
	 }

/* Debug info for all triggers for objects */
	 if(TEST > 1) {
	    printf("PSF: err= %d x0= %.2f y0= %.2f peak= %.1f sky= %.1f\n",
		   err, wpar.x0, wpar.y0, wpar.peak, wpar.sky);
	    printf("  maj= %.2f min= %.2f phi= %.3f chin= %.2f rms= %.2f resid= %.2f %.2f\n",
		   wpar.major, wpar.minor, wpar.phi*180/pi,
		   wpar.chin, wpar.rms, wpar.absresid, wpar.maxresid);
	    printf("Uncertainties: dx0= %.2f dy0= %.2f dpeak= %.1f dsky= %.1f\n",
		   wpar.dx0, wpar.dy0, wpar.dpeak, wpar.dsky);
	    printf("  dsx2= %.2f dsxy= %.2f dsy2= %.3f\n",
		   wpar.dsx2, wpar.dsxy, wpar.dsy2);
	    printf("Apflux: flux= %.1f dflux= %.1f sky=%.1f dsky= %.1f skyrms= %.1f\n",
		   wflux.flux, wflux.dflux, wflux.sky, wflux.dsky, 
		   wflux.skyrms);
	    if(failure) {
	       printf("FitFailCode= 0x%04x  ", failure);
	       failcode(failure);
	    }
	 }

	 if(err && okfit) failure |= FAIL_OKFIT;
	 if(failure) continue;

/* Tell us about a real object */
	 fprintf(fp, "%8.2f %8.2f  %7.1f %7.1f  %7.1f %7.1f %7.1f  %9.1f %8.1f  %6.2f %6.2f %6.1f  %2d %6.2f\n",
		 wpar.x0, wpar.y0, data[i+j*NX], wflux.sky, 
		 wpar.peak, wpar.dpeak, wpar.sky, wflux.flux, wflux.dflux, 
		 wpar.major, wpar.minor, wpar.phi*180/pi, err, wpar.chin);
	 nfind++;
      }

/* Subtract it from the image or residual image */
      if(residfile != NULL || subtract) {
	 rsub = sqrt(wpar.major*wpar.major+wpar.minor*wpar.minor) * RESIDRAD + 0.5;
	 if(TEST > 1) {
	    printf(" Subtract box x= %d %d y= %d %d\n",
		   i-rsub,i+rsub, j-rsub, j+rsub);
	 }
	 for(jj=j-rsub; jj<=j+rsub; jj++) {
	    if(jj < 0 || jj > ny-1) continue;
	    for(ii=i-rsub; ii<=i+rsub; ii++) {
	       if(ii < 0 || ii > nx-1) continue;
	       if(residfile != NULL) {
		  if(!trail) {
		     resid[ii+jj*NX] -= wauss(ii+0.5, jj+0.5, &wpar, &z2);
		  } else {
		     resid[ii+jj*NX] -= wtrail(ii+0.5, jj+0.5, &wpar, &z2);
		  }
	       }
	       if(subtract) {
		  if(!trail) {
		     data[ii+jj*NX] -= wauss(ii+0.5, jj+0.5, &wpar, &z2);
		  } else {
		     data[ii+jj*NX] -= wtrail(ii+0.5, jj+0.5, &wpar, &z2);
		  }
	       }
	    }
	 }
      }
   }

   if(TEST > 0) printf("Nfit= %d  Nfind= %d\n", nfit, nfind);

/* Write a residual image */
   if(residfile != NULL) wfitsreal(head, resid, residfile);

   exit(0);
}

/* Compute median and median RMS for an array */
int imgstat(int x0, int y0, int x1, int y1, int NX, int NY, 
	    float *data, double *med, double *rms, int *nok)
{
   int i, j, k, m, n, nsrt, err;
   double buf[MAXSRT];
   *nok = 0;
   nsrt = MIN((x1-x0+1)*(y1-y0+1), MAXSRT);

//   fprintf(stderr, "%d %d %d %d %d\n", x0, x1, y0, y1, nsrt);
   for(m=n=0; m<nsrt; m++) {
      k = (m/(double)nsrt) * (x1-x0+1)*(y1-y0+1);
      i = (k%(x1-x0+1)) + x0;
      j = (k/(x1-x0+1)) + y0;
      if(i<0 || i>=NX || j<0 || j>=NY) continue;
      buf[n++] = data[i+j*NX];
   }
   if(n < 1) {
      *med = *rms = 0.0;
      return(-1);
   }
   if(n < 2) {
      *med = buf[0];
      *rms = 0;
      return(0);
   }
   if( (err = qsort8(n, buf)) ) return(err);
   *med = 0.5 * (buf[(n-1)/2] + buf[n/2]);
   *rms = *med - buf[n/4];
   *nok = n;
   return(0);
}

/* Evaluate a 4 point linear interpolation of sky at image index k */
double skinterp(int k, int nx, int ny, int nxsky, int nysky, double *skymed,
		double *skyrms, double *rms)
{
   int is, js;
   double u, v, sky;
   if(nxsky == 0 && nysky == 0) {
      sky = skymed[0];
      *rms = skyrms[0];
   } else {
      is = (k%nx) / (nx/nxsky);
      js = (k/nx) / (ny/nysky);
/* Edge squares are calculated between 0:0.5, not 0:1 => extrapolation */
      if(is == 0) {
	 u = ((double)(k%nx) / (nx/nxsky) - is - 0.25) * (4.0/3.0);
      } else if(is == nxsky-1) {
	 u = ((double)(k%nx) / (nx/nxsky) - is) * (4.0/3.0);
      } else {
	 u = (double)(k%nx) / (nx/nxsky) - is;
      }
      if(js == 0) {
	 v = ((double)(k/nx) / (ny/nysky) - js - 0.25) * (4.0/3.0);
      } else if(js == nysky-1) {
	 v = ((double)(k/nx) / (ny/nysky) - js) * (4.0/3.0);
      } else {
	 v = (double)(k/nx) / (ny/nysky) - js;
      }
      sky = (1-u) * (1-v) * skymed[is+js*(nxsky+1)] +
	 u * (1-v) * skymed[is+1+js*(nxsky+1)] +
	 (1-u) * v * skymed[is+(js+1)*(nxsky+1)] +
	 u * v * skymed[is+1+(js+1)*(nxsky+1)];
      *rms = (1-u) * (1-v) * skyrms[is+js*(nxsky+1)] +
	 u * (1-v) * skyrms[is+1+js*(nxsky+1)] +
	 (1-u) * v * skyrms[is+(js+1)*(nxsky+1)] +
	 u * v * skyrms[is+1+(js+1)*(nxsky+1)];
   }
   return(sky);
}

/* Print reasons for failure */
int failcode(int failure)
{
   if(failure & FAIL_XBORDER) printf(" XBORDER");
   if(failure & FAIL_YBORDER) printf(" YBORDER");
   if(failure & FAIL_XMOVE) printf(" XMOVE");
   if(failure & FAIL_YMOVE) printf(" YMOVE");
   if(failure & FAIL_MAJBIG) printf(" MAJBIG");
   if(failure & FAIL_MINBIG) printf(" MINBIG");
   if(failure & FAIL_MAJSMALL) printf(" MAJSMALL");
   if(failure & FAIL_MINSMALL) printf(" MINSMALL");
   if(failure & FAIL_PEAK) printf(" PEAK");
   if(failure & FAIL_CHIN) printf(" CHIN");
   if(failure & FAIL_FLUXNEG) printf(" FLUXNEG");
   if(failure & FAIL_FERRNEG) printf(" FERRNEG");
   if(failure & FAIL_SNR) printf(" SNR");
   if(failure & FAIL_OKFIT) printf(" OKFIT");
   if(failure) printf("\n");
   return(0);
}

int syntax(char *prog)
{
   printf("%s image [options]\n\n", prog);
   printf("Finds objects in an image, reports on properties\n");
   printf("[options] include:\n\n");
   printf("  -out fname   Output file name (stdout)\n");
   printf("  -obj fname   Use x,y from fname for triggers (NULL)\n");
   printf("  -resid fname Write a residual FITS image (NULL)\n");
   printf("  -subtract    Subtract fits successively from image (default no)\n");
   printf("  -residsky    Subtract sky model from residual image (default no)\n");
   printf("  -trail       Fit a trailed PSF (default no)\n");
   printf("  -verb        Diagnostic verbosity\n");
   printf("  -VERB        More diagnostic verbosity\n\n");
   printf("  -bin N       Bin input image by a factor NxN (%d)\n", binfactor);
   printf("Trigger parameters:\n");
   printf("  -border N    Image border to avoid on all sides (%d)\n", x0border);
   printf("  -xlft N      Left x image border to avoid (%d)\n", x0border);
   printf("  -xrgt N      Right x image border to avoid (%d)\n", x1border);
   printf("  -ybot N      Bottom y image border to avoid (%d)\n", y0border);
   printf("  -ytop N      Top y image border to avoid (%d)\n", y1border);
   printf("  -badata B    Ignore data values of B (%.1f)\n", badata);
   printf("  -rad R       Local max required over radius r to trigger (%.1f)\n", srchrad);
   printf("  -sig S       Minimum value = sky+sig*rms to trigger (%.1f)\n", srchsig);
   printf("  -min Z       Minimum value above sky to trigger (%.1f)\n", srchcut);
   printf("  -max Z       Maximum peak to trigger (%.1f)\n\n", srchmax);
   printf("Fit acceptance parameters:\n");
   printf("  -fwmax W     Maximum FW to keep (%.2f)\n", fwmax);
   printf("  -fwmin W     Minimum FW to keep (%.2f)\n", fwmin);
   printf("  -net Z       Minimum peak-sky to keep (%.1f)\n", netmin);
   printf("  -move M      Max offset between peak and fit (%.1f)\n", maxmove);
   printf("  -chin C      Maximum chi/N to keep (%.2f)\n", chinmax);
   printf("  -snr S       Minimum flux SNR to keep (%.1f)\n", fitsig);
   printf("  -okfit I     Insist (0/1) on successful profile fit (%d)\n\n", okfit);
   printf("The fit is normally 7 parameters (position, peak, sky, PSF\n");
   printf("but if major, minor, and phi are set a 4 param fit is used:\n");
   printf("  -major M     Set major axis M for 4 param fit [pix]\n");
   printf("  -minor M     Set minor axis M for 4 param fit [pix]\n");
   printf("  -phi P       Set major axis angle to P for 4 param fit [deg]\n");
   printf("  -npar N      How many parameter fit (2/4/7)? (%d)\n\n", nwpar);
   printf("Noise, signal, and aperture parameters:\n");
   printf("  -eadu E      Electrons per ADU to estimate noise (%.1f)\n", eadu);
   printf("               if E is 'auto' eadu is set to mean/variance\n");
   printf("  -bias B      Image background exceeds sky by B (%.1f)\n", bias);
   printf("  -sat S       Saturation level is S (%.1f)\n", sat);
   printf("  -aprad R     Aperture radius for photometry (%d)\n", aprad);
   printf("  -skyrad R    Sky radius for photometry (%d)\n", skyrad);
   printf("  -nxsky N     Subdivide image by N in x for sky (%d)\n", nxsky);
   printf("  -nysky N     Subdivide image by N in y for sky (%d)\n", nysky);
   printf("  -flat F      Read flat file and divide by it\n");
   printf("  -flatnorm N  Divide flat file by N before application\n");
//   printf("\n  -nitf        Input is NITF, not FITS\n");
   return(0);
}
