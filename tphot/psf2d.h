/* psf2d.h -- header file for psf2d.c: fit a 2D Waussian profile */
/* 120929 v1.0 John Tonry */

#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define ABS(a) (((a) > 0) ? (a) : -(a))
#define NINT(x) (x<0?(int)((x)-0.5):(int)((x)+0.5))

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

/* Return error code 0: no error, 
 *                  -1: off image, 
 *                   1: fit didn't converge,
 *                N*10: N zero'ed pixels within FWHM
 */
int wpsf(int nx, int ny, float *a, double eadu, double sat, double badata,
          int aprad, int skyrad, int nwpar, 
          PSF2D_PARAM *wpar, PSF2D_FLUX *flux);

/* Fit a trailed PSF */
int trailpsf(int nx, int ny, float *a, double eadu, double sat, double badata,
	     int aprad, int skyrad, int nwpar, 
	     PSF2D_PARAM *wpar, PSF2D_FLUX *flux);
