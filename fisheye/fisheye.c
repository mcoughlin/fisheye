/* Program to fit fisheye lens parameters from a table of x,y,RA,Dec */
/* v1.01 140227 start evolving... */
/* v1.0  140206 John Tonry */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

#define ABS(a) (((a) > 0) ? (a) : -(a))
#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define NINT(a)  (((a) > 0) ? (int)((a)+0.5) : (int)((a)-0.5))

/* fitlm control parameters */
#define FITLM_VERBOSE 0x01	/* Verbose? */
#define FITLM_DERIV   0x02	/* Analytic first derivatives? */
#define FITLM_LINEAR  0x04	/* Linear function? */
#define FITLM_WGT     0x08	/* Data point weights provided? */

#define MAXPT 10000	/* Max number of data points */
#define MAXPAR 9	/* Number of parameters we're fitting */
#define MAXNEIGH 32	/* Max number of neighbors */

// #define FITXY;		/* Minimize x,y offsets instead of sky angles */

static int TEST=0;
static double LAT;		/* Make latitude a global variable */

typedef struct xyad {
   double x;			/* x posn on fisheye */
   double y;			/* y posn on fisheye */
   double ha;			/* HA or HA on sky */
   double dec;			/* Dec on sky */
   double xmap;			/* x posn calculated from HA,Dec */
   double ymap;			/* y posn calculated from HA,Dec */
   double thmap;		/* star theta posn calculated from HA,Dec */
   double hamap;		/* star HA posn calculated from x,y */
   double decmap;		/* star Dec posn calculated from x,y */
   double azimuth;		/* star azimuth calculated from HA,Dec */
   double alt;			/* star altitude calculated from HA,Dec */
   double jacobian;		/* star Jacobian [steradian/pix] */
   double dang;			/* angle between cat and calc sky posn */
   double dxnear;		/* <xmap-x> of neighbors */
   double dynear;		/* <ymap-y> of neighbors */
   struct xyad **near;		/* list of nearest neighbors */
   double *dist;		/* list of nearest neighbors */
   double resid;		/* residual wrt nearest neighbors */
} XYAD;


/* Fit a model to data, minimizing chi^2 using Levenberg-Marquardt */
typedef int LMFUNC(int kpix, void *dataptr, void *usr1, void *usr2, 
		   double *dataval, double *wgt, 
		   int npar, double *par, double *model, double *deriv);
int fitlm(int ndata, void *dataptr, void *usr1, void *usr2, LMFUNC func,
	  int npar, double *par, int *usepar, 
	  double *chi, double *cov, int *ctrl);
int printmat(int n, int *usepar, double *a, double *y, double chi);
/* Give user access to QUITFRAC LAMINIT, etc.  Not for the novice! */
int fitlm_param(int maxiter, int laminit, double quit, double parvary);

int fisheye(int kpix, void *dataptr, void *usr1, void *usr2, 
	    double *dataval, double *wgt, 
	    int npar, double *par, double *model, double *deriv);

/* Angle between HA1, Dec1 and HA2, Dec2 */
double skyangle(double ha1, double dec1, double ha2, double dec2);
/* Convert HA,Dec to az, alt, rotation */
void HD_to_AZ(double ha, double dec, double lat, 

	      double *alt, double *az, double *rot);
/* Map an image position back to the sky */
int fishsky(double ximg, double yimg, double *par, 
	    double *ha, double *dec, double *az, double *alt, double *jacobian);
/* Map a sky coordsinate to image position */
int skyfish(double ha, double dec, double *par, 
	    double *ximg, double *yimg, double *zz);

/* Solve for pole, azimuth, and scale from three points */
int threept(XYAD *xyad1, XYAD *xyad2, XYAD *xyad3, double cx, double cy, 
	    double *par, int *parity);
/* Project all triangles and estimate image parameters */
void project_triangles(int npt, XYAD *xyad, double cx, double cy, int parity);
/* Invert a matrix, return the determinant */
int invert(int n, double *a, double *det);

/* Find nearest neighbors to each point */
int nearest(int npt, XYAD *xyad, int nneigh);
/* Evaluate residual of this point */
void resid(int nneigh, XYAD *xyad);
/* heapsort */
void heapsort(int n, double *data, void **idx);
void siftdown(int root, int bottom, double *data, void **idx);
double median(int n, double *data);

void syntax(char *prog);

/* Create some data and fit a Gaussian to it */
int main(int argc, char **argv)
{
   int i, iter, err, npt, parity, eval, nneigh, fixctr, fixdist, maxiter;
   XYAD xyad[MAXPT];
   int npar, usepar[MAXPAR], ctrl;
   double par[MAXPAR], cov[MAXPAR*MAXPAR], chi;
   double cx, cy;
   double ha, dec, az, scale, quad, cube;
   double ang, x, y, z, u, v, r;
   double srtbuf[MAXPT];
   double dr=atan(1.0)/45.0, pi=4*atan(1.0);
   char *infile, *outfile, line[256];
   FILE *fp=stdin, *fpo=stdout;

/* Parse the arguments */
   LAT = 19.5362*dr;	/* Fisheye latitude required for zenith distance */

   ha = -10.0;		/* Pole of fisheye on sky */
   dec = -10.0;		/* Pole of fisheye on sky */
   az = -10.0;		/* Orientation of fisheye */
   scale = -1.0;	/* Scale (rad/pix) */
   quad = -0.018;	/* Quadratic correction to scale (rad^-1) */
   cube = -0.030;	/* Cubic correction to scale (rad^-2) */
   cx = 1421.0;		/* Center of fisheye */
   cy = 975.0;		/* Center of fisheye */
   fixctr = 0;		/* Fix cx,cy */
   fixdist = 0;		/* Fix distortion quad and cube */
   parity = +1;		/* HA (+1) vs RA (-1) */
   infile = "-";
   outfile = NULL;
   npt = 0;
   eval = 0;		/* Evaluate, not fit */
   nneigh = 10;		/* Number of neighbors from which to get residuals */
   maxiter = 20;	/* Max number of iterations */

/* Parse the arguments */
   for(i=1; i<argc; i++) {
      if (strcmp(argv[i], "-in") == 0) {		/* Input file */
	 infile = argv[++i];

      } else if (strcmp(argv[i], "-cx") == 0) {		/* image ctr in x */
	 sscanf(argv[++i], "%lf", &cx);

      } else if (strcmp(argv[i], "-cy") == 0) {		/* image ctr in y */
	 sscanf(argv[++i], "%lf", &cy);

      } else if (strcmp(argv[i], "-fixctr") == 0) {	/* fix cx,cy */
	 fixctr = 1;

      } else if (strcmp(argv[i], "-fixdist") == 0) {	/* fix quad,cube */
	 fixdist = 1;

      } else if (strcmp(argv[i], "-maxiter") == 0) {	/* max fit iterations */
	 sscanf(argv[++i], "%d", &maxiter);

      } else if (strcmp(argv[i], "-ha") == 0) {		/* Center HA [deg] */
	 sscanf(argv[++i], "%lf", &ha);
	 ha *= dr;

      } else if (strcmp(argv[i], "-ra") == 0) {		/* Center RA [deg] */
	 sscanf(argv[++i], "%lf", &ha);
	 parity = -1;
	 ha *= parity * dr;

      } else if (strcmp(argv[i], "-dec") == 0) {	/* Center Dec [deg] */
	 sscanf(argv[++i], "%lf", &dec);
	 dec *= dr;

      } else if (strcmp(argv[i], "-az") == 0) {		/* Fisheye rot [deg] */
	 sscanf(argv[++i], "%lf", &az);
	 az *= dr;

      } else if (strcmp(argv[i], "-lat") == 0) {	/* Obs latitude [deg] */
	 sscanf(argv[++i], "%lf", &LAT);
	 LAT *= dr;

      } else if (strcmp(argv[i], "-scale") == 0) {	/* scale [deg/pix] */
	 sscanf(argv[++i], "%lf", &scale);
	 scale *= dr;

      } else if (strcmp(argv[i], "-quad") == 0) {	/* quad */
	 sscanf(argv[++i], "%lf", &quad);

      } else if (strcmp(argv[i], "-cube") == 0) {	/* cubic term */
	 sscanf(argv[++i], "%lf", &cube);

      } else if (strcmp(argv[i], "-out") == 0) {	/* output file */
	 outfile = argv[++i];

      } else if (strcmp(argv[i], "-xy2sky") == 0) {	/* evaluate only */
	 eval = 1;

      } else if (strcmp(argv[i], "-sky2xy") == 0) {	/* evaluate only */
	 eval = 2;

      } else if (strcmp(argv[i], "-verb") == 0) {	/* moderate verbosity */
	 TEST = 1;

      } else if (strcmp(argv[i], "-VERB") == 0) {	/* high verbosity */
	 TEST = 2;

/* Deprecated options */
      } else if (strcmp(argv[i], "-parity") == 0) {	/* parity = +/-1 */
	 sscanf(argv[++i], "%d", &parity);

      } else if (strcmp(argv[i], "-eval") == 0) {	/* evaluate only */
	 eval = 2;

      } else if (i == 1) {				/* First arg is input */
	 infile = argv[i];

      } else {
	 fprintf(stderr, "Unknown arg `%s'\n", argv[i]);
	 syntax(argv[0]);
	 exit(1);
      }
   }

/* Did we get a specification for pole, azimuth, and scale? */
   if((ha<-2*pi || dec<-2*pi || az<-2*pi || scale<=0.0) && eval) {
      fprintf(stderr, "Error: cannot evaluate without -ha -dec -az -scale\n");
      exit(1);
   }

/* Initialize parameter array */
   par[0] = ha;
   par[1] = dec;
   par[2] = az;
   par[3] = log(MAX(1e-10, scale));
   par[4] = cx;		// center in x, note covariance with ha,dec
   par[5] = cy;		// center in y, note covariance with ha,dec
   par[6] = quad;	// quadratic component
   par[7] = cube;	// cubic component
   par[8] = parity;

   if(TEST > 1) {
      printf("Initial parameters: %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %d\n", 
	     par[0], par[1], par[2], par[3], par[4], par[5], par[6], par[7], (int)par[8]);
   }

/* Open input file */
   if(strcmp(infile, "-")) {
      if( (fp=fopen(infile, "r")) == NULL) {
	 fprintf(stderr, "Cannot open `%s' for reading\n", infile);
	 exit(1);
      }
   }

/* Open output file */
   if(outfile != NULL) {
      if(strcmp(outfile, "-")) {
	 if( (fpo=fopen(outfile, "w")) == NULL) {
	    fprintf(stderr, "Cannot open `%s' for writing\n", outfile);
	    exit(1);
	 }
      }
   }

/* Read the data */
   for(npt=0;   ; npt++) {
      if(npt == MAXPT || fgets(line, 256, fp) == NULL) break;
/* Normal use to fit data? */
      if(!eval) {
	 if( sscanf(line, "%lf %lf %lf %lf", &xyad[npt].x, &xyad[npt].y, 
		    &xyad[npt].ha, &xyad[npt].dec) != 4) {
	    fprintf(stderr, "Error reading `%s' at line %d\n", infile, npt+1);
	    exit(1);
	 }
	 if(TEST > 1) {
	    printf("%8.2f %8.2f %8.2f %8.2f\n", xyad[npt].x, xyad[npt].y, 
		   xyad[npt].ha, xyad[npt].dec);
	 }
	 xyad[npt].ha *= parity * dr;	// Angles are internally parity +1
	 xyad[npt].dec *= dr;

/*********************** EVALUATION *****************************/
/* Just evaluate x,y to HA,Dec do not save in xyad */
      } else {
	 if( sscanf(line, "%lf %lf", &x, &y) != 2) {
	    fprintf(stderr, "Error reading `%s' at line %d\n", infile, npt+1);
	    exit(1);
	 }
	 if(eval == 1) {	/* x,y to ha,dec */
	    fishsky(x, y, par, &u, &v, &az, &quad, &cube);
	    fprintf(fpo, "%7.2f %7.2f  %8.3f %8.3f\n", x, y, u*parity/dr, v/dr);
	 } else {		/* ha,dec to x,y */
	    x *= parity * dr;	// Angles are internally parity +1
	    y *= dr;
	    skyfish(x, y, par, &u, &v, &z);
	    fprintf(fpo, "%7.2f %7.2f  %8.3f %8.3f\n", u, v, x*parity/dr, y/dr);
	 }
      }
/*********************** EVALUATION *****************************/
   }


/*********************** EVALUATION *****************************/
/* If just map evaluation, we're done, exit */
   if(eval) exit(0);
/*********************** EVALUATION *****************************/


/*********************** BOOTSTRAP ******************************/
/* pole, azimuth or scale uninitialized -> three point estimate and exit */
   if(ha<-2*pi || dec<-2*pi || az<-2*pi || scale<=0.0) {
      project_triangles(npt, xyad, cx, cy, parity);
      exit(0);
   }
/*********************** BOOTSTRAP ******************************/


/* How many parameters does the function have? */
   npar = MAXPAR;
   for(i=0; i<npar; i++) usepar[i] = 1; // Fit all parameters
   if(fixctr) usepar[4] = usepar[5] = 0;	// Leave cx,cy fixed
   if(fixdist) usepar[6] = usepar[7] = 0;	// Leave quad,cube fixed
   usepar[8] = 0;			// Carry parity, but never fit

/* Find the nearest points (excluding self) for residual calculation */
   if(nneigh >= npt/2) nneigh = npt/2 - 1;
   nearest(npt, xyad, nneigh);

/* Fit the parameters to the data, no weights */
   for(iter=0; iter<maxiter; iter++) {
      ctrl = TEST > 0 ? FITLM_VERBOSE : 0;
      fitlm_param(20/*maxiter*/, 0/*lam*/, 1e-8/*quit*/, 1e-5/*parvary*/);
      err = fitlm(npt, (void *)xyad, NULL, NULL, fisheye,
		  npar, par, usepar, &chi, cov, &ctrl);

/* Correct angles into 0-360 */
      par[0] = fmod(par[0]+8*pi, 2*pi);
      par[2] = fmod(par[2]+8*pi, 2*pi);

      if(TEST > 0) {
	 printf("  Chi-scaled error/covariance matrix and params:\n");
	 printmat(npar, usepar, cov, par, chi);
	 printf("%cA= %6.2f Dec= %6.2f Az= %6.2f Scale= %7.5f  chi= %.5f, lam= %d, niter= %d\n", 
		parity>0?'H':'R', parity*par[0]/dr,
		par[1]/dr, par[2]/dr, exp(par[3])/dr,
	     chi, (ctrl&0xff)-128, ctrl/256);
      }
      if(!err) break;
   }

/* Compute the mapped coords for each point */
   for(i=0; i<npt; i++) {
      skyfish(xyad[i].ha, xyad[i].dec, par, 
	      &xyad[i].xmap, &xyad[i].ymap, &xyad[i].thmap);
      fishsky(xyad[i].x, xyad[i].y, par, 
	      &xyad[i].hamap, &xyad[i].decmap, 
	      &xyad[i].azimuth, &xyad[i].alt, &xyad[i].jacobian);
      xyad[i].dang = skyangle(xyad[i].hamap, xyad[i].decmap, 
			      xyad[i].ha, xyad[i].dec);
      srtbuf[i] = xyad[i].dang;
   }

/* Evaluate the median vector field of residuals for each point */
   for(i=0; i<npt; i++) resid(nneigh, xyad+i);

/* What is the median angle difference? */
   ang = median(npt, srtbuf);

/* Tell us about the fit results on stdout */
   printf("%cA= %6.2f Dec= %6.2f Az= %6.2f Scale= %7.5f quad= %8.5f cube= %8.5f cx= %6.1f cy= %6.1f  <dth>= %.3f %.3f N= %d lam= %d niter= %d\n", 
	  parity>0?'H':'R', parity*par[0]/dr,
	  par[1]/dr, par[2]/dr, exp(par[3])/dr, par[6], par[7], par[4], par[5],
	  sqrt(chi/npt)/dr, ang/dr, npt, (ctrl&0xff)-128, 20*(iter-1)+ctrl/256);

/* Write the results to an output file if requested */
   if(outfile != NULL) {
      fprintf(fpo, "# %cA= %7.3f Dec= %7.3f Az= %7.3f Scale= %7.5f quad= %8.5f cube= %8.5f cx= %6.1f cy= %6.1f  <dth>= %.3f %.3f N= %d lam= %d niter= %d\n", 
	      parity>0?'H':'R', parity*par[0]/dr,
	      par[1]/dr, par[2]/dr, exp(par[3])/dr, par[6], par[7], par[4], par[5],
	      sqrt(chi/npt)/dr, ang/dr, npt, (ctrl&0xff)-128, 20*(iter-1)+ctrl/256);
      fprintf(fpo, "#   x       y        %cA     Dec      %cAcalc  Decalc dAngle   Azi     Alt  ompix    xcalc   ycalc    r      theta   <dx>  <dy> resid\n", 
	      parity>0?'H':'R', parity>0?'H':'R');
      for(i=0; i<npt; i++) {
	 r = sqrt(pow(xyad[i].x-cx,2.0) + pow(xyad[i].y-cy,2.0));
	 u = parity*xyad[i].ha/dr;
	 if(u > 180.0) u -= 360.0;
	 if(u < -180.0) u += 360.0;
	 fprintf(fpo, "%7.1f %7.1f  %8.3f %7.3f  %8.3f %7.3f %5.3f %7.2f %7.2f %5.0f  %7.1f %7.1f %7.2f %8.5f %5.2f %5.2f %5.2f\n", 
		 xyad[i].x, xyad[i].y, u, xyad[i].dec/dr,
		 parity*xyad[i].hamap/dr, xyad[i].decmap/dr, xyad[i].dang/dr, 
		 xyad[i].azimuth/dr, xyad[i].alt/dr,
		 3600*3600*xyad[i].jacobian/(dr*dr),
		 xyad[i].xmap, xyad[i].ymap, r, xyad[i].thmap,
		 xyad[i].dxnear, xyad[i].dynear, xyad[i].resid);
      }
   }

   exit(0);
}


int fisheye(int kpix, void *dataptr, void *usr1, void *usr2, 
	     double *dataval, double *wgt,
	     int npar, double *par, double *model, double *deriv)
{
   XYAD *xyad=(XYAD *)dataptr;
   double ximg, yimg, ha, dec, xp, yp, az, alt, jac;

   ximg = xyad[kpix].x;
   yimg = xyad[kpix].y;
   ha   = xyad[kpix].ha;
   dec  = xyad[kpix].dec;

   *dataval = 0.0;

#ifdef FITXY
   skyfish(ha, dec, par, &xp, &yp);
   *model = sqrt((xp-ximg)*(xp-ximg) + (yp-yimg)*(yp-yimg));
# else
   fishsky(ximg, yimg, par, &xp, &yp, &az, &alt, &jac);
   *model = skyangle(xp, yp, ha, dec);
#endif

   if(TEST > 1) {
#ifdef FITXY
      printf("%7.1f  %7.1f %7.3f %7.3f  %8.3f %8.3f  %8.3f %8.3f\n",
	     ximg, yimg, ha, dec, xp, yp, ximg-xp, yimg-yp);
# else
      printf("%7.1f  %7.1f %7.3f %7.3f  %8.3f %8.3f  %8.3f %8.3f\n",
	     ximg, yimg, ha, dec, xp, yp, ha-xp, dec-yp);
#endif
   }

   return(0);
}

/* Map a sky coordinate to image position */
int skyfish(double ha, double dec, double *par, 
	    double *ximg, double *yimg, double *theta)
{
   double a0=par[0], d0=par[1], az=par[2], scale=exp(par[3]);
   double cx=par[4], cy=par[5], quad=par[6], cube=par[7];
   double x, y, z, phi, r;

/* Rotated coords */
   x =  cos(az)*(sin(d0)*cos(dec)*cos(a0-ha)-cos(d0)*sin(dec)) + 
        sin(az)*cos(dec)*sin(ha-a0);
   y = -sin(az)*(sin(d0)*cos(dec)*cos(a0-ha)-cos(d0)*sin(dec)) + 
        cos(az)*cos(dec)*sin(ha-a0);
   z =  cos(d0)*cos(dec)*cos(ha-a0) + sin(d0)*sin(dec);

/* Project onto image */
   phi = atan2(y, x);
   *theta = acos(z);
   r = (*theta) / scale * (1 + (*theta)*(quad + (*theta)*cube));
   *ximg = cx + r * cos(phi);
   *yimg = cy + r * sin(phi);

// Rotate by HA around z to bring HA,Dec to x-z plane
// Rotate by coDec around y to bring HA,Dec to N pole (N is in -u dir)
// Rotate by az around z to orient N in +u dir
//
//  x     caz saz  0      sd0  0 -cd0      ca0 sa0  0       ca*cd
//  y =  -saz caz  0  *    0   1   0   *  -sa0 ca0  0   *   sa*cd
//  z     0   0    1      cd0  0  sd0       0   0   1          sd
//
//         caz saz  0     sd0*cd*c(a0-a) - cd0*sd
//    =   -saz caz  0  *      cd*s(a-a0)
//         0   0    1     cd0*cd*c(a0-a) + sd0*sd
//
//         caz*(sd0*cd*c(a0-a) - cd0*sd) + saz*cd*s(a-a0)
//    =   -saz*(sd0*cd*c(a0-a) - cd0*sd) + caz*cd*s(a-a0)
//         cd0*cd*c(a0-a) + sd0*sd
//
//  phi = atan2(y,x)
//    r ~ A*sqrt(2) * sin(acos(z)/2) = A * sqrt(1-z)
//
//  ximg = cx +  r * cos(phi)
//  yimg = cy +  r * sin(phi)

   return(0);
}

/* Map an image position back to the sky */
int fishsky(double ximg, double yimg, double *par, 
	    double *ha, double *dec, double *azi, double *alt, double *jacobian)
{
   double a0=par[0], d0=par[1], az=par[2], scale=exp(par[3]);
   double cx=par[4], cy=par[5], quad=par[6], cube=par[7];
   double x, y, z, phi, r, u, v, w, theta, rot;

   r = sqrt((ximg-cx)*(ximg-cx)+(yimg-cy)*(yimg-cy));
   phi = atan2(yimg-cy, ximg-cx);

/* Rotated coords */
//   r*scale = theta + quad * theta^2 + cube * theta^3;
   if(ABS(quad) > 1e-6) {
      theta = (sqrt(1+4*quad*r*scale)-1) / (2*quad);
      theta += (r*scale - theta*(1+theta*(quad+theta*cube))) /
	       (1+theta*(2*quad+3*cube*theta));
      theta += (r*scale - theta*(1+theta*(quad+theta*cube))) /
	       (1+theta*(2*quad+3*cube*theta));
      theta += (r*scale - theta*(1+theta*(quad+theta*cube))) /
	       (1+theta*(2*quad+3*cube*theta));
   } else {
      theta = r*scale;
   }
   z = cos(theta);
//   z = 1 - r*r * (scale*scale);
   x = sqrt(1-z*z) * cos(phi);
   y = sqrt(1-z*z) * sin(phi);

/* Sky coords */
   u = (cos(az)*sin(d0)*cos(a0)-sin(az)*sin(a0))*x +
      (-sin(az)*sin(d0)*cos(a0)-cos(az)*sin(a0))*y +
        cos(d0)*cos(a0)*z;
   v = (cos(az)*sin(d0)*sin(a0)+sin(az)*cos(a0))*x +
      (-sin(az)*sin(d0)*sin(a0)+cos(az)*cos(a0))*y +
        cos(d0)*sin(a0)*z;
   w = -cos(az)*cos(d0)*x + sin(az)*cos(d0)*y + sin(d0)*z;

   *ha = atan2(v, u);
   *dec = asin(w);

/* Zenith distance */
   HD_to_AZ(*ha, *dec, LAT, alt, azi, &rot);

/* Jacobian (rad^2 per pix^2) */
   if(r < 1e-6) {
      *jacobian = scale * scale;
   } else if(ABS(quad) > 1e-6) {
      *jacobian = scale * sin(theta) / r / sqrt(1+4*quad*r*scale);
   } else {
      *jacobian = scale * sin(theta) / r;
   }

//  ximg = cx +  r * cos(phi)
//  yimg = cy +  r * sin(phi)
//
//  phi = atan2(y,x)
//    r ~ A*sqrt(2) * sin(acos(z)/2) = A * sqrt(1-z)
//
//  x     czsdca-szsa   czsdsa+szca  -czcd     u
//  y =  -szsdca-czsa  -szsdsa+czca   szcd  *  v
//  z            cdca   cdsa            sd     w

   return(0);
}

/* Project all triangles and estimate image parameters */
void project_triangles(int npt, XYAD *xyad, double cx, double cy, int parity)
{
   int i, j, k, n, p, err;
   double *threepar, *srtbuf, par[9];
   double dr=atan(1.0)/45;
   double ha, dec, azi, alt, jac, u, r, ang;

   if(npt > 100) {
      fprintf(stderr, "Error: too many points for triangle estimates\n");
      exit(1);
   }
   n = 0;
   threepar = (double *)calloc(4*npt*(npt-1)*(npt-2)/6, sizeof(double));
   srtbuf = (double *)calloc(npt*(npt-1)*(npt-2)/6, sizeof(double));
   if(TEST > 1) {
      printf("  k   j   i      %cA       Dec      Az     Scale    uXv.z    xyz.r     dot\n", parity>0?'H':'R');
   }
   for(k=2; k<npt; k++) {
      for(j=1; j<k; j++) {
	 for(i=0; i<j; i++) {
	    if(TEST > 1) printf("%3d %3d %3d  ", k, j, i);
	    if( (err = threept(xyad+k, xyad+j, xyad+i, cx, cy, 
			       threepar+4*n, &p)) == 0) {
	       srtbuf[n++] = p;	// parity
	    }
	 }
      }
   }
/* Populate parameters with median results, report, exit */
   par[8] = parity = median(n, srtbuf) > 0 ? +1 : -1;

   for(i=0; i<n; i++) srtbuf[i] = threepar[4*i];
   par[0] = parity * median(n, srtbuf);

   for(i=0; i<n; i++) srtbuf[i] = threepar[1+4*i];
   par[1] = median(n, srtbuf);

   for(i=0; i<n; i++) srtbuf[i] = threepar[2+4*i];
   par[2] = median(n, srtbuf);

   for(i=0; i<n; i++) srtbuf[i] = threepar[3+4*i];
   par[3] = log(median(n, srtbuf));

   par[4] = cx;
   par[5] = cy;

   if(TEST > 1) {
      printf("Bootstrap params: %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %d\n", 
	     par[0], par[1], par[2], par[3], par[4], par[5], par[6], par[7], (int)par[8]);
   }

/* Evaluate this fit */
   if(TEST > 0) {
      printf("#   x       y       %cA     Dec     %cAcalc  Decalc dAngle   Azi     Alt  ompix    xcalc   ycalc    r      theta   <dx>  <dy> resid\n", 
	     parity>0?'H':'R', parity>0?'H':'R');
   }
   for(i=0; i<npt; i++) {
      fishsky(xyad[i].x, xyad[i].y, par, &ha, &dec, &azi, &alt, &jac);
      ang = srtbuf[i] = skyangle(ha, dec, parity*xyad[i].ha, xyad[i].dec);
      skyfish(parity*xyad[i].ha, xyad[i].dec, par, 
	      &xyad[i].xmap, &xyad[i].ymap, &xyad[i].thmap);
      if(TEST > 0) {
	 r = sqrt(pow(xyad[i].x-cx,2.0) + pow(xyad[i].y-cy,2.0));
	 u = xyad[i].ha/dr;
	 if(u > 180.0) u -= 360.0;
	 if(u < -180.0) u += 360.0;
	 printf("%7.1f %7.1f  %7.2f %7.2f  %7.2f %7.2f %5.2f %7.2f %7.2f %5.0f  %7.1f %7.1f %7.2f %8.5f %5.2f %5.2f %5.2f\n", 
		xyad[i].x, xyad[i].y, u, xyad[i].dec/dr,
		parity*ha/dr, dec/dr, ang/dr, azi/dr, alt/dr, 
		3600*3600*jac/(dr*dr),
		xyad[i].xmap, xyad[i].ymap, r, xyad[i].thmap,
		xyad[i].dxnear, xyad[i].dynear, xyad[i].resid);
      }
   }
   ang = median(npt, srtbuf);

   printf("%cA=%.3f DEC=%.3f AZ=%.3f SCALE=%.5f PARITY=%d MED=%.3f\n", 
	  parity>0?'H':'R',
	  parity*par[0]/dr, par[1]/dr, par[2]/dr, exp(par[3])/dr, 
	  (int)par[8], ang/dr);
}


/* Solve for pole, azimuth, scale, and parity from three points */
int threept(XYAD *xyad1, XYAD *xyad2, XYAD *xyad3, double cx, double cy, 
	    double *par, int *parity)
{
   double dr=atan(1.0)/45.0, pi=4*atan(1.0);
   double x0[3], x1[3], x2[3], x3[3], u1[2], u2[2], u3[2];
   double xu2[3], xu3[3], uu2[2], uu3[2];
   double r, dot, norm, u12, u13, nx0, nx1, x01, c1, uvw, xyzr;

/* Image points */
   u1[0] = xyad1->x - cx;
   u1[1] = xyad1->y - cy;
   u2[0] = xyad2->x - cx;
   u2[1] = xyad2->y - cy;
   u3[0] = xyad3->x - cx;
   u3[1] = xyad3->y - cy;

/* Sky points on sphere */
   x1[0] = cos(xyad1->dec) * cos(xyad1->ha);
   x1[1] = cos(xyad1->dec) * sin(xyad1->ha);
   x1[2] = sin(xyad1->dec);
	                  
   x2[0] = cos(xyad2->dec) * cos(xyad2->ha);
   x2[1] = cos(xyad2->dec) * sin(xyad2->ha);
   x2[2] = sin(xyad2->dec);

   x3[0] = cos(xyad3->dec) * cos(xyad3->ha);
   x3[1] = cos(xyad3->dec) * sin(xyad3->ha);
   x3[2] = sin(xyad3->dec);


/* Radius of sphere: make cam distance equal to sky distance */
   r = sqrt((pow(u1[0]-u2[0],2.0)+pow(u1[1]-u2[1],2.0)) / 
	    (pow(x1[0]-x2[0],2.0)+pow(x1[1]-x2[1],2.0)+pow(x1[2]-x2[2],2.0)));

   x1[0] *= r;   x1[1] *= r;   x1[2] *= r;
   x2[0] *= r;   x2[1] *= r;   x2[2] *= r;
   x3[0] *= r;   x3[1] *= r;   x3[2] *= r;

   if(TEST > 2) {
      printf("%7.1f %7.1f  %7.1f %7.1f %7.1f\n", u1[0],u1[1], x1[0],x1[1],x1[2]);
      printf("%7.1f %7.1f  %7.1f %7.1f %7.1f\n", u2[0],u2[1], x2[0],x2[1],x2[2]);
      printf("%7.1f %7.1f  %7.1f %7.1f %7.1f\n", u3[0],u3[1], x3[0],x3[1],x3[2]);
   }

/* Unit vectors in plane, 3D coords, G-S orthonormalized */
   xu2[0] = x2[0] - x1[0];
   xu2[1] = x2[1] - x1[1];
   xu2[2] = x2[2] - x1[2];
   norm = sqrt(xu2[0]*xu2[0]+xu2[1]*xu2[1]+xu2[2]*xu2[2]);
   xu2[0] /= norm;  xu2[1] /= norm;  xu2[2] /= norm; 

   xu3[0] = x3[0] - x1[0];
   xu3[1] = x3[1] - x1[1];
   xu3[2] = x3[2] - x1[2];
   dot = xu3[0]*xu2[0] + xu3[1]*xu2[1] + xu3[2]*xu2[2];
   xu3[0] -= dot*xu2[0];  xu3[1] -= dot*xu2[1];  xu3[2] -= dot*xu2[2];
   norm = sqrt(xu3[0]*xu3[0]+xu3[1]*xu3[1]+xu3[2]*xu3[2]);
   xu3[0] /= norm;  xu3[1] /= norm;  xu3[2] /= norm; 

/* Unit vectors in plane, 2D coords, G-S orthonormalized */
   uu2[0] = u2[0] - u1[0];
   uu2[1] = u2[1] - u1[1];
   norm = sqrt(uu2[0]*uu2[0]+uu2[1]*uu2[1]);
   uu2[0] /= norm;  uu2[1] /= norm;

   uu3[0] = u3[0] - u1[0];
   uu3[1] = u3[1] - u1[1];
   dot = uu3[0]*uu2[0] + uu3[1]*uu2[1];
   uu3[0] -= dot*uu2[0];  uu3[1] -= dot*uu2[1];
   norm = sqrt(uu3[0]*uu3[0]+uu3[1]*uu3[1]);
   uu3[0] /= norm;  uu3[1] /= norm;

   if(TEST > 2) {
      printf("%7.3f %7.3f  %7.3f %7.3f %7.3f\n", uu2[0],uu2[1], xu2[0],xu2[1],xu2[2]);
      printf("%7.3f %7.3f  %7.3f %7.3f %7.3f\n", uu3[0],uu3[1], xu3[0],xu3[1],xu3[2]);
   }

/* Plane origin, 3D coords */
   u12 = u1[0]*uu2[0] + u1[1]*uu2[1];
   u13 = u1[0]*uu3[0] + u1[1]*uu3[1];
   x0[0] = x1[0] - u12*xu2[0] - u13*xu3[0];
   x0[1] = x1[1] - u12*xu2[1] - u13*xu3[1];
   x0[2] = x1[2] - u12*xu2[2] - u13*xu3[2];
/* Push it out to the surface of the sphere */
   norm = sqrt(x0[0]*x0[0]+x0[1]*x0[1]+x0[2]*x0[2]);
   x0[0] *= r/norm;  x0[1] *= r/norm;  x0[2] *= r/norm; 

   if(TEST > 2) {
      printf("%7.3f %7.3f  %7.1f %7.1f %7.1f\n", u12, u13, x0[0], x0[1], x0[2]);
   }

/* Angles and azimuth */
   nx0 = acos(x0[2]/r);		// Angle from x0 to N
   nx1 = acos(x1[2]/r);		// Angle from x1 to N
   x01 = acos((x0[0]*x1[0]+x0[1]*x1[1]+x0[2]*x1[2])/(r*r)); // Angle x0 to x1
/* Law of cosines: angle between x1-x0 to N-x0 */
   c1 = acos((cos(nx1)-cos(nx0)*cos(x01))/(sin(nx0)*sin(x01)));

/* Parameters: alpha, delta, azimuth, log(scale) [rad/pix] */
   par[0] = atan2(x0[1], x0[0]);
   par[1] = asin(x0[2] / r);
#if 0
   cw = atan2(x1[1],x1[0]) - atan2(x0[1],x0[0]);
   if(cw > 0 || cw < -pi) {
      par[2] = pi - atan2(u1[1], u1[0]) - c1;
   } else {
      par[2] = pi - atan2(u1[1], u1[0]) + c1;
   }
#endif
   par[2] = pi - atan2(u1[1], u1[0]) + c1;
   if(par[2] > pi) par[2] -= 2*pi;
   if(par[2] < -pi) par[2] += 2*pi;
   par[3] = skyangle(par[0], par[1], xyad1->ha, xyad1->dec) /
		sqrt(u1[0]*u1[0]+u1[1]*u1[1]);

/* Parity */
/* w component of (u2-u1)x(u3-u1) */
   uvw = uu2[0]*uu3[1] - uu2[1]*uu3[0];
/* radial component of (x2-x1)x(x3-x1) */
   xyzr = x1[0]/r*(xu2[1]*xu3[2]-xu2[2]*xu3[1]) + 
          x1[1]/r*(xu2[2]*xu3[0]-xu2[0]*xu3[2]) + 
          x1[2]/r*(xu2[0]*xu3[1]-xu2[1]*xu3[0]);

   if(TEST > 1) {
      printf("%8.2f %8.2f %8.2f %8.4f %8.2f %8.4f %8.4f\n", 
	     par[0]/dr, par[1]/dr, par[2]/dr, par[3]/dr, uvw, xyzr, uvw*xyzr);
   }

   *parity = uvw*xyzr > 0 ? +1 : -1;

   return(0);
}

/* Find nearest neighbors to each point */
int nearest(int npt, XYAD *xyad, int nneigh)
{
   int j, k;
   double dist[MAXPT];
   XYAD *idx[MAXPT];
   for(k=0; k<npt; k++) {
      for(j=0; j<npt; j++) {
	 dist[j] = (xyad[k].x-xyad[j].x) * (xyad[k].x-xyad[j].x) +
 	           (xyad[k].y-xyad[j].y) * (xyad[k].y-xyad[j].y);
	 idx[j] = xyad + j;
      }
      heapsort(npt, dist, (void **)idx);
      xyad[k].near = (XYAD **)calloc(nneigh, sizeof(XYAD *));
      xyad[k].dist = (double *)calloc(nneigh, sizeof(double));
      for(j=0; j<nneigh; j++) {
	 xyad[k].near[j] = idx[j+1];
	 xyad[k].dist[j] = sqrt(dist[j+1]);
      }
   }

   return(0);
}

/* Evaluate residual of this point */
void resid(int nneigh, XYAD *xyad)
{
   int i;
   double u[MAXNEIGH];

/* Local drift in x */
   for(i=0; i<nneigh; i++) u[i] = xyad->near[i]->xmap - xyad->near[i]->x;
   heapsort(nneigh, u, NULL);
   xyad->dxnear = 0.5*(u[nneigh/2] + u[(nneigh-1)/2]);

/* Local drift in y */
   for(i=0; i<nneigh; i++) u[i] = xyad->near[i]->ymap - xyad->near[i]->y;
   heapsort(nneigh, u, NULL);
   xyad->dynear = 0.5*(u[nneigh/2] + u[(nneigh-1)/2]);

/* Offset of this point wrt local drift */
   xyad->resid = sqrt(
      (xyad->xmap-xyad->x - xyad->dxnear)*(xyad->xmap-xyad->x - xyad->dxnear) +
      (xyad->ymap-xyad->y - xyad->dynear)*(xyad->ymap-xyad->y - xyad->dynear));
}

/* Angle between HA1, Dec1 and HA2, Dec2 */
double skyangle(double ha1, double dec1, double ha2, double dec2)
{
   return(acos(cos(dec1)*cos(ha1) * cos(dec2)*cos(ha2) + 
	       cos(dec1)*sin(ha1) * cos(dec2)*sin(ha2) + 
	       sin(dec1) * sin(dec2) ) );
}

/* Convert HA,Dec to az, alt, rotation */
void HD_to_AZ(double ha, double dec, double lat, 
	      double *alt, double *az, double *rot)
{
   double sind, sinh, cosh;

   sind = sin(dec) * sin(lat) + cos(dec) * cos(ha) * cos(lat);
   *alt  = asin(sind);

   sinh = - cos(dec) * sin(ha);
   cosh =   sin(dec) * cos(lat) - cos(dec) * cos(ha) * sin(lat);

   *az = atan2(sinh, cosh);

   sinh = -cos(*az) * sin(*alt) * sin(ha) * sin(lat) + sin(*az) *
      sin(*alt) * cos(ha) - sin(ha) * cos(*alt) * cos(lat);
   cosh = -sin(*az) * sin(ha) * sin(lat) - cos(*az) * cos(ha);
   *rot = -atan2(sinh, cosh);
}


double median(int n, double *data)
{
   heapsort(n, data, NULL);
   return(0.5*(data[(n-1)/2]+data[n/2]));
}
 
/* heapsort */
void heapsort(int n, double *data, void **idx)
{
   int i;
   double temp;
   void *itmp;
   for(i=n/2; i>=0; i--) siftdown(i, n-1, data, idx);
 
   for(i=n-1; i>=1; i--) {  // Swap
      temp = data[0];
      data[0] = data[i];
      data[i] = temp;
      if(idx != NULL) {
	 itmp = idx[0];
	 idx[0] = idx[i];
	 idx[i] = itmp;
      }
      siftdown(0, i-1, data, idx);
   }
}
 
void siftdown(int root, int bottom, double *data, void **idx)
{
  int maxchild = 2*root + 1;
  double temp;
  void *itmp;
  if(maxchild < bottom) {		// Find the biggest child
    int otherchild = maxchild + 1;	// Reversed for stability
    maxchild = (data[otherchild] > data[maxchild]) ? otherchild : maxchild;
  } else {				// Don't overflow
    if(maxchild > bottom) return;
  }
/* If we have the correct ordering, we are done. */
  if(data[root] >= data[maxchild]) return;
 
/* Swap */
  temp = data[root];
  data[root] = data[maxchild];
  data[maxchild] = temp;
  if(idx != NULL) {
     itmp = idx[root];
     idx[root] = idx[maxchild];
     idx[maxchild] = itmp;
  }
  siftdown(maxchild, bottom, data, idx);	// Tail recursion
}


/* Show us the covariance matrix */
int printmat(int n, int *usepar, double *a, double *y, double chi)
{
   int i, j, iu, ju, nuse=0;
   for(i=0; i<n; i++) nuse += usepar[i];
   for(i=iu=0; i<n; i++) {
      for(j=ju=0; j<n; j++) {
	 if(j < i) {
	    printf("             ");
	 } else if(usepar[i] == 0 || usepar[j] == 0) {
	    printf(" %12.4g", 0.0);
	 } else {
	    printf(" %12.4g", iu==ju?a[ju+iu*nuse]*sqrt(chi):a[ju+iu*nuse]);
	 }
	 ju += usepar[j];
      }
      printf("   | %12.4g\n", y[i]);
      iu += usepar[i];
   }
   return(0);
}

void syntax(char *prog)
{
   printf("Syntax: %s -in F [options]\n", prog);
   printf("Options include:\n");
   printf("    -in F          Input file of x y HA Dec (or HA Dec for -sky2xy)\n");
   printf("    -lat L         Observatory latitude\n");
   printf("    -cx X          Image center in x [pix]\n");
   printf("    -cy Y          Image center in y [pix]\n");
   printf("    -ha R          Fisheye pole HA (sets parity +1) [deg]\n");
   printf("    -ra R          Fisheye pole RA (sets parity -1) [deg]\n");
   printf("    -dec D         Fisheye pole Dec [deg]\n");
   printf("    -az Z          Fisheye azimuth [deg]\n");
   printf("    -scale S       Fisheye scale [pix/deg]\n");
   printf("    -quad Q        Fisheye quadratic term [rad^-1]");
   printf("    -cube C        Fisheye cubic term [rad^-2]\n");
   printf("    -xy2sky        Evaluate only x y to x y HA Dec\n");
   printf("    -sky2xy        Evaluate only HA Dec to x y HA Dec\n");
   printf("    -fixctr        Fix center cx,cy; do not fit\n");
   printf("    -fixdist       Fix distortion quad,cube; do not fit\n");
   printf("    -maxiter N     Max number of iterations (default 20)\n");
   printf("    -out F         Output file\n");
   printf("    -verb          Moderate verbosity\n");
   printf("    -VERB          High verbosity\n");
   printf(" deprecated:\n");
   printf("    -parity s      Image parity: HA(-1) vs HA(+1)\n");
}
