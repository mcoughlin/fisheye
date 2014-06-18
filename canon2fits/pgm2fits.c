/* Convert Canon CR2 or Nikon NEF PGM output from mdcraw to FITS */
/* 131204 Change name of Canon 5d MarkII and MarkIII to have a space! */
/* v.2 130908 Add Nikon */
/* 100614 - John Tonry */

/* Bundle up the R/G1/G2/B into a 3D FITS file with four or three planes */
/*
 * Syntax: pgm2fits [options] < IMG_0001.pgm > IMG_0001.fits
 *                  
 * e.g.
 * mdcraw -i -v IMG_0031.CR2 > metadata.031
 * mdcraw -4 -j -c -D -t 0 IMG_0001.CR2 | 
 *              pgm2fits -rgb -hdr metadata.031 > IMG_0001.fits
 *
 * where options include:
 *
 *  -rgb	Average G1 and G2, and write R, G, B instead of R, G1, G2, B
 *  -rggb       Write R, G1, G2, B (4 plane)
 *  -bw         Average R+G1+G2+B, write monochrome (1 plane)
 *  -r          Write R (1 plane)
 *  -g          Write 0.5*(G1+G2) (1 plane)
 *  -b          Write B (1 plane)
 *  -g1         Write G1 (1 plane)
 *  -g2         Write G2 (1 plane)
 *  -div2	Divide the result by 2
 *  -bayer	Write RGGB pixels as full-sized image
 *  -hdr        Read metadata file from mdcraw; add to header
 *  -notrim     Don't trim to thumb size from metadata file
 *  -image=x0,y0:x1,y1    Designation of what's image
 *  -bias=x0,y0:x1,y1     Designation of what's bias
 *  -zone=H     Correct reported time to UTC (H=reported-UTC e.g. H=-10 for HST)
 *
 * Note that this require (m)dcraw, the modified version of Dave Coffin's
 * library to convert camera raw files to something useful
 */
#define _XOPEN_SOURCE
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>

#define NFITS 2880	/* FITS magic size */
#define BORDER 20	/* Bias border (raw pixels) */

static int _ENDIAN_TEST = 1;
#define LOWENDIAN() (*(char*)&_ENDIAN_TEST)

#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define ABS(a) (((a) > 0) ? (a) : -(a))

typedef struct imgfmt {
   char *camname;	/* Camera name */
   int ix0;		/* Image start (Bayer) pixel in x */
   int inx;		/* Number of image (Bayer) pixels in x */
   int iy0;		/* Image start pixel in y */
   int iny;		/* Number of image pixels in y */
   int bx0;		/* Bias start pixel in x */
   int bnx;		/* Number of bias pixels in x */
   int by0;		/* Bias start pixel in y */
   int bny;		/* Number of bias pixels in y */
   int bias;		/* Nominal bias level */
} IMGFMT;

/* Various types of cameras we know about */
IMGFMT knowncam[]={
//  "Name reported by mdcraw" x0   nx  y0   ny    bx0 bnx by0  bny   bias
   {"Canon EOS 5D Mark II",   79, 2817, 0, 1876,   41, 35, 10, 1847, 1024},
   {"Canon EOS 5D Mark III",  61, 2897, 0, 1935,    2, 55, 10, 1905, 2048},
   {"Canon EOS 6D",           40, 2743, 0, 1828,    2, 16, 10, 1808, 2048},
   {"NIKON D700",              1, 2142, 0, 1422,    0,  1, 10, 1402, 0},
   {"NIKON D800",              0, 3688, 0, 2462, 3689,  6, 10, 2442, 0}
};
static int nknown=sizeof(knowncam)/sizeof(IMGFMT);

static int TEST=0;		/* Define for debug output */

int nominalhead(int nx, int ny, int nz, int sx, int sy, int zone, int fd, char *hdr,
                int *nhead, char *header, int *imx, int *imy, IMGFMT *imfmt);
int metahead(char *fname, int zone, int *nhead, char *header, int *imx, int *imy, IMGFMT *imfmt);
int dmodhead(int *nhead, char *header, char *key, double val, char *cmt);
int imodhead(int *nhead, char *header, char *key, int val, char *cmt);
double jdate(int year, int month, int day, double ut);
void syntax(char *prog);

int main(int argc, char *argv[])
{
   int i, j, k, qp, qm, med, nread, nbuf, nhead, nbias;
   int div2, trim, zone;
   char line[NFITS], header[NFITS];
   IMGFMT imfmt, hdrfmt;
   unsigned short *rawdata, *data;
   unsigned short *clrptr[4], *rpt, *g1pt, *g2pt, *bpt, *dpt;
   unsigned short buf[NFITS/sizeof(unsigned short)];
   int count[65536];
   int NX, NY;			// Number of raw pixels in image
   int thx, thy;		// Thumb size in Bayer pixels
   int imx, imy, sx, sy;	// Desired image area in Bayer pixels
   int nclr, maxval, clr;
   int fd=STDOUT_FILENO;
   double bias, noise;
   char *hdr;

/* Defaults */
   nclr = 3;		/* Write RGB by default */
   clr = 0;		/* Which of various 1 colors? */
   div2 = 0;		/* Divide by 2? */
   hdr = NULL;		/* Info file from dcraw with metadata */
   trim = 1;		/* Trim to thumb size? */
   imx = imy = sx = sy = 0;	/* Requests for specific areas */
   imfmt.inx = imfmt.bnx = imfmt.bias = 0;	/* Requests for specific areas */
   zone = 0;		/* Error in camera date wrt UTC */

/* Parse the arguments */
   for(i=1; i<argc; i++) {

      if(strcmp(argv[i], "-rgb") == 0) {
	 nclr = 3;

      } else if(strcmp(argv[i], "-rggb") == 0) {
	 nclr = 4;

      } else if(strcmp(argv[i], "-bw") == 0) {
	 nclr = 1;
	 clr = 0;

      } else if(strcmp(argv[i], "-r") == 0) {
	 nclr = 1;
	 clr = 1;

      } else if(strcmp(argv[i], "-g") == 0) {
	 nclr = 1;
	 clr = 2;

      } else if(strcmp(argv[i], "-b") == 0) {
	 nclr = 1;
	 clr = 3;

      } else if(strcmp(argv[i], "-g1") == 0) {
	 nclr = 1;
	 clr = 21;

      } else if(strcmp(argv[i], "-g2") == 0) {
	 nclr = 1;
	 clr = 22;

      } else if(strcmp(argv[i], "-bayer") == 0) {
	 nclr = 1;
	 clr = -1;

      } else if(strcmp(argv[i], "-div2") == 0) {
	 div2 = 1;

      } else if(strcmp(argv[i], "-notrim") == 0) {
	 trim = 0;

      } else if(strncmp(argv[i], "-bias=", 6) == 0) {
	 if(sscanf(argv[i]+6,"%d,%d:%d,%d", &sx, &sy, &imx, &imy) != 4) {
	    fprintf(stderr, "Cannot scan 'x0,y0:x1,y1' from %s\n", argv[i]);
	    exit(1);
	 }
	 imfmt.bx0 = sx;
	 imfmt.by0 = sy;
	 imfmt.bnx = imx - sx + 1;
	 imfmt.bny = imy - sy + 1;

      } else if(strncmp(argv[i], "-image=", 7) == 0) {
	 if(sscanf(argv[i]+7,"%d,%d:%d,%d", &sx, &sy, &imx, &imy) != 4) {
	    fprintf(stderr, "Cannot scan 'x0,y0:x1,y1' from %s\n", argv[i]);
	    exit(1);
	 }
	 imfmt.bx0 = sx;
	 imfmt.by0 = sy;
	 imfmt.bnx = imx - sx + 1;
	 imfmt.bny = imy - sy + 1;

      } else if(strncmp(argv[i], "-zone=", 6) == 0) {
	 sscanf(argv[i]+6, "%d", &zone);

      } else if(strcmp(argv[i], "-hdr") == 0) {
	 hdr = argv[++i];

      } else if(strcmp(argv[i], "-test") == 0) {
	 TEST = 1;

      } else {
	 fprintf(stderr, "Unrecognized argument `%s'\n", argv[i]);
	 syntax(argv[0]);
	 exit(1);
      }   
   }

/* Read data header and perform some sanity checks */
   if( (i=scanf("%s %d %d %d", line, &NX, &NY, &maxval)) != 4) {
      fprintf(stderr, "Cannot scan header of input\n");
      exit(1);
   }
   if(strcmp(line, "P5") != 0) {
      fprintf(stderr, "Magic number is %s, not P5 required\n", line);
      exit(1);
   }
   if(maxval != 65535) {
      fprintf(stderr, "maxval = %d, does not indicate 16 bit FITS\n", maxval);
      exit(1);
   }

/* Skip white space -- mdcraw writes newlines */
   nread = fread(line, 1, 1, stdin);
   if(line[0] != '\n') {
      fprintf(stderr, "Not synced!  Did not get two newlines after maxval\n");
      exit(1);
   }

   if(TEST) {
      fprintf(stderr, "RawSize= %d x %d", NX, NY);
   }

/* We're now positioned at the data */

/* We know how big the image is; allocate space for the whole thing */
   rawdata = (unsigned short *)calloc(NX*NY, sizeof(unsigned short));
/* Addresses of each Bayer colors: R-G1 / G2-B, top row first */
   clrptr[0] = rawdata+NX;			/* R */
   clrptr[1] = rawdata+NX+1;			/* G1 */
   clrptr[2] = rawdata;				/* G2 */
   clrptr[3] = rawdata+1;			/* B */

/* Allocate space for each extracted plane from Bayer pixels */
   if(clr >= 0) {
      data = (unsigned short *)calloc(NX/2*NY/2, sizeof(unsigned short));
   } else {
      data = (unsigned short *)calloc(NX*NY, sizeof(unsigned short));
   }

/* Read all the data */
   nread = fread(rawdata, sizeof(short int), NX*NY, stdin);
   if(nread < NX*NY) {
      fprintf(stderr, "Short read, needed %d, got %d, %d short\n",
              NX*NY, nread, NX*NY-nread);
      fprintf(stderr, "continuing anyway...\n");
   }

/* PGM and FITS are big-endian, swab for arithmetic */
   if(LOWENDIAN()) swab(rawdata, rawdata, NX*NY*sizeof(unsigned short));

/* Create a nominal FITS header; read any metadata */
   nominalhead(NX/2, NY/2, nclr, 0, 0, zone, fd, hdr, &nhead, header, &thx, &thy, &hdrfmt);

/* Set imfmt to the desired size */
   if(!trim) {						// Priority 1: no trim!
      imfmt.ix0 = 0;
      imfmt.inx = NX/2;
      imfmt.iy0 = 0;
      imfmt.iny = NY/2;
      imfmt.bx0 = imfmt.bnx = imfmt.by0 = imfmt.bny = 0;
      imfmt.bias = 0;

   } else if(imfmt.inx > 0 && imfmt.bnx > 0) {		// Priority 2: command line
      if(hdrfmt.inx > 0) imfmt.bias = hdrfmt.bias;

   } else if(hdrfmt.inx > 0 && hdrfmt.bnx > 0) {	// Priority 3: known camera
      imfmt = hdrfmt;

   } else if(hdr != NULL) {				// Priority 4: thumb from header
      imfmt.ix0 = (NX - thx) / 2;	// Assume image right
      imfmt.inx = thx / 2;
      imfmt.iy0 = 0;			// Assume image bottom
      imfmt.iny = thy / 2;
      imfmt.bx0 = 0;			// Assume bias left
      imfmt.bnx = (NX-thx) / 2;
      imfmt.by0 = 0;			// Assume bias bottom
      imfmt.bny = thy / 2;
      imfmt.bias = 0;

   } else {						// Priority 5: treat like no trim
      imfmt.ix0 = 0;
      imfmt.inx = NX/2;
      imfmt.iy0 = 0;
      imfmt.iny = NY/2;
      imfmt.bx0 = imfmt.bnx = imfmt.by0 = imfmt.bny = 0;
      imfmt.bias = 0;
   }

/* Update header with size of image we're going to write */
   imodhead(&nhead, header, "NAXIS1  ", imfmt.inx, "Image size in x");
   imodhead(&nhead, header, "NAXIS2  ", imfmt.iny, "Image size in y");
   imodhead(&nhead, header, "CNPIX1  ", 0,  "Image offset in x");
   imodhead(&nhead, header, "CNPIX2  ", 0,  "Image offset in y");

/* Add the image and bias area parameters to the header */
   imodhead(&nhead, header, "IMSPIX1 ", imfmt.ix0,  "Image start in x");
//   imodhead(&nhead, header, "IMNPIX1 ", imfmt.inx,  "Image npix in x");
   imodhead(&nhead, header, "IMSPIX2 ", imfmt.iy0,  "Image start in y");
//   imodhead(&nhead, header, "IMNPIX2 ", imfmt.iny,  "Image npix in y");

   imodhead(&nhead, header, "BSPIX1  ", imfmt.bx0,  "Bias start in x");
   imodhead(&nhead, header, "BNPIX1  ", imfmt.bnx,  "Bias npix in x");
   imodhead(&nhead, header, "BSPIX2  ", imfmt.by0,  "Bias start in y");
   imodhead(&nhead, header, "BNPIX2  ", imfmt.bny,  "Bias npix in y");

/* Actually dumping Bayer pixels as a 2x2 array */
   if(clr == -1) {
      imodhead(&nhead, header, "NAXIS1  ", 2*imfmt.inx, "Image size in x");
      imodhead(&nhead, header, "NAXIS2  ", 2*imfmt.iny, "Image size in y");
   }

/* Iterate over all the planes requested */
   nbuf = 0;
   for(k=0; k<nclr; k++) {
/* Collect all the desired pixels from the raw ones */
      for(j=0; j<NY; j+=2) {
	 rpt =  clrptr[0] + (NY-2-j)*NX;
	 g1pt = clrptr[1] + (NY-2-j)*NX;
	 g2pt = clrptr[2] + (NY-2-j)*NX;
	 bpt =  clrptr[3] + (NY-2-j)*NX;
	 dpt =  data + (j/2)*(NX/2);
	 if(clr == -1) {		// all Bayer pixels
	    for(i=0; i<NX; i+=2) {
	       data[i+0+(j+0)*NX] = rpt[i];
	       data[i+1+(j+0)*NX] = g1pt[i];
	       data[i+0+(j+1)*NX] = g2pt[i];
	       data[i+1+(j+1)*NX] = bpt[i];
	    }
	    if(div2) {
	       for(i=0; i<NX; i++) data[i+j*NX] /= 2;
	       for(i=0; i<NX; i++) data[i+(j+1)*NX] /= 2;
	    }
	 } else if(nclr == 1) {
	    if(clr == 0) {	// BW=(R+G1+G2+B)/4
	       for(i=0; i<NX/2; i++) *dpt++ = 
					((int)rpt[2*i] + (int)g1pt[2*i] +
					 (int)g2pt[2*i] + (int)bpt[2*i]) / 4;
	    } else if(clr == 1) {	// R
	       for(i=0; i<NX/2; i++) *dpt++ = rpt[2*i];
	    } else if(clr == 2) {	// G=(G1+G2)/2
	       for(i=0; i<NX/2; i++) *dpt++ = ((int)g1pt[2*i] + (int)g2pt[2*i]) / 2;
	    } else if(clr == 21) {	// G1
	       for(i=0; i<NX/2; i++) *dpt++ = g1pt[2*i];
	    } else if(clr == 22) {	// G2
	       for(i=0; i<NX/2; i++) *dpt++ = g2pt[2*i];
	    } else if(clr == 3) {	// B
	       for(i=0; i<NX/2; i++) *dpt++ = bpt[2*i];
	    }

	    if(div2) {
	       dpt = data + (j/2)*(NX/2);
	       for(i=0; i<NX/2; i++) *dpt++ /= 2;
	    }

	 } else if(nclr == 3) {		// R, (G1+G2)/2, B
	    if(k != 1) {
	       rpt =  clrptr[k] + (NY-2-j)*NX;
	       for(i=0; i<NX/2; i++) *dpt++ = rpt[2*i];
	    } else {
	       g1pt = clrptr[1] + (NY-2-j)*NX;
	       g2pt = clrptr[2] + (NY-2-j)*NX;
	       for(i=0; i<NX/2; i++) *dpt++ = ((int)g1pt[2*i] + (int)g2pt[2*i]) / 2;
	    }

	    if(div2) {
	       dpt =  data + (j/2)*(NX/2);
	       for(i=0; i<NX/2; i++) *dpt++ /= 2;
	    }

	 } else {		// R, G1, G2, B
	    rpt =  clrptr[k] + (NY-2-j)*NX;
	    for(i=0; i<NX/2; i++) *dpt++ = rpt[2*i];

	    if(div2) {
	       dpt =  data + (j/2)*(NX/2);
	       for(i=0; i<NX/2; i++) *dpt++ /= 2;
	    }

	 }
      }

/* Calculate bias and noise between +/-3 * quartile range to avoid hot pix */
      if(k == 0) {
	 nbias = 0;
/* Zero bias count buffer */
	 for(i=0; i<65536; i++) count[i] = 0;
	 for(j=imfmt.by0; j<imfmt.by0+imfmt.bny; j++) {
	    for(i=imfmt.bx0; i<imfmt.bx0+imfmt.bnx; i++) {
	       count[data[i+j*(NX/2)]] += 1;
	       nbias++;
	    }
	 }
	 for(i=j=qm=med=qp=0; i<65536; i++) {
	    j += count[i];
	    if(qm == 0 && j >= nbias/4) qm = i;
	    if(med == 0 && j >= nbias/2) med = i;
	    if(j >= (3*nbias)/4) {
	       qp = i;
	       break;
	    }
	 }
	 bias = noise = 0.0;
	 for(i=MAX(0, med-3*(qp-qm)), j=0; i<MIN(65536, med+3*(qp-qm)); i++) {
	    bias += i * (double)count[i];
	    noise += i*i * (double)count[i];
	    j += count[i];
	 }
	 bias /= MAX(1,j);
	 noise = noise/MAX(1,j) - bias*bias;
	 if(noise > 0) noise = sqrt(noise);

	 if(TEST) {
//	    fprintf(stderr, " median %d %d %d %d %d %d \n", nbias, qm, med, qp,
//		 med-3*(qp-qm), med+3*(qp-qm));
	    fprintf(stderr, " Bias= %8.2f Noise= %8.2f\n", bias, noise);
	 }

/* Add to header */
/* Canon-style leave the bias */
	 if(hdrfmt.inx == 0 || imfmt.bias > 0) {
	    dmodhead(&nhead, header, "BIAS    ", bias, "Bias level");
	 } else {
/* Nikon-style subtract the bias */
	    imodhead(&nhead, header, "BIAS    ", 0, "Residual bias level");
	    dmodhead(&nhead, header, "BIASBLK ", bias, "Bias level in black");
	 }
	 dmodhead(&nhead, header, "NOISE   ", noise, "Noise level");

/* Write the FITS header */
	 i = write(fd, header, NFITS);
      }

/* Write the desired pixels from the processed pixel array */
      if(clr >= 0) {
	 for(j=imfmt.iy0; j<imfmt.iy0+imfmt.iny; j++) {
	    for(i=imfmt.ix0; i<imfmt.ix0+imfmt.inx; i++) {
	       buf[nbuf++] = data[i+j*(NX/2)];
	       if(nbuf == NFITS/sizeof(unsigned short) ) {
		  if(LOWENDIAN()) swab(buf, buf, NFITS);
		  nread = fwrite(buf, sizeof(unsigned short), 
				 NFITS/sizeof(unsigned short), stdout);
		  nbuf = 0;
	       }
	    }
	 }
      } else {
	 for(j=2*imfmt.iy0; j<2*imfmt.iy0+2*imfmt.iny; j++) {
	    for(i=2*imfmt.ix0; i<2*imfmt.ix0+2*imfmt.inx; i++) {
	       buf[nbuf++] = data[i+j*NX];
	       if(nbuf == NFITS/sizeof(unsigned short) ) {
		  if(LOWENDIAN()) swab(buf, buf, NFITS);
		  nread = fwrite(buf, sizeof(unsigned short), 
				 NFITS/sizeof(unsigned short), stdout);
		  nbuf = 0;
	       }
	    }
	 }
      }
   }

/* Pad out the last bit and write it */
   if(nbuf != 0) {
      bzero(buf+nbuf, NFITS-nbuf*sizeof(unsigned short));
      if(LOWENDIAN()) swab(buf, buf, NFITS);
      nread = fwrite(buf, sizeof(unsigned short), 
                     NFITS/sizeof(unsigned short), stdout);
   }

   exit(0);
}

/* Create a nominal header */
int nominalhead(int nx, int ny, int nz, int sx, int sy, int zone, int fd, char *meta,
                int *nhead, char *header, int *imx, int *imy, IMGFMT *imfmt)
{
   int i, nh;
   for(i=0; i<NFITS; i++) header[i] = ' ';
   nh = 0;
   sprintf(header+80*nh++, "%s", "SIMPLE  =                    T");
   sprintf(header+80*nh++, "%s%20d",   "BITPIX  = ", 16);
   sprintf(header+80*nh++, "%s%20d",   "NAXIS   = ", (nz==1) ? 2 : 3);
   sprintf(header+80*nh++, "%s%20d",   "NAXIS1  = ",nx);
   sprintf(header+80*nh++, "%s%20d",   "NAXIS2  = ",ny);
   sprintf(header+80*nh++, "%s%20d",   "NAXIS3  = ",nz);
   sprintf(header+80*nh++, "%s%20d",   "CNPIX1  = ",sx);
   sprintf(header+80*nh++, "%s%20d",   "CNPIX2  = ",sy);
   sprintf(header+80*nh++, "%s%20.3f", "BSCALE  = ",1.0);
   sprintf(header+80*nh++, "%s%20.3f", "BZERO   = ",0.0);
   sprintf(header+80*nh++, "%s",       "END     ");
   for(i=0; i<NFITS; i++) if(header[i] == '\0') header[i] = ' ';
   *nhead = nh;

   *imx = *imy = 0;
   imfmt->inx = 0;
/* Parse the mdcraw metadata file and create FITS cards */
   if(meta != NULL) metahead(meta, zone, nhead, header, imx, imy, imfmt);

   return(0);
}

/* Add a double entry to header */
int dmodhead(int *nhead, char *header, char *key, double val, char *cmt)
{
   int i;
   for(i=0; i<*nhead-1; i++) if(strncmp(header+80*i, key, 8) == 0) break;
   sprintf(header+80*i, "%-8s= %20.3f // %s", key, val, cmt);
   if(i == *nhead-1) {
      sprintf(header+80*(i+1), "%s",       "END     ");
      *nhead += 1;
   }
   for(i=0; i<NFITS; i++) if(header[i] == '\0') header[i] = ' ';
   return(0);
}

/* Add an integer entry to header */
int imodhead(int *nhead, char *header, char *key, int val, char *cmt)
{
   int i;
   for(i=0; i<*nhead-1; i++) if(strncmp(header+80*i, key, 8) == 0) break;
   sprintf(header+80*i, "%-8s= %20d // %s", key, val, cmt);
   if(i == *nhead-1) {
      sprintf(header+80*(i+1), "%s",       "END     ");
      *nhead += 1;
   }
   for(i=0; i<NFITS; i++) if(header[i] == '\0') header[i] = ' ';
   return(0);
}

/* Add metadata to the FITS header */
int metahead(char *fname, int zone, int *nhead, char *header, int *imx, int *imy, IMGFMT *imfmt)
{
   int i, nh;
   FILE *fp;
   char line[1024], filename[1024], timestamp[80], camera[80];
   char daymult[80], cammult[80];
   double iso, fno, etime, persec, flen, mjd;
   double cam[4];
   int day, month, year, hour, min, sec;
   char mon[20];
   char *moname[]={"Jan", "Feb", "Mar", "Apr", "May", "Jun", 
		    "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"};

// FILENAME=Filename: /local/jt/CanonUnix/5D2/IMG_0031.CR2
// TIMESTMP=Timestamp: Mon Jun 14 05:45:07 2010
// DATEOBS = '2010-06-14T05:45:07'
// MJD-OBS =  55361.23966
// INSTRUME=Camera: Canon EOS 5D Mark II
// ISO     =ISO speed: 200
// EXPTIME =Shutter: 1/24.7 sec
// APERTURE=Aperture: f/1.4
// FOCLEN  =Focal length: 24.0 mm
// DAYMULT =Embedded ICC profile: no
// CAMMULT =Number of raw images: 1
// THUMB1  =Thumb size:  5616 x 3744
// THUMB2  =
//          Full size:   1448 x 3804
//          Image size:  1448 x 3804
//          Output size: 1448 x 3804
//          Raw colors: 3
// DAYMULT =Daylight multipliers: 2.224558 0.928662 1.164364
// CAMMULT =Camera multipliers: 2093.000000 1024.000000 1736.000000 1024.000000

   for(nh=0; nh<NFITS/80; nh++) if(strncmp(header+80*nh, "END     ", 8) == 0 || 
                            strncmp(header+80*nh, "        ", 8) == 0) break;
   if(nh == NFITS/80) nh = 0;

   if( (fp=fopen(fname, "r")) == NULL) {
      fprintf(stderr, "Cannot read metadata from `%s'\n", fname);
      return(-1);
   }

   while(1) {
      if(fgets(line, 1024, fp) == NULL) break;
      if(strncmp(line, "Filename:", 9) == 0) {
	 sscanf(line, "Filename: %s", filename);
	 sprintf(header+80*nh++, "%-8s= '%s' // %s", 
	         "FILENAME", filename, "Original file name");

      } else if(strncmp(line, "Timestamp:", 10) == 0) {
	 sscanf(line, "Timestamp: %[^\n]", timestamp);
	 sprintf(header+80*nh++, "%-8s= '%s' // %s", 
	         "TIMESTMP", timestamp, "Camera timestamp");

	 sscanf(timestamp, "%*s %s %d %d:%d:%d %d", 
	        mon, &day, &hour, &min, &sec, &year);
	 for(month=0; month<12; month++) 
	    if(strncasecmp(mon, moname[month], 3) == 0) break;
	 sprintf(header+80*nh++, "%-8s= '%04d-%02d-%02dT%02d:%02d:%02d' // %s", 
	 "DATE-OBS", year, month+1, day, hour, min, sec, "Exposure date");

	 mjd = jdate(year, month+1, day, (double)(hour-zone)+min/60.0+sec/3600.0)
	    - 2400000.5;              /* Convert from JD to MJD */
	 sprintf(header+80*nh++, "%-8s= %20.6f // %s", 
	         "MJD-OBS", mjd, "Exposure modified Julian date");

	 sprintf(header+80*nh++, "%-8s= %20d // %s", 
	         "ZONEHOUR", zone, "[hr] Date minus UTC");

      } else if(strncmp(line, "Camera:", 7) == 0) {
	 sscanf(line, "Camera: %[^\n]", camera);
	 sprintf(header+80*nh++, "%-8s= '%s' // %s", 
	         "INSTRUME", camera, "Camera name");

/* Add some bias and image info if we recognize the camera */
	 for(i=0; i<nknown; i++) {
	    if(strcmp(camera, knowncam[i].camname) == 0) {
	       *imfmt = knowncam[i];
//	       sprintf(header+80*nh++, "%-8s= %20d // %s", "IMX0", knowncam[i].ix0, "Image start in x");
//	       sprintf(header+80*nh++, "%-8s= %20d // %s", "IMNX", knowncam[i].inx, "Image npix in x");
//	       sprintf(header+80*nh++, "%-8s= %20d // %s", "IMY0", knowncam[i].iy0, "Image start in y");
//	       sprintf(header+80*nh++, "%-8s= %20d // %s", "IMNY", knowncam[i].iny, "Image npix in y");
//	       sprintf(header+80*nh++, "%-8s= %20d // %s", "BIASX0", knowncam[i].bx0, "Bias start in x");                              
//	       sprintf(header+80*nh++, "%-8s= %20d // %s", "BIASNX", knowncam[i].bnx, "Bias npix in x");
//	       sprintf(header+80*nh++, "%-8s= %20d // %s", "BIASY0", knowncam[i].by0, "Bias start in y");
//	       sprintf(header+80*nh++, "%-8s= %20d // %s", "BIASNY", knowncam[i].bny, "Bias npix in y");
	       sprintf(header+80*nh++, "%-8s= %20d // %s", "BIASNOM", knowncam[i].bias, "Nominal bias level");
	       break;
	    }
	 }

      } else if(strncmp(line, "ISO speed:", 10) == 0) {
	 sscanf(line, "ISO speed: %lf", &iso);
	 sprintf(header+80*nh++, "%-8s= %20.1f // %s", 
	         "ISO", iso, "ISO speed");
 
      } else if(strncmp(line, "Shutter: 1/", 11) == 0) {
	 sscanf(line, "Shutter: 1/%lf", &persec);
	 sprintf(header+80*nh++, "%-8s= %20.6f // %s", 
	         "EXPTIME", 1/persec, "[sec] Exposure time");

      } else if(strncmp(line, "Shutter:", 8) == 0) {
	 sscanf(line, "Shutter: %lf", &etime);
	 sprintf(header+80*nh++, "%-8s= %20.6f // %s", 
	         "EXPTIME", etime, "[sec] Exposure time");

      } else if(strncmp(line, "Aperture: f/", 12) == 0) {
	 sscanf(line, "Aperture: f/%lf", &fno);
	 sprintf(header+80*nh++, "%-8s= %20.1f // %s", 
	         "APERTURE", fno, "Aperture f-number");

      } else if(strncmp(line, "Focal length:", 13) == 0) {
	 sscanf(line, "Focal length: %lf", &flen);
	 sprintf(header+80*nh++, "%-8s= %20.0f // %s", 
	         "FLEN", flen, "[mm] Lens focal length");

      } else if(strncmp(line, "Thumb size:", 11) == 0) {
	 sscanf(line, "Thumb size: %d x %d", imx, imy);
	 sprintf(header+80*nh++, "%-8s= %20d // %s", 
	         "THUMB1", *imx, "Actual image size (x)");
	 sprintf(header+80*nh++, "%-8s= %20d // %s", 
	         "THUMB2", *imy, "Actual image size (y)");

      } else if(strncmp(line, "Daylight multipliers:", 21) == 0) {
	 sscanf(line, "Daylight multipliers: %[^\n]", daymult);
	 sprintf(header+80*nh++, "%-8s= '%s' // %s", 
	         "DAYMULT", daymult, "Daylight multiplier factors");

      } else if(strncmp(line, "Camera multipliers:", 19) == 0) {
	 sscanf(line, "Camera multipliers: %[^\n]", cammult);
	 sscanf(cammult, "%lf %lf %lf %lf", &cam[0], &cam[1], &cam[2], &cam[3]);
	 sprintf(header+80*nh++, "%-8s= '%.1f %.1f %.1f %.1f' // %s", 
	 "CAMMULT", cam[0], cam[1], cam[2], cam[3], "Camera multiplier factors");

      } else {
      }
   }
   sprintf(header+80*nh++, "%s",       "END     ");
   for(i=0; i<NFITS; i++) if(header[i] == '\0') header[i] = ' ';
   *nhead = nh;
   fclose(fp);
   return(0);
}

double jdate(int year, int month, int day, double ut)
{
   double tmp;
   tmp = 367.0*year + 0.5 + ut/24.0;
   tmp = tmp - ((7*(year+((month+9)/12)))/4) + ((275*month)/9) + day + 1721013;
   return(tmp);
}


void syntax(char *prog)
{
   fprintf(stderr, "Syntax %s [options]  < IMG_0001.pgm > IMG_0001.fits\n", prog);
   fprintf(stderr, " where [options] include:\n");
   fprintf(stderr, "  -hdr     Read metadata file from dcraw; add to header\n");
   fprintf(stderr, "  -rgb     Write R, (G1+G2)/2, B (3 planes)\n");
   fprintf(stderr, "  -bw      Write (R+G1+G2+B)/4 (1 plane)\n");
   fprintf(stderr, "  -r       Write R (1 plane)\n");
   fprintf(stderr, "  -g       Write (G1+G2)/2 (1 plane)\n");
   fprintf(stderr, "  -b       Write B (1 plane)\n");
   fprintf(stderr, "  -g1      Write G1 (1 plane)\n");
   fprintf(stderr, "  -g2      Write G2 (1 plane)\n");
   fprintf(stderr, "  -rggb    Write R, G1, G2, B (4 planes)\n");
   fprintf(stderr, "  -bayer   Write full RGGB (1 plane)\n");
   fprintf(stderr, "  -notrim  Don't trim to image size from metadata file\n");
   fprintf(stderr, "  -div2    Divide the result by 2\n");
   fprintf(stderr, "  -image=x0,y0:x1,y1    Designation of what's image\n");
   fprintf(stderr, "  -bias=x0,y0:x1,y1     Designation of what's bias\n");
   fprintf(stderr, "  -zone=H     Correct reported time to UTC (H=reported-UTC e.g. H=-10 for HST)\n");
   fprintf(stderr, "  -test    Enable a bit of debug info\n");
}
