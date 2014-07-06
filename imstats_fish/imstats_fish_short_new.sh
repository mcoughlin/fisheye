#!/bin/bash
# Script to convert a fisheye camera image to FITS and provide some info
# Syntax: ucam_fish.sh obsname [options]
# E.g.    ucam_fish.sh 56754.3993171 [options]
# 
# It converts a CR2 image to monochrome FITS, computes a bunch of
# statistics including solving for the stars
#
# This creates a variety of files.
#
# 56754.3993171.fits.fz = monochrome FITS image
# 56754.3993171.fish    = star list with fisheye information
# 56754.3993171.mch     = star list with photometry information
# 56754.3993171_bw.jpg  = sky image cropped and stretched to show stars
# 56754.3993171_zp.png  = sky image showing a ZP point for each star
# 56754.3993171_sky.png = sky image showing a sky point for each star
#
# Note that PRODIR must point to ql_fish.pro and fishsky.pl
#

set -e
set -x

# sky brightness towards the end of twilight
DAYBRGT=100

# Observatory longitude and latitude (REQUIRED)
#LNG=-155.5763           # MLO
#LAT=19.5362             # MLO

LAT=-30.67833
LNG=-70.73669

# Bootstrap file of x,y,RA,Dec (REQUIRED if pointing not provided)
HA=3.08
DEC=-27.83
AZ=95.55
SCALE=0.06898
PARITY=1
MED=2.566

# Where can hip.dat.dat be found: $CATDIR/hip.dat (REQUIRED)
CATDIR=/lsst/home/coughlin/git-repo/fisheye
PRODIR=/lsst/home/coughlin/git-repo/fisheye

# Other default parameters
VERB=0

# Minimum number of tphot detections to proceed
MINSTAR=100

# Match tolerance [pix]
MCHTOL1=30           # First match tolerance [pix] at m<3 and N<100
MCHTOL=5             # Subsequence match tolerances [pix]
RSIG=3               # Final pruning sigma

MCHTOL=6             # Subsequence match tolerances [pix]
RSIG=5               # Final pruning sigma

#MCHTOL=10             # Subsequence match tolerances [pix]
#RSIG=8               # Final pruning sigma

# Clean up afterwards?
CLEAN=1

if [[ $# -lt 1 ]]
then
    echo "Usage: ucam3_imstats.sh obsname [options]" >&2
    exit 1
fi

obs=$1 ; shift 1

if [[ $obs =~ cr2$ ]]; then
	obs=`basename $obs .cr2`
fi

eval $@

hipcat=$CATDIR/hip.dat

# A few sanity checks
if [[ ! -e ${obs}.cr2 ]] ; then
  echo Observation ${obs}.cr2 does not exist
fi

if [[ ! -e $hipcat ]] ; then
  echo Hipparcos catalog $hipcat does not exist
  exit 1
fi

if [[ ! -e $PRODIR/ql_fish.pro ]] ; then
  echo '$PRODIR/ql_fish.pro' does not exist.  Is PRODIR set?
  exit 1
fi

# Create the FITS file
#raw2fits -bw $obs.cr2

if [[ ! -e $obs.obs ]] ; then
  printf "FILENAME= '%s'     / Observation ID\n" $obs > $obs.obs
fi
etime=`fitshdr $obs.fits | grep EXPTIME | awk '{printf "%.6f", $3}'`

# f-ratio, used to correct for vignetting
frat=`fitshdr $obs.fits | grep APERTURE | awk '{printf "%.6f", $2}'`
frat=4

# Calculate some statistics for the image
#stats=(`monsta $PRODIR/ql_fish.pro $obs`)

# Report all the image stats from ql_ucam.pro
#printf "QL_NSAT = %20d / Number of saturated pixels\n" ${stats[3]} >> $obs.obs
#printf "QL_SKY  = %20.1f / [ADU] median sky\n" ${stats[0]} >> $obs.obs
#printf "QL_RMS  = %20.1f / [ADU] mRMS sky\n" ${stats[1]} >> $obs.obs
#printf "QL_BRGT = %20.3e / [ADU/sec] sky brightness at f/4\n" ${stats[2]} >> $obs.obs

# If the background is too bright it's daytime and don't look for stars
night=`echo ${stats[2]} | awk -v f=$DAYBRGT '{nite=$1>f?0:1; print nite}'`

if [[ $night -eq 1 ]] ; then

# tphot the image, lots o' hardwired parameters...
  TRAD=4 TSIG=2 TMIN=50 TBIAS=2048 TAPRAD=4 
  TAPRAD=8 SKYRAD=40
  #TAPRAD=4

  #tphot $obs.fits -bias $TBIAS -rad $TRAD -sig $TSIG -min $TMIN -aprad $TAPRAD -okfit 0 -chin 100000 -out ${obs}_old.tph

  tphot $obs.fits -bias $TBIAS -rad $TRAD -sig $TSIG -min $TMIN -aprad $TAPRAD -skyrad $SKYRAD -okfit 0 -chin 100000 -out ${obs}_old.tph
  #tphot $obs.fits -bias $TBIAS -rad $TRAD -sig $TSIG -min $TMIN -aprad $TAPRAD -okfit 0 -chin 100000 -out $obs.tph -nxsky 10 -nysky 10

  #grep -v "#" ${obs}_old.tph | awk '{ if (($8/$9) > 10) print $0}' > ${obs}.tph
  grep -v "#" ${obs}_old.tph | awk '{ if (($8/$9) > 1) print $0}' > ${obs}.tph

  #python make_tph_majorminor.py --input ${obs}.tph --output ${obs}.tph

  ntph=`wc $obs.tph | awk '{print $1}'`
# How many stars did we find?
  printf "QL_NTPH = %20d / Number of stars in tphot\n" $ntph >> ${obs}.obs

  if [[ $ntph -ge $MINSTAR ]] ; then
    echo Observation $obs tphot found $ntph stars

# Do the fisheye astrometry fit

# First iteration: MLIM=3, NSTAR=100, matched fit then pruned fit:
    if [[ $VERB -gt 0 ]] ; then echo First iteration: MLIM=3, NSTAR=100 ; fi

    fishiter.sh $obs CATDIR=$CATDIR MLIM=3 NSTAR=100 MCHTOL=$MCHTOL1 HA=$HA DEC=$DEC AZ=$AZ SCALE=$SCALE LNG=$LNG LAT=$LAT FRAT=$frat CLEAN=$CLEAN VERB=$VERB

# Second iteration: MLIM=5, NSTAR=1000, matched fit then pruned fit:
    if [[ $VERB -gt 0 ]] ; then echo Second iteration: MLIM=5, NSTAR=1000 ; fi
    if [[ -e $obs.fish ]] ; then
      parm=(`grep "HA=" $obs.fish`)
      HA=${parm[2]} DEC=${parm[4]} AZ=${parm[6]} SCALE=${parm[8]} QUAD=${parm[10]} CUBE=${parm[12]} CX=${parm[14]} CY=${parm[16]}
      fishiter.sh $obs CATDIR=$CATDIR MLIM=5 NSTAR=1000 MCHTOL=$MCHTOL HA=$HA DEC=$DEC AZ=$AZ SCALE=$SCALE CX=$CX CY=$CY QUAD=$QUAD CUBE=$CUBE LNG=$LNG LAT=$LAT FRAT=$frat CLEAN=$CLEAN VERB=$VERB

# Third iteration: MLIM=6, NSTAR=2000, matched fit then pruned fit:
      if [[ $VERB -gt 0 ]] ; then echo Third iteration: MLIM=6, NSTAR=2000 ; fi
      parm=(`grep "HA=" $obs.fish`)
      HA=${parm[2]} DEC=${parm[4]} AZ=${parm[6]} SCALE=${parm[8]} QUAD=${parm[10]} CUBE=${parm[12]} CX=${parm[14]} CY=${parm[16]}
      fishiter.sh $obs CATDIR=$CATDIR MLIM=6 NSTAR=2000 MCHTOL=5 HA=$HA DEC=$DEC AZ=$AZ SCALE=$SCALE CX=$CX CY=$CY QUAD=$QUAD CUBE=$CUBE RSIG=$RSIG LNG=$LNG LAT=$LAT FRAT=$frat CLEAN=$CLEAN VERB=$VERB

# Insert a few extras that may be useful
      parm=(`grep "HA=" $obs.fish`)
      HA=${parm[2]} DEC=${parm[4]} AZ=${parm[6]} SCALE=${parm[8]} QUAD=${parm[10]} CUBE=${parm[12]} CX=${parm[14]} CY=${parm[16]}
      printf "HA-OBS  = %20.3f / [deg] hour angle\n" $HA >> $obs.obs
      printf "DEC-OBS = %20.3f / [deg] declination\n" $DEC >> $obs.obs
      printf "AZ-OBS  = %20.3f / [deg] azimuth\n" $AZ >> $obs.obs
      printf "SCALE   = %20.5f / [deg/pix] plate scale\n" $SCALE >> $obs.obs

      nstar=`wc $obs.fish | awk '{print $1-1}'`
      zp=(`awk -v dt=$etime 'NR>1{print $11-$18-2.5/log(10)*log(dt)}' $obs.mch | median verb`)
      mu=(`awk 'NR>1{print $11-$18+$20}' $obs.mch | median verb`)

      printf "NSTAR   = %20d / Number of stars in fisheye\n" $nstar >> ${obs}.obs
      printf "ZP-OBS  = %20.2f / [mag] zeropoint 1ADU/sec\n" ${zp[0]} >> $obs.obs
      printf "DZP-OBS = %20.2f / [mag] zeropoint scatter\n" ${zp[1]} >> $obs.obs
      printf "MU-OBS  = %20.2f / [mag/sec^2] sky brightness\n" ${mu[0]} >> $obs.obs

# Make the zeropoint and sky pngs
#      if [[ $nstar -gt $MINSTAR ]] ; then
#        cat > $obs.pl <<EOF
#input $PRODIR/fishsky.mongo
#jpg $obs.zp.fig
#zpmap $obs
#hardcopy wait
#jpg $obs.sky.fig
#skymap $obs
#hardcopy wait
#end
#EOF
#        mongo $obs.pl
#        fig2dev -L png -g black $obs.zp.fig ${obs}_zp.png
#        fig2dev -L png -g black $obs.sky.fig ${obs}_sky.png

#        rm $obs.pl $obs.zp.fig $obs.sky.fig
#      fi
    fi


  else   # Too few tphot stars to attempt a fisheye fit
    printf "NSTAR   = %20d / Number of stars in fisheye\n" 0 >> ${obs}.obs
  fi

else   # Not night time
  printf "QL_NTPH = %20d / Number of stars in tphot\n" 0 >> ${obs}.obs

# Make a smaller picture if daytime
#  convert -crop 3100x2560+375+0 -scale 25% $obs.jpg ${obs}_clr.jpg
fi


# Make a smaller picture increasing constrast for night time
#mdcraw -e -c $obs.cr2 | convert - -crop 4656x3840+540+0 -scale 20\% -contrast-stretch 2\%x0.1\% -gravity south -fill white -undercolor '#00000080' -annotate +0+5 "$obs" $(basename $obs .cr2)_clr.jpg


# Insert the .obs info  into FITS file
#fitshdr $obs.fits $obs.obs

# fpack the FITS file and delete original
#fpack -D -Y $obs.fits


# Copy the file to our weather page.
#if [[ "$SUMMIT_CRON" ]]; then
#	cd /home/atlas/fisheye
#	convert latest.jpg -resize 400 latest400.jpg
#	convert latest.jpg -resize 1000 latest1000.jpg
#	scp latest400.jpg latest1000.jpg atlas@atlas.ifa-instruments.org:/var/www/weather/mlo/
#fi

#exit 0
