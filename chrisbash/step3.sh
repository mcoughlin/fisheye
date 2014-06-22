#!/bin/bash
#  9. Obtain and refine camera pointing
# 10. Fit RA, Dec, az, alt, and PA for all exposures

# Where are the extra binaries, bash scripts and monsta scripts?
export PRODIR=~/src/geo

# Where are the night's files? /basedir/sitecam/YYMMDD
dir=/local/jason/s13/geo/aa/130715

filter=M        # Which filter?  B, G, R, M (onochrome)

dt=16.0		# Offset from camera time to center of observation

major=23.8      # tphot streak length [pix] 32sec*(15"/sec)*cos(dec)/(20"/pix)
minor=1.0       # tphot cross-streak width [pix]
minpeak=10      # minimum star brightness to consider

gfrac=0.42	# What's the g fraction of this BW?  Good for BGRI
rfrac=0.58	# What's the r fraction of this BW?  Good for BGRI

im0=0		# First good image
im1=0	        # Last good image

xctr=1408	# Image center pixel (Canon 5D!)
yctr=938	# Image center pixel (Canon 5D!)

d3=-4.5e-05    # lens distortion (Canon 135mm f/2)
scale=19.87    # plate scale ["/pix] (Canon 135mm f/2)

phi=-1.2

# gfrac=0.60	# What's the g fraction of this BW?  Good for BGR
# rfrac=0.40	# What's the r fraction of this BW?  Good for BGR

lng=-117.22858 # Observatory E longitude [deg] (Mt. Joyce)
lat=33.03989   # Observatory latitude [deg] (Mt. Joyce)
elev=45.0      # Observatory elevation [m] (Mt. Joyce)

az=""          # Pre-determined pointing azimuth [deg]
alt=""         # Pre-determined pointing altitude [deg]
pa=""          # Pre-determined pointing position angle [deg]
astrobs=""     # Image selected for astrometry

# Load any new variables from the command line
eval $@
export PATH=$PRODIR:$PATH

# Load variables from the observation file
cd $dir/$filter
prefix=`ls ??1[0-9]????_*.flat | head -1 | tr '_' ' ' | awk '{print $1}'`
source ../$prefix.obs

# Override all variables from the command line
eval $@

echo Appending to the observation parameter file ../$prefix.obs
echo "" >> ../$prefix.obs
echo "# Step3:" >> ../$prefix.obs
echo "unixdatestep3=`date +%s`       # Date step3 was run" >> ../$prefix.obs
echo "datestep3=\"`date`\"       # Date step3 was run" >> ../$prefix.obs
echo "gfrac=$gfrac       # Filter fraction from g" >> ../$prefix.obs
echo "rfrac=$rfrac       # Filter fraction from g" >> ../$prefix.obs
echo "dt=$dt          # Exposure center - recorded time [sec]" >> ../$prefix.obs

# Are we provided with an external az,alt,pa?  If not work from WCS
if [ "$az" == "" -o  "$alt" == "" -o  "$pa" == "" -o "$astrobs" == "" ] ; then

# Analyze the WCS file
  wcs=`ls *.wcs | awk '{print $1}'`
  obs=`basename $wcs .wcs`

# Re-tphot this file so we go reasonably deep!
  ntph=`wc $obs.tfh | awk '{print $1}'`
  echo $ntph
  if [ $ntph -lt 500 ] ; then
    echo tphoting image $obs to deeper limits...
    tphot $obs.flat -trail -aprad 4 -rad 15 -move 15 -min 50 -sig 2 -major $major -minor $minor -phi $phi -out $obs.tfh
  fi

  wc $obs.tfh

  echo ""
  echo Analyzing the WCS for image $obs

  a0=`fitshdr $obs.wcs | grep CRVAL1 | awk '{print $3}'`
  d0=`fitshdr $obs.wcs | grep CRVAL2 | awk '{print $3}'`
  x0=`fitshdr $obs.wcs | grep CRPIX1 | awk '{print $3}'`
  y0=`fitshdr $obs.wcs | grep CRPIX2 | awk '{print $3}'`
  CD11=`fitshdr $obs.wcs | grep CD1_1 | awk '{print $3}'`
  CD12=`fitshdr $obs.wcs | grep CD1_2 | awk '{print $3}'`
  CD21=`fitshdr $obs.wcs | grep CD2_1 | awk '{print $3}'`
  CD22=`fitshdr $obs.wcs | grep CD2_2 | awk '{print $3}'`

  wcscale=`echo $CD11 $CD21 $CD12 $CD22 | awk '{det=$1*$4-$2*$3; s=det>0?sqrt(det):sqrt(-det); scale=3600*s; printf "%.2f", scale}'`

  echo "Plate scale from WCS is $wcscale"
  echo "wcscale=$wcscale     # [arcsec/pix] Plate scale from WCS" >> ../$prefix.obs

  radec=`echo $xctr $yctr | awk -v a0=$a0 -v d0=$d0 -v x0=$x0 -v y0=$y0 \
     -v cd11=$CD11 -v cd21=$CD21 -v cd12=$CD12 -v cd22=$CD22 \
     '{x=$1-x0; y=$2-y0; a=cd11*x+cd12*y+a0; d=cd21*x+cd22*y+d0; printf "%7.2f %7.2f %9.4f %9.4f\n", $1,$2,a,d}'`
  ra=`echo $radec | awk '{print $3}'`
  dec=`echo $radec | awk '{print $4}'`
  pa=`echo $CD11 $CD12 | awk '{print 90+atan2($1,$2)*57.296}'`
  mjd=`grep $obs sniff.out | awk '{print $2}'`

  azalt=`echo $xctr $yctr | cam2sky ra=$ra dec=$dec pa=$pa mjd=$mjd dt=$dt x0=$xctr y0=$yctr lng=$lng lat=$lat elev=$elev azalt`
  az=`echo $azalt | awk '{print $1}'`
  alt=`echo $azalt | awk '{print $2}'`

# Working off of externally provided coords, don't need WCS
else
  obs=$astrobs
  mjd=`grep $obs sniff.out | awk '{print $2}'`
  radec=`echo $xctr $yctr | cam2sky az=$az alt=$alt pa=$pa mjd=$mjd dt=$dt x0=$xctr y0=$yctr lng=$lng lat=$lat elev=$elev`
  ra=`echo $radec | awk '{print $1}'`
  dec=`echo $radec | awk '{print $2}'`
fi

echo ""
echo First pass:
echo ""
echo ra= $ra dec= $dec pa= $pa
echo az= $az alt= $alt mjd= $mjd dt= $dt

# Wide open tolerance in case astrometry.net isn't exactly on
first=`geostars.sh obs=$obs az=$az alt=$alt pa=$pa dt=$dt d3=$d3 scale=$scale xctr=$xctr yctr=$yctr lng=$lng lat=$lat elev=$elev TOL1=0.1 minpeak=$minpeak clean=0`

geostars.sh obs=$obs az=$az alt=$alt pa=$pa dt=$dt d3=$d3 scale=$scale xctr=$xctr yctr=$yctr lng=$lng lat=$lat elev=$elev TOL1=0.1 minpeak=$minpeak clean=0 TEST=1

# exit 0

# Look for the mode in the coord offsets
da=`awk '{da=($3-$1)*100; i=da<0?da-0.5:da+0.5; print int(i)}' $obs.fitlist | sort -g | awk 'BEGIN{prev=1000; best=0; nbest=0} $1!=prev{prev=$1; n=1} $1==prev{n++} n>nbest{best=prev; nbest=n} END{print 0.01*best}'`
dd=`awk '{da=($4-$2)*100; i=da<0?da-0.5:da+0.5; print int(i)}' $obs.fitlist | sort -g | awk 'BEGIN{prev=1000; best=0; nbest=0} $1!=prev{prev=$1; n=1} $1==prev{n++} n>nbest{best=prev; nbest=n} END{print 0.01*best}'`

azinew=`echo $az $alt $da $dd $pa | awk '{printf "%.4f", $1-($3*cos($5/57.296)-$4*sin($5/57.296))/cos($2/57.296)}'`
altnew=`echo $az $alt $da $dd $pa | awk '{printf "%.4f", $2+($3*sin($5/57.296)+$4*cos($5/57.296))}'`

echo ""
echo Offsets in RA and Dec, new pointing az,alt :
echo ""
echo da= $da dd= $dd az= $azinew alt= $altnew

# geostars.sh obs=$obs az=$azinew alt=$altnew pa=$pa dt=$dt d3=$d3 scale=$scale xctr=$xctr yctr=$yctr lng=$lng lat=$lat elev=$elev TEST=1 TOL1=0.02

second=`geostars.sh obs=$obs az=$azinew alt=$altnew pa=$pa dt=$dt d3=$d3 scale=$scale xctr=$xctr yctr=$yctr lng=$lng lat=$lat elev=$elev minpeak=$minpeak TOL1=0.02`

az=`echo $second | awk '{print $5}'`
alt=`echo $second | awk '{print $6}'`
pa=`echo $second | awk '{print $7}'`

echo ""
echo Final az,alt,pa for geostars:
echo ""
echo az=$az alt=$alt pa=$pa dt=$dt

echo "az=$az         # Azimuth of astrobs pointing" >> ../$prefix.obs
echo "alt=$alt        # Altitude of astrobs pointing" >> ../$prefix.obs
echo "pa=$pa         # PA of astrobs pointing" >> ../$prefix.obs
echo "ra=$ra         # RA of astrobs pointing" >> ../$prefix.obs
echo "dec=$dec        # Dec of astrobs pointing" >> ../$prefix.obs
echo "astrmjd=$mjd      # MJD of astrobs pointing" >> ../$prefix.obs

echo ""
echo Now processing all tphot output...

rm -f geostars.flat 2>/dev/null

for ((i=im0; i<=im1; i++)) ; do 
  i4=`printf "%04d" $i`
  ofile=`ls *_$i4.flat`
  obs=`basename $ofile .flat`
  geo=`geostars.sh obs=$obs clean=2 az=$az alt=$alt pa=$pa dt=$dt d3=$d3 scale=$scale xctr=$xctr yctr=$yctr lng=$lng lat=$lat elev=$elev gfrac=$gfrac rfrac=$rfrac minpeak=$minpeak`
  echo "$geo"
# Did we get anything believable?
  nmch=`echo $geo | awk '{print $8}'`
  if [ "$nmch" == "" -o $nmch -lt 100 ] ; then
    echo Got a problem, Houston... Too few stars matched!
    echo Is it time for geomode with large tolerance?
  fi

  echo "$geo" >> geostars.flat
done
