#!/bin/bash
# Script to find stars and fit them in a trailed image
# Syntax: geostars.sh \
#  obs=20130629_1000
#  [ alt=51.2951
#    az=170.0690
#    pa=8.7437
#    d3=-4.5e-05
#    scale=19.87
#    tphot={0:1} ]

obs=20130629_0000   # (REQUIRED)

tphot=0        # Run tphot?
trailen=23.5   # (fixed) trail length ~(texp)*(15"/s)*cos(dec)/(19.87"/pix)
trailwid=1.0   # (fixed) trail width
trailphi=-9.0  # (fixed) trail phi
minpeak=100    # minimum peak above sky of trail to keep

xctr=1408	# Image center pixel (Canon 5D!)
yctr=938	# Image center pixel (Canon 5D!)

camwid=10.0    # Deg FOV for gscat
camhi=8.0      # Deg FOV for gscat

az=189.926     # camera azimuth [130629]
alt=51.2951    # camera altitude [130629]
pa=8.7437      # camera PA [130629]
d3=-4.5e-05    # lens distortion
scale=19.87    # plate scale ["/pix]

lng=-117.22858 # Observatory E longitude [deg] (Mt. Joyce)
lat=33.03989   # Observatory latitude [deg] (Mt. Joyce)
elev=45.0      # Observatory elevation [m] (Mt. Joyce)

gfrac=0.6      # fraction of g in Canon BW 
rfrac=0.4      # fraction of r in Canon BW 

dt=16          # Difference between exposure center and recorded time

TOL1=0.02      # First lock on tolerance
mlim1=10       # First limiting catalog magnitude
TOL2=0.005     # Final lock on tolerance
mlim=10        # Final limiting catalog magnitude
# VIG=0.014      # Vignetting correction (unflattened) [mag/deg^2]
VIG=0.0        # Vignetting correction (flattened) [mag/deg^2]

clean=1
TEST=0

if [ "$1" == "HDR" ] ; then
  echo "   Image           MJD       RA     Dec      Azi     Alt    PA     Nstar  dRA    +/-   dDec    +/-    ZP     +/-   ap-fit  +/-"
  exit 0
fi

eval $@

# Which extension?  tph for unflattened or tfh for flattened.
if [ $VIG == "0.0" ] ; then ext=tfh ; else ext=tph ; fi

# Process image with tphot to get a star list?
if [ $tphot -eq 1 ] ; then
  bias=`fitshdr $obs.flat | grep BIAS | awk '{print $3}'`
  if [ $VIG == "0.0" ] ; then
    tphot $obs.flat -trail -aprad 4 -rad 15 -move 15 -min 50 -sig 2 -major $trailen -minor $trailwid -phi $trailphi -out $obs.$ext
  else
    tphot $obs.fits -trail -rad 15 -move 15 -min 50 -sig 2 -bias $bias \
       -major $trailen -minor $trailwid -phi $trailphi -nxsky 25 -nysky 20 \
       -residsky -out $obs.$ext -resid $obs.resid
  fi
fi

if [ ! -e $obs.$ext ] ; then
  echo No $ext output?  Cannot proceed, sorry.
  exit 1
fi

# What's the time of the exposure?
mjd=`fitshdr $obs.flat | grep MJD-OBS | awk -v dt=$dt '{printf "%.7f", $3+dt/86400}'`

# I'm pissed, just insist on geomode to tune up the az,alt
if [ $TEST -gt 0 ] ; then
  geomode.sh obs=$obs az=$az alt=$alt pa=$pa lng=$lng lat=$lat elev=$elev d3=$d3 scale=$scale xctr=$xctr yctr=$yctr dt=$dt clean=0 TEST=1
fi
geom=`geomode.sh obs=$obs az=$az alt=$alt pa=$pa lng=$lng lat=$lat elev=$elev d3=$d3 scale=$scale xctr=$xctr yctr=$yctr dt=$dt clean=0`
az=`echo $geom | awk '{print $9}'`
alt=`echo $geom | awk '{print $10}'`

# Get a preliminary map of RA, Dec
awk -v pk=$minpeak 'NR>1 && $5>pk{print $0}' $obs.$ext | cam2sky mjd=$mjd dt=0 alt=$alt az=$az pa=$pa d3=$d3 scale=$scale x0=$xctr y0=$yctr lng=$lng lat=$lat elev=$elev > $obs.rd
  
# Where are we, approximately?
a0=`awk '{s+=$1}END{printf "%.2f", s/NR}' $obs.rd`
d0=`awk '{s+=$2}END{printf "%.2f", s/NR}' $obs.rd`

# Get the star catalog
if [ ! -e $obs.gscat ] ; then
  gscat $a0 $d0 $camwid $camhi filt=r mlim=$mlim1 > $obs.gscat
fi

if [ $TEST -gt 0 ] ; then
  echo "Obs= $obs  MJD= $mjd   a0= $a0 d0= $d0  az= $az  alt=$alt"
fi

# Match up the initial star coords with the catalog.
# Center approximately in order to decouple offset from rotation
colmerge 1,2 $obs.rd 1,2 $obs.gscat -tol $TOL1 -nobar | awk -v a0=$a0 -v d0=$d0 '{printf "%9.5f %9.5f %9.5f %9.5f\n", $1-a0,$2-d0,$3-a0,$4-d0}' > $obs.fitlist
npt=`wc $obs.fitlist | awk '{print $1}'`

if [ $TEST -gt 0 ] ; then echo "First match npt= $npt" ; fi

# Uh oh, nothing good. Try to get in range with geomode
if [ $npt -lt 500 ] ; then
  geom=`geomode.sh obs=$obs az=$az alt=$alt pa=$pa dt=$dt lng=$lng lat=$lat elev=$elev d3=$d3 scale=$scale xctr=$xctr yctr=$yctr`
  if [ $TEST -gt 0 ] ; then 
    echo ""
    echo "Too few points, fall back on geomode.sh to update az,alt:"
    echo "$geom"
    echo ""
  fi
  az=`echo $geom | awk '{print $9}'`
  alt=`echo $geom | awk '{print $10}'`

# Update the map of RA, Dec
  awk -v pk=$minpeak 'NR>1 && $5>pk{print $0}' $obs.$ext | cam2sky mjd=$mjd dt=0 alt=$alt az=$az pa=$pa d3=$d3 scale=$scale x0=$xctr y0=$yctr lng=$lng lat=$lat elev=$elev > $obs.rd

# Center approximately in order to decouple offset from rotation
  colmerge 1,2 $obs.rd 1,2 $obs.gscat -tol $TOL1 -nobar | awk -v a0=$a0 -v d0=$d0 '{printf "%9.5f %9.5f %9.5f %9.5f\n", $1-a0,$2-d0,$3-a0,$4-d0}' > $obs.fitlist
fi

offset="`fitrot -if $obs.fitlist -rotate -prune -$((npt/10))`"

if [ $TEST -gt 0 ] ; then
  fitrot -if $obs.fitlist -scale -prune -$((npt/10))
fi

# Offsets in coordinates
dra=`echo $offset | awk '{print $2}'`
ddec=`echo $offset | awk '{print $3}'`
dpa=`echo $offset | awk '{print $5}'`

if [ $TEST -gt 0 ] ; then
  echo "npt= $npt   dRA= $dra   dDec= $ddec   dPA= $dpa"
fi

# Compute a better alt, az, pa  (FIXME! CHEESY!)
az2=`echo $az $alt $dra $ddec $pa | awk '{printf "%.5f", $1-($3*cos($5/57.296)-$4*sin($5/57.296))/cos($2/57.296)}'`
alt2=`echo $az $alt $dra $ddec $pa | awk '{printf "%.5f", $2+($3*sin($5/57.296)+$4*cos($5/57.296))}'`
pa2=`echo $pa $dpa | awk '{printf "%.4f", $1-$2}'`

if [ $TEST -gt 0 ] ; then
  echo "Az= $az $az2   Alt= $alt $alt2   PA= $pa $pa2"
fi

# Get a better map of RA, Dec
awk -v pk=$minpeak 'NR>1 && $5>pk{print $0}' $obs.$ext | cam2sky mjd=$mjd dt=0 alt=$alt2 az=$az2 pa=$pa2 d3=$d3 scale=$scale x0=$xctr y0=$yctr lng=$lng lat=$lat elev=$elev > $obs.rd2
  
# Merge it with the catalog and test its fit
colmerge 1,2 $obs.rd2 1,2 $obs.gscat -tol $TOL2 -nobar | awk -v a0=$a0 -v d0=$d0 '{printf "%9.5f %9.5f %9.5f %9.5f\n", $1-a0,$2-d0,$3-a0,$4-d0}' > $obs.fitlist2
npt=`wc $obs.fitlist2 | awk '{print $1}'`
offset="`fitrot -if $obs.fitlist2 -rotate -prune -$((npt/10))`"

if [ $TEST -gt 0 ] ; then
  echo "$offset"
fi

# Extract x,y peak,dpeak flux,dflux from tph file, merge with RA,Dec
awk -v pk=$minpeak 'NR>1 && $5>pk{printf "%8.2f %8.2f %8.1f %6.1f %8.0f %8.0f\n",$1,$2,$5,$6,$8,$9}' $obs.$ext | colmerge 0 $obs.rd2 0 - -nobar > $obs.radec

# Get the star catalog again if we want to push deeper
if [ $mlim -gt $mlim1 ] ; then
  gscat $a0 $d0 $camwid $camhi filt=r mlim=$mlim > $obs.gscat2
else
  ln -s -f $obs.gscat $obs.gscat2
fi

# Merge tph radec file with the catalog
echo "    RA[tph]     Dec[tph]       x         y       peak   dpeak    flux     dflux    RA[cat]   Dec[cat]    g      r      i      z      y" > $obs.mch
colmerge 1,2 $obs.radec 1,2 $obs.gscat -tol $TOL2 -nobar > $obs.mch

# Some statistics from the match
medra=`awk '{print $1-$9}' $obs.mch | median verb`
medec=`awk '{print $2-$10}' $obs.mch | median verb`
# Canon 5D is 2816x1876;
# Canon standard BW ~ 0.6*g+0.4*r+0.0*i
# Canon 135mm lens vignettes at about 0.014*[deg]^2
medzp=`awk -v vig=$VIG -v g=$gfrac -v r=$rfrac '{i=1-g-r; m=g*$11+r*$12+i*$13; rad=(($3-1408)*($3-1408)+($4-938)*($4-938))/180; dm=vig*r2; minst=-2.5*log($5)/log(10); print m-minst+dm}' $obs.mch | median verb`
apfit=`awk '{ap=$7>0?$7:1; map=-2.5*log(ap)/log(10); mpk=-2.5*log($5)/log(10); print map-mpk}' $obs.mch | median verb`

if [ $TEST -gt 1 ] ; then
  echo $medra
  echo $medec
  echo $medzp
  echo $apfit
fi

# Tell us about it
echo $obs $mjd $a0 $d0 $az2 $alt2 $pa2 $medra $medec $medzp $apfit | \
   awk '{printf "%s  %13.7f %6.2f %6.2f  %9.5f %8.5f %6.4f  %4d %6.2f %6.2f %6.2f %6.2f %7.3f %6.3f %7.3f %6.3f\n", \
    $1,$2,$3,$4,$5,$6,$7,$10,$8*3600,$9*3600,$11*3600,$12*3600,$14,$15,$17,$18}'

# Clean up
if [ $clean -gt 0 ] ; then
  rm $obs.rd $obs.rd2 $obs.fitlist $obs.fitlist2
fi
if [ $clean -gt 1 ] ; then
  rm $obs.gscat $obs.gscat2 $obs.radec
fi

exit 0
