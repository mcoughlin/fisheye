#!/bin/bash
# Script to robustly find stars in a tphot data file
# Syntax: geomode.sh \
#  obs=20130629_1000
#  [ alt=51.2951
#    az=170.0690
#    pa=8.7437
#    d3=-4.5e-05
#    scale=19.87
#    dt=16 ]

obs=20130629_0000   # (REQUIRED)

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

dt=16          # Difference between exposure center and recorded time

TOL1=0.1       # First lock on tolerance
mlim=10	       # Limiting mag in catalog (brighten in plane?)
degbin=0.01    # Bin size [deg] for mode

clean=0
TEST=0

eval $@

# Which extension?  tph for unflattened or tfh for flattened.
ext=tfh

if [ ! -e $obs.$ext ] ; then
  echo No $ext output?  Cannot proceed, sorry.
  exit 1
fi

# What's the time of the exposure?
mjd=`fitshdr $obs.flat | grep MJD-OBS | awk '{print $3}'`

# Get a preliminary map of RA, Dec
awk -v pk=$minpeak 'NR>1 && $5>pk{print $0}' $obs.$ext | cam2sky mjd=$mjd dt=$dt alt=$alt az=$az pa=$pa d3=$d3 scale=$scale x0=$xctr y0=$yctr lng=$lng lat=$lat elev=$elev > $obs.rd
  
# Where are we, approximately?
a0=`awk '{s+=$1}END{printf "%.2f", s/NR}' $obs.rd`
d0=`awk '{s+=$2}END{printf "%.2f", s/NR}' $obs.rd`

# Get the star catalog
if [ ! -e $obs.gscat ] ; then
  gscat $a0 $d0 $camwid $camhi filt=r mlim=$mlim > $obs.gscat
fi

if [ $TEST -gt 0 ] ; then
  echo "Geomode: Obs= $obs  MJD= $mjd   a0= $a0   d0= $d0" 
fi

# Match up the initial star coords with the catalog.
# Center approximately in order to decouple offset from rotation
colmerge 1,2 $obs.rd 1,2 $obs.gscat -tol $TOL1 -nobar | awk -v a0=$a0 -v d0=$d0 '{printf "%9.5f %9.5f %9.5f %9.5f\n", $1-a0,$2-d0,$3-a0,$4-d0}' > $obs.fitlist
npt=`wc $obs.fitlist | awk '{print $1}'`

# Look for the mode in the coord offsets, bins of degbin (0.01 deg)
apeak=`awk -v bin=$degbin '{da=($1-$3)/bin; i=da<0?da-0.5:da+0.5; print int(i)}' $obs.fitlist | sort -g | awk -v bin=$degbin 'BEGIN{prev=1000; best=0; nbest=0} $1!=prev{prev=$1; n=1} $1==prev{n++} n>nbest{best=prev; nbest=n} END{printf "%.3f %d", bin*best, nbest}'`
dpeak=`awk -v bin=$degbin '{da=($2-$4)/bin; i=da<0?da-0.5:da+0.5; print int(i)}' $obs.fitlist | sort -g | awk -v bin=$degbin 'BEGIN{prev=1000; best=0; nbest=0} $1!=prev{prev=$1; n=1} $1==prev{n++} n>nbest{best=prev; nbest=n} END{printf "%.3f %d", bin*best, nbest}'`

da=`echo $apeak | awk '{print $1}'`
na=`echo $apeak | awk '{print $2}'`
dd=`echo $dpeak | awk '{print $1}'`
nd=`echo $dpeak | awk '{print $2}'`
ntot=`wc $obs.fitlist | awk '{print $1}'`

asig=`echo $na $ntot $degbin $TOL1 | awk '{meanocc=$2/(2*$4/$3); m=meanocc>0?meanocc:1; printf "%.1f", $1/m}'`
dsig=`echo $nd $ntot $degbin $TOL1 | awk '{meanocc=$2/(2*$4/$3); m=meanocc>0?meanocc:1; printf "%.1f", $1/m}'`

azinew=`echo $az $alt $da $dd $pa | awk '{printf "%.4f", $1+($3*cos($5/57.296)-$4*sin($5/57.296))/cos($2/57.296)}'`
altnew=`echo $az $alt $da $dd $pa | awk '{printf "%.4f", $2-($3*sin($5/57.296)+$4*cos($5/57.296))}'`

if [ $TEST -gt 0 ] ; then
  echo ""
  echo Offsets in RA and Dec, new pointing azi,alt :
  echo ""
  echo dra=$da ddec=$dd az=$azinew alt=$altnew ntot=$ntot rasig=$asig decsig=$dsig
fi

# echo "  Image     MJD         RA     Dec   dRA  dDec  sigR sigD   Azi       Alt      PA  naz  nalt  ntot"
printf "%s %10.6f %6.2f %6.2f %5.2f %5.2f %4.1f %4.1f %8.4f %8.4f %8.4f %4d %4d %4d\n" $obs $mjd $a0 $d0 $da $dd $asig $dsig $azinew $altnew $pa $na $nd $ntot

# Clean up
if [ $clean -gt 0 ] ; then
  rm $obs.rd $obs.fitlist
fi
if [ $clean -gt 1 ] ; then
  rm $obs.gscat
fi

exit 0
