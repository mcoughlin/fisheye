#!/bin/bash
# 10a. Fix up the times in geostars.flat, create geostars.dat
# This is the place that the entire run of times in geostars.flat
# is examined and the times are corrected assuming that the frames are
# taken on a strict (but non-integer) cadence.
# The pointing and star match files are updated.

# Where are the extra binaries, bash scripts and monsta scripts?
export PRODIR=~/src/geo

# Where are the night's files? /basedir/sitecam/YYMMDD
dir=/local/jason/s13/geo/aa/130715

filter=M        # Which filter?  B, G, R, M (onochrome)

xctr=1408	# Image center pixel (Canon 5D!)
yctr=938	# Image center pixel (Canon 5D!)
minpeak=10      # minimum star brightness to consider

d3=-4.5e-05    # lens distortion
scale=19.87    # plate scale ["/pix]

dt=16.0		# Offset from camera time to center of observation

gfrac=0.42	# What's the g fraction of this BW?  Good for BGRI
rfrac=0.58	# What's the r fraction of this BW?  Good for BGRI

lng=-117.22858 # Observatory E longitude [deg] (Mt. Joyce)
lat=33.03989   # Observatory latitude [deg] (Mt. Joyce)
elev=45.0      # Observatory elevation [m] (Mt. Joyce)

# Load any new variables from the command line
eval $@
export PATH=$PRODIR:$PATH

# Load variables from the observation file
cd $dir/$filter
prefix=`ls *.fits | head -1 | tr '_' ' ' | awk '{print $1}'`
source ../$prefix.obs

# Override all variables from the command line
eval $@

echo Appending to the observation parameter file ../$prefix.obs
echo "" >> ../$prefix.obs
echo "# Step3a:" >> ../$prefix.obs
echo "unixdatestep3a=`date +%s`       # Date step3a was run" >> ../$prefix.obs
echo "datestep3a=\"`date`\"       # Date step3a was run" >> ../$prefix.obs

# Check that we have geostars.flat!
if [ ! -e geostars.flat ] ; then
  echo "geostars.flat does not exist.  Please run step3.sh (or geostars)"
  exit 1
fi

# Fit time interval between frames [sec]
dtlin=`awk 'NR==1{mjd0=$2} {printf "%d %.6f\n", NR-1, ($2-mjd0)*86400}' geostars.flat | linear - 2 2 0 1`
dtrms=`echo $dtlin | awk '{printf "%d", $3}'`

# Whoopsie, must be that stupid 130629 glitch
if [ $dtrms -gt 0 ] ; then
  echo Bad fit: $dtlin
  dtlin=`awk 'NR==1{mjd0=$2} NR>100{printf "%d %.6f\n", NR-1, ($2-mjd0)*86400}' geostars.flat | linear - 2 2 0 1`
fi
dtrms=`echo $dtlin | awk '{print $3}'`
t0obs=`echo $dtlin | awk '{print $9}'`
dtobs=`echo $dtlin | awk '{print $10}'`

echo "Time(frame) fit:  t0= $t0obs tobs= $dtobs rms= $dtrms"

echo "t0=$t0obs       # Time offset wrt camera" >> ../$prefix.obs
echo "tcad=$dtobs       # Observation cadence" >> ../$prefix.obs
echo "trms=$dtrms          # RMS [sec]" >> ../$prefix.obs

# Update geostars.flat into geostars.dat
rm -f geostars.dat 2>/dev/null
mjd0=`head -1 geostars.flat | awk '{print $2}'`
n=0
exec 10<geostars.flat
while read LINE<& 10 ; do
  mjd=`echo $LINE | awk '{print $2}'`
  dtcorr=`echo $mjd $mjd0 $t0obs $n $dtobs $dt | awk '{dtfit=($2-$1)*86400+$3+$4*$5; dtnint=dtfit>0?int(dtfit+0.5):int(dtfit-0.5); crap=0} dtnint>3 || dtnint<-3{crap=dtnint} {printf "%.3f", dtfit-crap+$6}'`
  let n++

  obs=`echo $LINE | awk '{print $1}'`
  az=`echo $LINE | awk '{print $5}'`
  alt=`echo $LINE | awk '{print $6}'`
  pa=`echo $LINE | awk '{print $7}'`
#  echo "$obs $n $mjd $az $alt $pa $dtcorr"

  geo=`geostars.sh obs=$obs clean=2 az=$az alt=$alt pa=$pa dt=$dtcorr d3=$d3 scale=$scale xctr=$xctr yctr=$yctr lng=$lng lat=$lat elev=$elev gfrac=$gfrac rfrac=$rfrac minpeak=$minpeak`
  echo "$geo"
  echo "$geo" >> geostars.dat
done
