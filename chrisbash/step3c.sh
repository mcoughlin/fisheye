#!/bin/bash
# 10b.  Calculate how this filter depends on gr and determine gfrac and rfrac

# Where are the extra binaries, bash scripts and monsta scripts?
export PRODIR=~/src/geo

# Where are the night's files? /basedir/sitecam/YYMMDD
dir=/local/jason/s13/geo/aa/130715

filter=M        # Which filter?  B, G, R, M (onochrome)

bluest=0.0
reddest=1.2
clrpair=gr

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
echo "# Step3c:" >> ../$prefix.obs
echo "unixdatestep3c=`date +%s`       # Date step3c was run" >> ../$prefix.obs
echo "datestep3c=\"`date`\"       # Date step3c was run" >> ../$prefix.obs

nmch=`ls *_????.mch | wc | awk '{print $1}'`
if [ $nmch -lt 5 ] ; then
  echo Fewer than 5 mch files?  Have you run step3 and geostars?
  exit 1
fi

# Loop over (g-r) color in steps of 0.1, dump out median r band ZP
# Output is (g-r)   r+2.5log(fit_flux)   RMS   Npt
if [ $clrpair == gr ] ; then
  gr0=`echo $bluest | awk '{printf "%d", $1*10}'`
  gr1=`echo $reddest | awk '{printf "%d", $1*10+0.5}'`
  for ((i=$gr0; i<=$gr1; i++)) ; do
    clr=`echo $i | awk '{printf "%4.1f", 0.1*$1}'`
    m=`cat *_????.mch | awk -v clr=$i '{f=$5; df=$6; g=$11; r=$12; i=$13; ok=(f>0 && df<0.1*f)} ok && g-r>=0.1*(clr-0.5) && g-r<0.1*(clr+0.5){printf "%7.3f\n", r+2.5*log(f)/log(10)}' | median verb`
    echo "$clr $m"
  done > medzp.gr

elif [ $clrpair == ri ] ; then
  ri0=-2
  ri1=`echo $reddest | awk '{printf "%d", $1*10+0.5}'`
  for ((i=$ri0; i<=$ri1; i++)) ; do
    clr=`echo $i | awk '{printf "%4.1f", 0.1*$1}'`
    m=`cat *_????.mch | awk -v clr=$i '{f=$5; df=$6; g=$11; r=$12; i=$13; ok=(f>0 && df<0.1*f)} ok && r-i>=0.1*(clr-0.5) && r-i<0.1*(clr+0.5){printf "%7.3f\n", r+2.5*log(f)/log(10)}' | median verb`
    echo "$clr $m"
  done > medzp.ri

elif [ $clrpair == gi ] ; then
  gi0=-2
  gi1=`echo $reddest | awk '{printf "%d", $1*10+0.5}'`
  for ((i=$gi0; i<=$gi1; i++)) ; do
    clr=`echo $i | awk '{printf "%4.1f", 0.1*$1}'`
    m=`cat *_????.mch | awk -v clr=$i '{f=$5; df=$6; g=$11; r=$12; i=$13; ok=(f>0 && df<0.1*f)} ok && g-i>=0.1*(clr-0.5) && g-i<0.1*(clr+0.5){printf "%7.3f\n", g+2.5*log(f)/log(10)}' | median verb`
    echo "$clr $m"
  done > medzp.gi
fi

if [ $clrpair == gr ] ; then
  lin=`linear medzp.gr 2 2 0 1`
  echo r-minst fit as a function of "(g-r)": $lin
elif [ $clrpair == ri ] ; then
  lin=`linear medzp.ri 2 2 0 1`
  echo r-minst fit as a function of "(r-i)": $lin
elif [ $clrpair == gi ] ; then
  lin=`linear medzp.gi 2 2 0 1`
  echo g-minst fit as a function of "(g-i)": $lin
fi

# Get the sky brightness
skyadu=`grep -v skyval *.tfh | median col=4 verb | awk '{print $1}'`
apfit=`median col=15 geostars.flat verb | awk '{print $1}'`

rms=`echo $lin | awk '{printf "%.3f", $3}'`
rzp=`echo $lin | awk '{print $9}'`
clr=`echo $lin | awk '{print $10}'`
zp=`echo $rzp | awk '{printf "%.3f", $1}'`

# Pick a typical color for the sky...
skymag=`echo $rzp $clr $skyadu $scale $apfit | awk '{printf "%.2f", $1-$5+0.4*$2-2.5*log($3/$4/$4)/log(10)}'`

# echo $rzp $clr $skyadu $scale $skymag $apfit

if [ $clrpair == gr ] ; then
#  r-minst = 15.653 - 0.410*(g-r)
#  minst = -15.653 + 0.41*g + (1-0.41)*r
  rfrac=`echo $clr | awk '{printf "%.2f", 1+$1}'`
  gfrac=`echo $clr | awk '{printf "%.2f", -$1}'`
  echo "minst = -$zp + $gfrac*g + $rfrac*r   RMS ~ $rms"
elif [ $clrpair == ri ] ; then
#  r-minst = 15.653 + 0.363*(r-i)
#  minst = r -15.653 - 0.36*r - -0.36*i
#  minst = -15.653 + 0.64*r + 0.36*i
  rfrac=`echo $clr | awk '{printf "%.2f", 1-$1}'`
  ifrac=`echo $clr | awk '{printf "%.2f", $1}'`
  echo "minst = -$zp + 0*g + $rfrac*r + $ifrac*i   RMS ~ $rms"
elif [ $clrpair == gi ] ; then
#  g-minst = 15.653 - 0.070*(g-i)
#  minst = g - 15.653 + 0.070*(g-i)
#  minst = -15.653 + 1.07*g + -0.07*i
  gfrac=`echo $clr | awk '{printf "%.2f", 1-$1}'`
  ifrac=`echo $clr | awk '{printf "%.2f", $1}'`
  echo "minst = -$zp + $gfrac*g + 0*r + $ifrac*i   RMS ~ $rms"
fi
echo "skymag = $skymag"

echo "gfrac=$gfrac       # Filter fraction from g" >> ../$prefix.obs
echo "rfrac=$rfrac       # Filter fraction from r" >> ../$prefix.obs
echo "zp=$zp         # Star fit zeropoint" >> ../$prefix.obs
echo "skyadu=$skyadu       # Sky brightness [ADU/pix]" >> ../$prefix.obs
echo "skymag=$skymag       # Sky brightness [ABmag/arcsec]" >> ../$prefix.obs

# Actually, keep this around...
# rm medzp.gr
