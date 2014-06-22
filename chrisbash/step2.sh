#!/bin/bash
# 5. Make a set of flatfields (immedian)
# 6. Map hot pixels.  Just copy canon5d_??.hotpix from 2013_07_09.
# 7. Check the time in the headers against GPS time of exposure start:
# 8. Flatten and tphot all the images (flatten.pro, tphot)

# Where are the extra binaries, bash scripts and monsta scripts?
export PRODIR=~/src/geo

# Where are the night's files? /basedir/sitecam/YYMMDD
dir=/local/jason/s13/geo/aa/130715

im0=755		# First good image
im1=1410	# Last good image

major=23.8      # tphot streak length [pix] 32sec*(15"/sec)*cos(dec)/(20"/pix)
minor=1.0       # tphot cross-streak width [pix]

anetcount=300   # Count of stars for astrometry.net

filter=M        # Which filter?  B, G, R, M (onochrome)

dotphot=1       # Do tphot?

skipdone=0       # Skip files that already have flat and tphot

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
echo "# Step2:" >> ../$prefix.obs
echo "unixdatestep2=`date +%s`       # Date step2 was run" >> ../$prefix.obs
echo "datestep2=\"`date`\"       # Date step2 was run" >> ../$prefix.obs
echo "filter=$filter         # Bayer filter processed" >> ../$prefix.obs
echo "im0=$im0         # First sequence number" >> ../$prefix.obs
echo "im1=$im1         # Last sequence number" >> ../$prefix.obs

# echo Copying hot pixel mask from 2013_07_09...
# echo
# echo BEWARE, copying hot pixel mask is really bogus, needs HELP!...
# echo Probably best to do this with a median of a bunch of darks.
# echo
# echo
# echo
# cp ../../2013_07_09/BW/canon5d_??.hotpix .

echo Creating median images by hundreds for flatfields...

h0=`echo $im0 | awk '{print int($1/100)}'`
h1=`echo $im1 | awk '{print int($1/100)}'`

for ((i=$h0; i<=$h1; i++)) ; do
  twodigit=`printf "%02d" $i`
  files=(`ls *_${twodigit}??.fits`)
  nfile=${#files[@]} 
  mid=${files[$((nfile/2))]} 
  bias=`fitshdr $mid | grep BIAS | awk '{print $3}'`
  if [ $skipdone -eq 1 -a -e  median.${twodigit}00 ] ; then
    echo median.${twodigit}00 done, skipping...
    continue
  fi
  immedian -out median.${twodigit}00 -norm 600 -bias $bias ${files[@]}
  echo median.${twodigit}00 created
done

mid=`echo $im0 $im1 | awk '{printf "%04d", ($1+$2)/2}'`
midfile=`ls *_$mid.fits`
obs=`basename $midfile .fits`
astrobs=$obs

echo Doing tphot on the mid image $obs to get the trail angle...

med=`echo $im0 $im1 | awk '{printf "%04d", 100*int(($1+$2)/200)}'`

monsta $PRODIR/flatten.pro $obs.{fits,flat} median.$med
tphot $obs.flat -aprad 4 -rad 15 -move 15 -min 500 -sig 3 -out $obs.tfh

medphi=`median col=12 $obs.tfh verb`
phi=`echo $medphi | awk '{phi=$1<90?$1:$1-180; printf "%.1f", phi}'`

echo "astrobs=$obs     # Obs selected for astrometry.net" >> ../$prefix.obs
echo "phi=$phi         # Median angle of star streaks" >> ../$prefix.obs
echo "Median angle $phi"

# Remove this non-forced PSF tfh file
rm $obs.tfh

echo Now flattening and tphoting all images...

for ((i=im0; i<=im1; i++)) ; do 
  num=`echo $i | awk '{printf "%04d", $1}'`
  midfile=`ls *_$num.fits`
  obs=`basename $midfile .fits`
  if [ $skipdone -eq 1 -a -e $obs.flat -a -e $obs.tfh ] ; then
    echo $obs already done, skipping
    continue
  fi
  flat=`echo $i | awk '{printf "median.%04d", 100*int($1/100)}'`
  skynoise=`monsta $PRODIR/flatten.pro $obs.fits $obs.flat $flat`
  if [ $dotphot -eq 1 ] ; then
# Cut off search at 2.5 sigma
    cutoff=`echo $skynoise | awk -v nsig=2.5 '{print $4*nsig}'`
    tphot $obs.flat -trail -aprad 4 -rad 15 -move 15 -min $cutoff -sig 2 -major $major -minor $minor -phi $phi -out $obs.tfh
    echo "$skynoise  `wc $obs.tfh`"
  else
    echo $obs.flat flattened
  fi
done

# Write a set of x,y from a selected to submit to astrometry.net
sort -g -k 8 -r $astrobs.tfh | uniq | grep -v skyfit | awk '$5>100{printf "%s,%s\n", $1,$2}' | head -$anetcount > $astrobs.txt
cp $astrobs.txt /tmp

echo $astrobs.txt copied to /tmp for convenience in uploading
echo Please submit to nova.astrometry.net with image width 8-24 deg.
echo Download the wcs.fits and rename $astrobs.wcs
