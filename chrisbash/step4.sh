#!/bin/bash
# 11. Suppress stars from each image and detect sats (stardiff.pro, tphot)
# 12. Link detections into satellites (tracks)

# Where are the extra binaries, bash scripts and monsta scripts?
export PRODIR=~/src/geo

# Where are the night's files? /basedir/sitecam/YYMMDD
dir=/local/jason/s13/geo/aa/130715

filter=M        # Which filter?  B, G, R, M (onochrome)

im0=755		# First good image
im1=1410	# Last good image

d3=-4.5e-05    # lens distortion
scale=19.87    # plate scale ["/pix]

hotpixmask=canon5d_jt.hotpix    # Hot pixel mask

dodiff=1       # regenerate the stardiff files?
deldiff=1       # delete the stardiff files?

skipdone=0      # Skip files that already have diff and tphot

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
echo "# Step4:" >> ../$prefix.obs
echo "unixdatestep4=`date +%s`       # Date step4 was run" >> ../$prefix.obs
echo "datestep4=\"`date`\"       # Date step4 was run" >> ../$prefix.obs
echo "hotpixmask=$hotpixmask       # Name of hot pixel mask" >> ../$prefix.obs

# Get the leading part of the filename, assume of the form XXX_nnnn
pre=`ls *.wcs | tr '_' ' ' | awk '{print $1}'`

# Create star-diffed images by subtracting prev and post; tphot result
if [ $dodiff -eq 1 ] ; then
  for ((i=$((im0+1)); i<=$((im1-1)); i++)) ; do 
    this=`printf "%04d" $i`
    prev=`printf "%04d" $((i-1))`
    next=`printf "%04d" $((i+1))`

    if [ $skipdone -eq 1 -a -e ${pre}_$this.td ] ; then
      echo ${pre}_$this.td already done, skipping
      continue
    fi

    skynoise=`monsta $PRODIR/stardiff.pro ${pre}_{$this,$prev,$next}.flat ${pre}_${this}.diff $hotpixmask`
# Cut off search at 4 sigma
    cutoff=`echo $skynoise | awk -v nsig=4 '{print $4*nsig}'`
    tphot ${pre}_${this}.diff -min $cutoff -sig 2 -aprad 4 -okfit 0 -out ${pre}_$this.td
    echo "$skynoise  `wc ${pre}_$this.td`"
    if [ $deldiff == 1 ] ; then rm ${pre}_${this}.diff ; fi
  done
fi

# Extract the tracks from all these output files
tracks meta=geostars.dat d3=$d3 scale=$scale ${pre}_????.td > tracks.out

# Just for honks...
# grep -v skyval *.td > /tmp/foo0.cs0
# awk '{printf "%4d %10.3f %10.6f\n", $1, ($3-56482)*86400, $6-(($3-56482)*86400-15000)*15.0411/3600.0}' tracks.out > /tmp/foo
# awk '{printf "%4d %4d %9.6f %7.2f %7.2f %8.4f %8.4f %8.4f %8.4f %7.2f %6.2f %6.2f %6.2f %6.1f\n", $1, $2, ($3-56482)*24, $4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14}' tracks.out > /tmp/foo.jt1
