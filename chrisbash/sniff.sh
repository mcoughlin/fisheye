#!/bin/bash
# Script to dump out header and other quantities from a Canon image
# Syntax: sniff.sh filename

eval $@

if [ $sniffile == "HDR" ] ; then
  echo "    Image          MJD        Sky      RMS     Max     Noise    Bias"
  exit 0
fi

bias=`fitshdr $sniffile | grep "BIAS    =" | awk '{print $3}'`
noise=`fitshdr $sniffile | grep "NOISE   =" | awk '{print $3}'`
mjd=`fitshdr $sniffile | grep "MJD-OBS =" | awk '{printf "%.6f", $3}'`

sniff=`monsta $PRODIR/biasky.pro $sniffile`

printf "%s  %12.6f %8.1f %7.1f %8.1f %7d %7.1f %8.1f\n" `basename $sniffile .fits` $mjd $sniff $noise $bias
