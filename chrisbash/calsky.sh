#!/bin/bash
# Pick apart an .ics file from CalSky
# Syntax calsky.sh Calsky-Event.ics

echo " Long    mag    Azi    Alt    Dist     RA       Dec     Sat         Sat"

grep DESCRIPTION $1 | sed 's/\\n/\n/g' | egrep 'DESCRIPT|Azi|Orb|RA' | tr '\n' ' ' | sed 's/DESCRIPTION:/\n/g' | tr -d '\302' | tr '\260' d | tr \' m | sed 's/= /=/g' | awk 'NR>1{print $0}' > /tmp/foo
exec 10</tmp/foo
while read LINE<&10 ; do 
  Satellite=`echo $LINE | sed 's/Satellite //g' | awk -F\( '{print $1}'`
  Position="0.0" Magnitude="0.0"
  west=1
  for a in $LINE ; do
    e=`echo $a | grep =`
    if [ "$e" != "" ] ; then eval $a ; fi
    if [ $a == "West" ] ; then west=-1; fi
    id=`echo $a | grep '[1-2]...-...-[A-Z])'`
    if [ "$id" != "" ] ; then desig=`echo $id | tr -d \)` ; fi
  done
  West=`echo $LINE | grep West`
  long=`echo $Position | tr -d d | awk -v w=$west '{print $1*w}'`
  mag=`echo $Magnitude | sed 's/mag//g'`
  az=`echo $Azimuth | tr -d d`
  alt=`echo $Altitude | tr -d d`
  dist=`echo $Distance | sed 's/km//g'`
  ra=`echo $RA | tr "[a-z]" ' ' | awk '{printf "%.3f", 15*($1+$2/60.0)}'`
  neg=`echo $Dec | grep "-"`
  if [ "$neg" != "" ] ; then south=-1 ; else south=1 ; fi
  dec=`echo $Dec | tr "[a-z-]" ' ' | awk -v s=$south '{printf "%.3f", s*($1+$2/60.0)}'`

  printf "%6.1f %6.1f %6.1f %6.1f %8.1f %7.3f %8.3f  %s  \"%s\"\n" $long $mag $az $alt $dist $ra $dec $desig "$Satellite"
done

rm /tmp/foo
