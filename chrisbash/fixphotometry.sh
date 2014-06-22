#!/bin/bash

# performs quick photometry on each image, using tphot
# modified Jan 13 2014 to remove bias from each image to get appropriate sky levels. 

cd $fixdirpath/$fixdirname/M
rm *.phot
rm *.nstars
for i in *.M.fits; do imagebias=`gethead $i BIAS`; tphot $i -out `basename $i .fits`.phot -bias $imagebias; done
# how many stars in each image? 
for i in *.M.phot; do wc $i  | awk '{print ($1-1)}'  >> $fixdirname.nstars ; done
mv $fixdirname.nstars ..

cd $fixdirpath/$fixdirname/B
rm *.phot
for i in *.B.fits; do imagebias=`gethead $i BIAS`; tphot $i -out `basename $i .fits`.phot; done

cd $fixdirpath/$fixdirname/G
rm *.phot
for i in *.G.fits; do imagebias=`gethead $i BIAS`; tphot $i -out `basename $i .fits`.phot; done

cd $fixdirpath/$fixdirname/R
rm *.phot
for i in *.R.fits; do imagebias=`gethead $i BIAS`; tphot $i -out `basename $i .fits`.phot; done

cd ..
rm temp
# add extra line to nstars file to accommodate header line coming up
echo " " > temp2
cat $fixdirname.nstars >> temp2
mv temp2 $fixdirname.nstars
rm temp2
# add another column to the obslog file
paste $fixdirname.obslog $fixdirname.nstars >> temp
mv temp $fixdirname.obslog

# separate out long and short exposures
echo "#image     BIAS EXPTIME MJD-OBS                       M        B       G       R nstars" > $fixdirname.short.obslog
grep short $fixdirname.obslog >> $fixdirname.short.obslog
echo "#image     BIAS EXPTIME MJD-OBS                       M        B       G       R nstars" > $fixdirname.long.obslog
grep long $fixdirname.obslog >> $fixdirname.long.obslog
 

