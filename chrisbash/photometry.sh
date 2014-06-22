#!/bin/bash

# performs quick photometry on each image, using tphot
# modified Jan 13 2014 to remove bias from each image to get appropriate sky levels. 

cd $dirpath/$dirname/M
rm *.phot
rm *.nstars
for i in *.M.fits; do imagebias=`gethead $i BIAS`; tphot $i -out `basename $i .fits`.phot -bias $imagebias; done
# how many stars in each image? 
for i in *.M.phot; do wc $i  | awk '{print ($1-1)}'  >> $dirname.nstars ; done
mv $dirname.nstars ..

cd $dirpath/$dirname/B
rm *.phot
for i in *.B.fits; do imagebias=`gethead $i BIAS`; tphot $i -out `basename $i .fits`.phot; done

cd $dirpath/$dirname/G
rm *.phot
for i in *.G.fits; do imagebias=`gethead $i BIAS`; tphot $i -out `basename $i .fits`.phot; done

cd $dirpath/$dirname/R
rm *.phot
for i in *.R.fits; do imagebias=`gethead $i BIAS`; tphot $i -out `basename $i .fits`.phot; done

cd ..
rm temp
# add extra line to nstars file to accommodate header line coming up
echo " " > temp2
cat $dirname.nstars >> temp2
mv temp2 $dirname.nstars
rm temp2
# add another column to the obslog file
paste $dirname.obslog $dirname.nstars >> temp
mv temp $dirname.obslog

# separate out long and short exposures
echo "#image     BIAS EXPTIME MJD-OBS                       M        B       G       R nstars" > $dirname.short.obslog
grep short $dirname.obslog >> $dirname.short.obslog
echo "#image     BIAS EXPTIME MJD-OBS                       M        B       G       R nstars" > $dirname.long.obslog
grep long $dirname.obslog >> $dirname.long.obslog
 

