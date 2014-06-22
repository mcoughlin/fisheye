#!/bin/bash

# for a rough measure of sky brightness, extract mean for a region near the center of sensor, in each band. Uses header to subtract bias scalar

cd $dirpath/$dirname/M
rm *.sky.*.dat
rm *.bias.*.dat
rm *debias*.dat
for i in *.M.fits; do getpix $i 800-1200 800-1200 -m | grep Mean | awk '{print $2}' >> $dirname.sky.M.dat ; done
gethead *.fits BIAS EXPTIME MJD-OBS  > $dirname.M.obslog 
paste $dirname.M.obslog $dirname.sky.M.dat >> temp
awk '{print ($5-$2)/$3}' temp >> $dirname.skydebiased.M.dat
# put sky rate ADU/sec/pix into obslog file
paste $dirname.M.obslog $dirname.skydebiased.M.dat >> temp2
mv temp2 $dirname.M.obslog
rm temp2
rm temp

cd $dirpath/$dirname/B
rm *.sky.*.dat
rm *.bias.*.dat
rm *debias*.dat
for i in *.B.fits; do getpix $i 800-1200 800-1200 -m | grep Mean | awk '{print $2}' >> $dirname.sky.B.dat ; done
gethead *.fits BIAS EXPTIME MJD-OBS  > $dirname.B.obslog 
paste $dirname.B.obslog $dirname.sky.B.dat >> temp
awk '{print ($5-$2)/$3}' temp >> $dirname.skydebiased.B.dat
# put sky rate ADU/sec/pix into obslog file
paste $dirname.B.obslog $dirname.skydebiased.B.dat >> temp2
mv temp2 $dirname.B.obslog
rm temp2
rm temp

cd $dirpath/$dirname/G
rm *.sky.*.dat
rm *.bias.*.dat
rm *debias*.dat
for i in *.G.fits; do getpix $i 800-1200 800-1200 -m | grep Mean | awk '{print $2}' >> $dirname.sky.G.dat ; done
gethead *.fits BIAS EXPTIME MJD-OBS  > $dirname.G.obslog 
paste $dirname.G.obslog $dirname.sky.G.dat >> temp
awk '{print ($5-$2)/$3}' temp >> $dirname.skydebiased.G.dat
# put sky rate ADU/sec/pix into obslog file
paste $dirname.G.obslog $dirname.skydebiased.G.dat >> temp2
mv temp2 $dirname.G.obslog
rm temp2
rm temp

cd $dirpath/$dirname/R
rm *.sky.*.dat
rm *.bias.*.dat
rm *debias*.dat
for i in *.R.fits; do getpix $i 800-1200 800-1200 -m | grep Mean | awk '{print $2}' >> $dirname.sky.R.dat ; done
gethead *.fits BIAS EXPTIME MJD-OBS  > $dirname.R.obslog 
paste $dirname.R.obslog $dirname.sky.R.dat >> temp
awk '{print ($5-$2)/$3}' temp >> $dirname.skydebiased.R.dat
# put sky rate ADU/sec/pix into obslog file
paste $dirname.R.obslog $dirname.skydebiased.R.dat >> temp2
mv temp2 $dirname.R.obslog
rm temp2
rm temp

cd $dirpath/$dirname/R
ls *.R.fits > listing
mv listing ..
cd ..
              
echo "#image     BIAS EXPTIME MJD-OBS                       M        B       G       R nstars " > $dirname.obslog
paste M/$dirname.M.obslog  B/*.skydebiased.B.dat G/*.skydebiased.G.dat R/*.skydebiased.R.dat >> $dirname.obslog
