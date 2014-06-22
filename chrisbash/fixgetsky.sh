#!/bin/bash

# for a rough measure of sky brightness, extract mean for a region near the center of sensor, in each band. Uses header to subtract bias scalar

cd $fixdirpath/$fixdirname/M
rm *.sky.*.dat
rm *.bias.*.dat
rm *debias*.dat
for i in *.M.fits; do getpix $i 800-1200 800-1200 -m | grep Mean | awk '{print $2}' >> $fixdirname.sky.M.dat ; done
gethead *.fits BIAS EXPTIME MJD-OBS  > $fixdirname.M.obslog 
paste $fixdirname.M.obslog $fixdirname.sky.M.dat >> temp
awk '{print ($5-$2)/$3}' temp >> $fixdirname.skydebiased.M.dat
# put sky rate ADU/sec/pix into obslog file
paste $fixdirname.M.obslog $fixdirname.skydebiased.M.dat >> temp2
mv temp2 $fixdirname.M.obslog
rm temp2
rm temp

cd $fixdirpath/$fixdirname/B
rm *.sky.*.dat
rm *.bias.*.dat
rm *debias*.dat
for i in *.B.fits; do getpix $i 800-1200 800-1200 -m | grep Mean | awk '{print $2}' >> $fixdirname.sky.B.dat ; done
gethead *.fits BIAS EXPTIME MJD-OBS  > $fixdirname.B.obslog 
paste $fixdirname.B.obslog $fixdirname.sky.B.dat >> temp
awk '{print ($5-$2)/$3}' temp >> $fixdirname.skydebiased.B.dat
# put sky rate ADU/sec/pix into obslog file
paste $fixdirname.B.obslog $fixdirname.skydebiased.B.dat >> temp2
mv temp2 $fixdirname.B.obslog
rm temp2
rm temp

cd $fixdirpath/$fixdirname/G
rm *.sky.*.dat
rm *.bias.*.dat
rm *debias*.dat
for i in *.G.fits; do getpix $i 800-1200 800-1200 -m | grep Mean | awk '{print $2}' >> $fixdirname.sky.G.dat ; done
gethead *.fits BIAS EXPTIME MJD-OBS  > $fixdirname.G.obslog 
paste $fixdirname.G.obslog $fixdirname.sky.G.dat >> temp
awk '{print ($5-$2)/$3}' temp >> $fixdirname.skydebiased.G.dat
# put sky rate ADU/sec/pix into obslog file
paste $fixdirname.G.obslog $fixdirname.skydebiased.G.dat >> temp2
mv temp2 $fixdirname.G.obslog
rm temp2
rm temp

cd $fixdirpath/$fixdirname/R
rm *.sky.*.dat
rm *.bias.*.dat
rm *debias*.dat
for i in *.R.fits; do getpix $i 800-1200 800-1200 -m | grep Mean | awk '{print $2}' >> $fixdirname.sky.R.dat ; done
gethead *.fits BIAS EXPTIME MJD-OBS  > $fixdirname.R.obslog 
paste $fixdirname.R.obslog $fixdirname.sky.R.dat >> temp
awk '{print ($5-$2)/$3}' temp >> $fixdirname.skydebiased.R.dat
# put sky rate ADU/sec/pix into obslog file
paste $fixdirname.R.obslog $fixdirname.skydebiased.R.dat >> temp2
mv temp2 $fixdirname.R.obslog
rm temp2
rm temp

cd $fixdirpath/$fixdirname/R
ls *.R.fits > listing
mv listing ..
cd ..
              
echo "#image     BIAS EXPTIME MJD-OBS                       M        B       G       R nstars " > $fixdirname.obslog
paste M/$fixdirname.M.obslog  B/*.skydebiased.B.dat G/*.skydebiased.G.dat R/*.skydebiased.R.dat >> $fixdirname.obslog
