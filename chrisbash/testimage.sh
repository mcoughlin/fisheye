#!/bin/bash
# modified Dec 12 2013 to allow N frames set from command line, overwrites what was set in initialize.sh

rm test.cr2
rm test.fits

# adjust these to get appropriate  images
gphoto2 --set-config=/main/capturesettings/shutterspeed=1/8000
gphoto2 --set-config aperture=11

gphoto2 --capture-image-and-download --filename test.cr2

cr2fits -bw test.cr2
getpix test.fits 800-1200 800-1200 -m | grep Mean 
echo BIAS =`gethead BIAS test.fits`
