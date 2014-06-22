#!/bin/bash
# modified Dec 12 2013 to allow N frames set from command line, overwrites what was set in initialize.sh
export maxframes=$1
echo first argument is $1

source ~/dayinitialize.sh

echo maxframes is $1
echo framecounter is set to $2

export pauseinterval=600
export framecounter=$2

while [ $framecounter -lt $maxframes ]; do

echo sleeping for $pauseinterval
sleep $pauseinterval

# adjust these to get appropriate sun images
gphoto2 --set-config=/main/capturesettings/shutterspeed=1/8000
gphoto2 --set-config aperture=11

# this next trick pads out the framecounter to 4 digits for decent file names
export paddedcounter=`printf  "%04i\n" $framecounter`
gphoto2 --capture-image-and-download --filename $dirpath/$dirname/$dirname.daycal.$paddedcounter.cr2

let framecounter=framecounter+1
export paddedcounter=`printf  "%04i\n" $framecounter`

done

# at this stage, we're done collecting images for the night. Wrap things up...
ls $dirpath/$dirname/*.cr2 | wc | awk '{print $1," images collected"}' >> $dirpath/$dirname/$dirname.daycal.log

# now do image conversion to fits, and run some initial analysis.

# echo making fits files at `date`
#source ~/makefits.sh

# gzip CR2/*.cr2
# mkdir /Volumes/2TB/$dirname
# mv CR2/*.gz /Volumes/2TB/$dirname/

# df -h

#moving the log file to appropriate directory at `date`

# mv ~/cronlog $dirpath/$dirname/$dirname.cronlog

