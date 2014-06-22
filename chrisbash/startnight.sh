#!/bin/bash
# modified Dec 12 2013 to allow N frames set from command line, overwrites what was set in initialize.sh
export maxframes=$1
echo first arguent is $1

source ~/initialize.sh

echo maxframes is $1

while [ $framecounter -lt $maxframes ]; do

echo sleeping for $pauseinterval
sleep $pauseinterval

gphoto2 --set-config=/main/capturesettings/shutterspeed=3
# this next trick pads out the framecounter to 4 digits for decent file names
export paddedcounter=`printf  "%04i\n" $framecounter`
gphoto2 --capture-image-and-download --filename $dirpath/$dirname/$dirname.$paddedcounter.short.cr2

let framecounter=framecounter+1
export paddedcounter=`printf  "%04i\n" $framecounter`
gphoto2 --set-config=/main/capturesettings/shutterspeed=30
export paddedcounter=`printf  "%04i\n" $framecounter`
gphoto2 --capture-image-and-download --filename $dirpath/$dirname/$dirname.$paddedcounter.long.cr2

let framecounter=framecounter+1
done

# at this stage, we're done collecting images for the night. Wrap things up...
ls $dirpath/$dirname/*.cr2 | wc | awk '{print $1," images collected"}' >> $dirpath/$dirname/$dirname.log

# now do image conversion to fits, and run some initial analysis.
echo making fits files at `date`
source ~/makefits.sh
echo extracting sky values at `date`
source ~/getsky.sh
echo doing photometry at `date`
source ~/photometry.sh
echo making sky brightness plot
source ~/makeskyplot.sh

# copy result files to Amazon Web Services machine
echo copying files to aws at `date`

scp -i ~/aws/aws1.pem.txt -r $dirpath/$dirname/$dirname.* ec2-user@54.200.60.175:~/data/

echo compressing and removing images at `date`
cd $dirpath/$dirname
rm B/*.fits
rm G/*.fits
rm R/*.fits
rm M/*.fits
gzip CR2/*.cr2
mkdir /Volumes/2TB/$dirname
mv CR2/*.gz /Volumes/2TB/$dirname/

df -h

# moving the log file to appropriate directory at `date`

mv ~/cronlog $dirpath/$dirname/$dirname.cronlog

