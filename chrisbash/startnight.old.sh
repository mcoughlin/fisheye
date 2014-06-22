#!/bin/bash
source ~/initialize.sh

while [ $framecounter -lt $maxframes ]; do

echo sleeping for $pauseinterval
sleep $pauseinterval

gphoto2 --set-config=/main/capturesettings/shutterspeed=1
# this next trick pads out the framecounter to 4 digits for decent file names
export paddedcounter=`printf  "%04i\n" $framecounter`
gphoto2 --capture-image-and-download --filename $dirpath/$dirname/$dirname.$paddedcounter.short.cr2

let framecounter=framecounter+1
export paddedcounter=`printf  "%04i\n" $framecounter`
gphoto2 --set-config=/main/capturesettings/shutterspeed=10
export paddedcounter=`printf  "%04i\n" $framecounter`
gphoto2 --capture-image-and-download --filename $dirpath/$dirname/$dirname.$paddedcounter.long.cr2

let framecounter=framecounter+1
done

# at this stage, we're done collecting images for the night. Wrap things up...
ls $dirpath/$dirname/*.cr2 | wc | awk '{print $1," images collected"}' >> $dirpath/$dirname/$dirname.log

# now do image conversion to fits, and run some initial analysis.

source ~/makefits.sh
source ~/getsky.sh
source ~/photometry.sh

# stick cronlog file in that night's directory
mv ~/cronlog $dirpath/$dirname/$dirname.cronlog

# copy result files to Amazon Web Services machine
scp -i ~/aws/aws1.pem.txt -r $dirpath/$dirname/$dirname.* ec2-user@54.200.60.175:~/data/

