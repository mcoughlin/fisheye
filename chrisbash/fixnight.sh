#!/bin/bash
# modified March 1 2014 to allow us to process data taken but not reduced. 
# Takes one command line argument, namely directory to fix, e.g. ut022514
# C. Stubbs 

export fixdirpath=/Users/christopherstubbs/data
export fixdirname=$1

# at this stage, we're done collecting images for the night. Wrap things up...
ls $fixdirpath/$fixdirname/*.cr2 | wc | awk '{print $1," images collected"}' >> $fixdirpath/$fixdirname/$fixdirname.log

# now do image conversion to fits, and run some initial analysis.
echo making fits files at `date`
source ~/fixmakefits.sh
echo extracting sky values at `date`
source ~/fixgetsky.sh
cd $fixdirpath/$fixdirname
source ~/ticker.sh & 
echo doing photometry at `date`
source ~/fixphotometry.sh
echo making sky brightness plot
source ~/fixmakeskyplot.sh

# copy result files to Amazon Web Services machine
echo copying files to aws at `date`

scp -i ~/aws/aws1.pem.txt -r $fixdirpath/$fixdirname/$fixdirname.* ec2-user@54.200.60.175:~/data/

echo compressing and removing images at `date`
cd $fixdirpath/$fixdirname
rm B/*.fits
rm G/*.fits
rm R/*.fits
rm M/*.fits
gzip CR2/*.cr2
mkdir /Volumes/2TB/$fixdirname
mv CR2/*.gz /Volumes/2TB/$fixdirname/

df -h

# moving the log file to appropriate directory at `date`

mv ~/cronlog $fixdirpath/$fixdirname/$fixdirname.cronlog

