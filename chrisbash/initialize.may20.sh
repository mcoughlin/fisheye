#!/bin/bash
# try to fix up user issues for crontab
source /Users/christopherstubbs/.bash_profile

# set up directory names and data paths
export datadisk=/dev/disk1s2
export dirpath=/Users/christopherstubbs/data
export dirname=`date +"ut%m%d%y"`
export thismonth=`date +"%m"`
# create tonight's directory
mkdir $dirpath/$dirname
cd $dirpath/$dirname
date | awk '{print "started at "$0 }'>> $dirname.log

df -H | grep $datadisk | sed s/G//g | awk '{print $4," GB left on data disk"}' >> $dirname.log
export spaceleft=`df -H | grep $datadisk | sed s/G//g | awk '{print $4}'`

# fix this later
#if ($spaceleft>20); then
#	mail -s "allsky camera lots space!" stubbs@physics.harvard.edu
#fi

# wipe out any existing connections to camera
killall PTPCamera

# set up environment variables for camera
export framecounter=0001
export xctr=2875
export yctr=1920

# prefix for images, gets framecounter.cr2 appended
export imageprefix=$dirname.
# interval after end of last image before starting next one, in seconds
export pauseinterval=50 

# configure aspects of camera that won't change
# synch camera datetime to ntp-served acquisition computer value, in UT
gphoto2 --set-config syncdatetime=1
gphoto2 --set-config iso=1600
# edited May 20 2014 from aperture of 2.8 for longer focal length lens to aperture=4.0
gphoto2 --set-config aperture=4

cd ~

