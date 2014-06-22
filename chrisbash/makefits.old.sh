#!/bin/bash
# Convert CR2 to FITS (cr2fits)

# Assumptions:

# Where are the extra binaries, bash scripts and monsta scripts?
export PRODIR=/usr/local/bin

# What is the from-root directory where these files should go?
dirto=$dirpath/$dirname
dirfrom=$dirto

# What filter are we extracting into a FITS file?  B, G, R, M (onochrome)
filter=M

eval $@

export PATH=$PRODIR:$PATH

# Load OS specific functions
source stepOS.sh

yymmdd=`echo $dirname | sed s/ut//`

#
# assume canon format will have prefix.nnnn.imtype.cr2
cd $dirto
prefix=$dirname

echo " " >> $prefix.log
echo "# Step1:" >> $prefix.log
echo "datestep1=\"`date`\"       # Date step1 was run" >> $prefix.log
echo "dirfrom=$dirfrom       # source directory of CR2" >> $prefix.log
echo "dirto=$dirto         # root directory of results" >> $prefix.log
echo "obsdate=$dirname        # UT date of observations" >> $prefix.log

echo Separating files into CR2 and color subdirectories
mkdir CR2 M B G R 
rm M/*.fits
rm B/*.fits
rm G/*.fits
rm R/*.fits

chmod 644 *.CR2 *.cr2 
mv -f *.CR2 *.cr2 CR2

cd $dirto/CR2

filter=M 
color=bw
echo "Extracting $filter FITS from CR2, raw2fits -$color ..."
cr2fits -dir ../$filter -$color  *.cr2

filter=B 
color=b
echo "Extracting $filter FITS from CR2, raw2fits -$color ..."
cr2fits -dir ../$filter -$color  *.cr2

filter=G 
color=g
echo "Extracting $filter FITS from CR2, raw2fits -$color ..."
cr2fits -dir ../$filter -$color  *.cr2

filter=R 
color=r
echo "Extracting $filter FITS from CR2, raw2fits -$color ..."
cr2fits -dir ../$filter -$color  *.cr2

echo "fixing up fits file names"
cd $dirto/M
for i in *.fits; do mv "$i" "${i/.fits}".M.fits; done

cd $dirto/B
for i in *.fits; do mv "$i" "${i/.fits}".B.fits; done

cd $dirto/G
for i in *.fits; do mv "$i" "${i/.fits}".G.fits; done

cd $dirto/R
for i in *.fits; do mv "$i" "${i/.fits}".R.fits; done


