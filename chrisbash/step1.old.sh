#!/bin/bash
# 1. Copy data to disk
# 2. Convert CR2 to FITS (cr2fits)
# 3. Sniff all the FITS files to get a sense of what's what (sniff.sh)
# 4. Make a movie
# 5. Evaluate which data are worth processing.

# Assumptions:

# data live somewhere with CR2 files:
#   dirfrom
# for example dirfrom=/media/KINGSTON/Stubbs/2013_07_15

# data are copied to destination directory with the same date name, under a
# a camera-site code subdirectory:
#   dirto/CAMSITE/YYYY_MM_DD
# for example dirto=/local/jason/s13/geo

# These are the six arguments.  If you do not like these defaults,
# supply alternates on the command line in exactly this format,
# e.g. "step1.sh sitecam=ab"

# Where are the extra binaries, bash scripts and monsta scripts?
export PRODIR=~/src/geo

# What is the directory from which these files come?
dirfrom=/media/KINGSTON/Stubbs

# What is the from-root directory where these files should go?
dirto=/local/jason/s13/geo

# What is the site-camera code?
sitecam=aa

# What filter are we extracting into a FITS file?  B, G, R, M (onochrome)
filter=M

# Create movie and jpegs if not present?
dojpg=0

# Reported time minus UTC [hr]
zone=0

eval $@

export PATH=$PRODIR:$PATH

# Load OS specific functions
source stepOS.sh

echo Copying site/camera $sitecam data from $dirfrom

if [ ! -e $dirto/$sitecam ] ; then 
  echo Please create directory $dirto/$sitecam
  echo and create there the site-cam parameter file $sitecam.sitecam
  exit 0
fi

if [ ! -e $dirto/$sitecam/$sitecam.sitecam ] ; then 
  echo Please create in $dirto/$sitecam the site-cam parameter file $sitecam.sitecam
  exit 0
fi

# Sniff out the night from one of the CR2 files
cr2=`ls $dirfrom/*.CR2 | head -1`
if [ "$cr2" == "" ] ; then cr2=`ls $dirfrom/*.cr2 | head -1` ; fi
timestamp=`mdcraw -i -v $cr2 | grep Timestamp | sed "s/Timestamp://g"`
# Assume Canon timestamp of the form "DayofWeek Month day HH:MM:SS year"
yymmdd=`datetoymd "$timestamp" $zone`

echo Copying data...
if [ ! -e $dirto/$sitecam/$yymmdd ] ; then mkdir $dirto/$sitecam/$yymmdd ; fi
cp -pvf $dirfrom/*.{cr2,CR2,JPG} $dirto/$sitecam/$yymmdd

#
# assume canon format will have prefix_nnnn.CR2
cd $dirto/$sitecam/$yymmdd

im1=`ls *.CR2 | head -1`
prein=`echo $im1 | tr "_" ' ' | awk '{print $1}'`
prefix=$sitecam$yymmdd

echo Renaming all file prefixes from $prein to $prefix
for f in *.CR2 *.cr2 ; do
  nn=`echo $f | tr '_' ' ' | awk '{print $2}'`
  mv -v $f ${prefix}_$nn
done
for f in *.JPG *.jpg ; do
  nn=`echo $f | tr '_' ' ' | awk '{print $2}'`
  mv -v $f ${prefix}_$nn
done

echo Creating the observation parameter file $prefix.obs
cp ../$sitecam.sitecam ./$prefix.obs
source ./$prefix.obs
echo "" >> $prefix.obs
echo "# Step1:" >> $prefix.obs
echo "unixdatestep1=`date +%s`       # Date step1 was run" >> $prefix.obs
echo "datestep1=\"`date`\"       # Date step1 was run" >> $prefix.obs
echo "dirfrom=$dirfrom       # source directory of CR2" >> $prefix.obs
echo "dirto=$dirto         # root directory of results" >> $prefix.obs
echo "sitecam=$sitecam         # Sitecam code" >> $prefix.obs
echo "camdate=$yymmdd         # UT date of observations" >> $prefix.obs
echo "filter=$filter         # Bayer filter processed" >> $prefix.obs

xctr=`echo $imsec | tr '[,:]' ' ' | awk '{printf "%d", ($3-$1+1)/2}'`
yctr=`echo $imsec | tr '[,:]' ' ' | awk '{printf "%d", ($4-$2+1)/2}'`
echo "xctr=$xctr         # Image center in x [pix]" >> $prefix.obs
echo "yctr=$yctr         # Image center in y [pix]" >> $prefix.obs

echo Separating files into CR2 and JPG
mkdir CR2 $filter JPG

chmod 644 *.CR2 *.JPG *.cr2 *.jpg
mv -f *.CR2 *.cr2 CR2
mv -f *.JPG *.jpg JPG

cd $dirto/$sitecam/$yymmdd/CR2
if [ $filter == "M" ] ; then
  color=bw
elif [ $filter == "B" ] ; then
  color=b
elif [ $filter == "G" ] ; then
  color=g
elif [ $filter == "R" ] ; then
  color=r
else
  echo Unrecognized filter $filter 
  exit 1
fi

echo "Extracting $filter FITS from CR2, canon2fits -$color ..."

cr2fits -dir ../$filter -$color -image=$imsec -bias=$biasec -zone=$zone *.CR2

cd $dirto/$sitecam/$yymmdd/$filter
echo Sniffing all images...
for f in *.fits ; do sniff.sh sniffile=$f PRODIR=$PRODIR ; done > sniff.out

if [ $dojpg -eq 1 ] ; then 
  echo Making a movie of all the images...

  cd $dirto/$sitecam/$yymmdd/JPG
  njpg=`ls *.JPG | wc | awk '{print $1}'`

  if [ $njpg -gt 10 ] ; then
    n=0
    for f in *.JPG ; do
      imnum=`echo $f | tr "_." ' ' | awk '{print $2}'`
      convert $f -scale 20\%x20\%  -gravity northeast -stroke "#FFF" -annotate 0 "$imnum" `printf "mov-%04d.jpg" $n`
      let n++
    done
    ffmpeg -f image2 -i mov-%04d.jpg -sameq -r 12 movie.mp4

  else
    cd $dirto/$sitecam/$yymmdd/CR2
    n=0
    for f in *.CR2 ; do
      imnum=`echo $f | tr "_." ' ' | awk '{print $2}'`
      mdcraw -e -c $f | convert - -scale 20\%x20\%  -gravity northeast -stroke "#FFF" -annotate 0 "$imnum" ../JPG/`printf "mov-%04d.jpg" $n`
      let n++
    done
    cd ../JPG
    ffmpeg -f image2 -i mov-%04d.jpg -sameq -r 12 movie.mp4
  fi

  if [ -e ../JPG/movie.mp4 ] ; then
    echo Movie is in $dirto/$sitecam/$yymmdd/JPG/movie.mp4
  fi
fi

echo Evaluating results from night $yymmdd from camera $sitecam

cd $dirto/$sitecam/$yymmdd/$filter

cat sniff.out | tr '_' ' ' | awk '{printf "%d %.1f %.1f\n", $2,$4,$7}' > foo.sniff

mongo $PRODIR/sniff.pl

echo Please evaluate the range of good frame numbers from sniff.ps...
echo Pass these on to step2.sh as im0=n1 im1=n2
