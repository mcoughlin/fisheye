#!/bin/bash
# Process a fisheye observation
# Syntax: fishy obs.fits [args]
# E.g.    fishy 56713.2500231.fits [args]
#
# NOTE: fisheye and fishiter.sh must be in the path
#
# Writes $obs.fish and $obs.mch
#
# Examine results with mongo fishy.pl ; astrometry $obs ; photometry $obs
#
# 140413 v1.1: uses Hipparcos, can bootstrap off of x,y,RA,Dec

# First argument is always the FITS file
fits=$1
shift 1
obs=`basename $fits .fits`

# Observatory longitude and latitude (REQUIRED)
LNG=-155.5763           # MLO
LAT=19.5362             # MLO

# Bootstrap file of x,y,RA,Dec (REQUIRED if pointing not provided)
BOOT=""

# Initial pointing parameters (REQUIRED if bootstrap file not provided)
#
## Note that these have to be quite accurate, the initial match of
## MCHTOL1, default 30 pix, will fail unless HA, DEC, and AZ are 
## accurate to dANG < 1/2 MCHTOL1*SCALE (i.e. 1 deg or so) and SCALE is
## accurate to dSCALE/SCALE < 1/3 MCHTOL1/1000pix (i.e. 1% or better)
#
HA=""
DEC=""
AZ=""
SCALE=""

# Where can hip.dat.dat be found: $CATDIR/hip.dat (REQUIRED)
CATDIR=..

# Other default parameters
VERB=0

# Minimum number of tphot detections to proceed
MINSTAR=100

# Match tolerance [pix]
MCHTOL1=30           # First match tolerance [pix] at m<3 and N<100
MCHTOL=5             # Subsequence match tolerances [pix]
RSIG=3               # Final pruning sigma

# Clean up afterwards?
CLEAN=1

eval $@

# Bootstrap?
if [[ ! -z $BOOT ]] ; then
  nboot=`wc $BOOT | awk '{print $1}'`
  if [[ -z $nboot || $nboot -lt 3 ]] ; then
    echo Too few stars in boot file N=$nboot, aborting
    exit 1
  fi
  utcobs=`fitshdr $obs.fits | grep "MJD-OBS =" | awk '{print $3}'`
  if [[ -z $utcobs ]] ; then
    echo MJD-OBS not found in FITS header, aborting
    exit 1
  fi
  gmst=`echo $utcobs | awk '{gmst=(280.460618+360.985647366*($1-51544.5))%360.0; gmst=gmst<0?gmst+360:gmst; printf "%.6f", gmst}'`
  LST=`echo $gmst $LNG | awk '{lst=$1+$2; lst=lst>0?lst:lst+360; printf "%.6f\n", lst}'`

# Convert RA to HA
  awk -v lst=$LST '{print $1,$2,lst-$3,$4,$5}' $BOOT > $BOOT.ha

# Show us the results?
  if [[ $VERB -gt 0 ]] ; then echo "UTC= $utcobs LST= $LST triangles: `fisheye $BOOT.ha`" ; fi

# Capture this triangle estimate of pole, az, and scale,
  eval `fisheye $BOOT.ha`
# and fit the bootstrap data
  if [[ $VERB -gt 0 ]] ; then
    echo Bootstrap fit:
  fi
  fisheye $BOOT.ha -ha $HA -dec $DEC -az $AZ -scale $SCALE -lat $LAT -fixdist -fixctr -out $obs.fishb

# Load up the pole, az, and scale
  boot=(`grep "HA=" $obs.fishb`)
  HA=${boot[2]} DEC=${boot[4]} AZ=${boot[6]} SCALE=${boot[8]}

  if [[ $CLEAN -gt 0 ]] ; then rm $BOOT.ha ; fi

else
  if [[ -z $HA || -z $DEC || -z $AZ || -z $SCALE ]] ; then
    echo "Error: must specify either BOOT= or else HA= DEC= AZ= SCALE="
    exit 1
  fi
fi

if [[ $VERB -gt 0 ]] ; then 
  echo Initial params: HA= $HA DEC= $DEC AZ= $AZ SCALE= $SCALE LNG= $LNG LAT= $LAT
fi

# tphot the image, lots o' hardwired parameters that need tuning...
TRAD=4 TSIG=2 TMIN=50 TBIAS=2048 TAPRAD=4 

tphot $obs.fits -bias $TBIAS -rad $TRAD -sig $TSIG -min $TMIN -aprad $TAPRAD -okfit 0 -chin 100000 -out $obs.tph

nstar=`wc $obs.tph | awk '{print $1}'`
if [[ $nstar -lt $MINSTAR ]] ; then
  echo Error: $obs $nstar stars found, too few - aborting...
  exit 1
else
  echo Observation $obs tphot found $nstar stars
fi

# First iteration: MLIM=3, NSTAR=100, matched fit then pruned fit:
if [[ $VERB -gt 0 ]] ; then
  echo First iteration: MLIM=3, NSTAR=100
fi
fishiter.sh $obs CATDIR=.. MLIM=3 NSTAR=100 MCHTOL=$MCHTOL1 HA=$HA DEC=$DEC AZ=$AZ SCALE=$SCALE LNG=$LNG LAT=$LAT CLEAN=$CLEAN VERB=$VERB

# Second iteration: MLIM=5, NSTAR=1000, matched fit then pruned fit:
if [[ $VERB -gt 0 ]] ; then
  echo Second iteration: MLIM=5, NSTAR=1000
fi
parm=(`grep "HA=" $obs.fish`)
HA=${parm[2]} DEC=${parm[4]} AZ=${parm[6]} SCALE=${parm[8]} QUAD=${parm[10]} CUBE=${parm[12]} CX=${parm[14]} CY=${parm[16]}
fishiter.sh $obs CATDIR=.. MLIM=5 NSTAR=1000 MCHTOL=$MCHTOL HA=$HA DEC=$DEC AZ=$AZ SCALE=$SCALE CX=$CX CY=$CY QUAD=$QUAD CUBE=$CUBE LNG=$LNG LAT=$LAT CLEAN=$CLEAN VERB=$VERB

# Third iteration: MLIM=6, NSTAR=2000, matched fit then pruned fit:
if [[ $VERB -gt 0 ]] ; then
  echo Third iteration: MLIM=6, NSTAR=2000
fi
parm=(`grep "HA=" $obs.fish`)
HA=${parm[2]} DEC=${parm[4]} AZ=${parm[6]} SCALE=${parm[8]} QUAD=${parm[10]} CUBE=${parm[12]} CX=${parm[14]} CY=${parm[16]}
fishiter.sh $obs CATDIR=.. MLIM=6 NSTAR=2000 MCHTOL=5 HA=$HA DEC=$DEC AZ=$AZ SCALE=$SCALE CX=$CX CY=$CY QUAD=$QUAD CUBE=$CUBE RSIG=$RSIG LNG=$LNG LAT=$LAT CLEAN=$CLEAN VERB=$VERB
