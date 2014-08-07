#!/bin/bash
# Do one fisheye fit iteration of tph file against Hipparcos
# Syntax: fishiter.sh obs HA=HA DEC=dec AZ=az SCALE=sc [arguments]
# e.g.    fishiter.sh IMG_9082 
# Expects $obs.fits, $obs.tph, Hipparcos catalog in CATDIR
# Produces $obs.fish and $obs.mch

obs=$1
shift

#LNG=-155.576300       # Site longitude
#LAT=19.5362           # Site latitude
#LST=""                # Local sidereal time

LAT=-30.67833
LNG=-70.73669
#LST=159.75

#CATDIR=..             # Directory where star catalog is found
CATDIR=/lsst/home/coughlin/git-repo/fisheye
PRODIR=/lsst/home/coughlin/git-repo/fisheye

MLIM=3                # Limiting magnitude from star catalog
NSTAR=100             # Number of stars to collect from tph file
MCHTOL=14             # Match tolerance [pix]
RSIG=5                # Prune residuals worse than RSIG * median
#MCHTOL=100

FRAT=0.0              # Vignetting correction f-ratio (0 means no correction)

#HA=0.731                 # Observation hour angle [required]
#DEC=-27.825                # Observation declination [required]
# MWC: Azimuth is for north to the left, so we are at 90
#AZ=93.550                 # Observation azimuth [required]
#SCALE=0.04594
#PARITY=-1
#MED=4.552

#HA=0.015
#DEC=-27.092
#AZ=95.082
#SCALE=0.04705
#PARITY=1
#MED=1.054

echo $HA
echo "HA = " $HA
echo "DEC = " $DEC
echo "AZ = " $AZ
echo "SCALE = " $SCALE
echo "PARITY = " $PARITY
echo "MED = " $MED

#CX=1420               # Observation image center [default, can be overridden]
#CY=975                # Observation image center [default]
#QUAD=-0.02            # Observation quadratic scale term [default]
#CUBE=-0.03            # Observation cubic scale term [default]

# After running fishiter.sh once
#QUAD=0.00870
#CUBE=-0.04598
#CX=1434.1
#CY=969.9

FIXDIST=0             # Fix distortion?
FIXCTR=0              # Fix center?

VFRAC=0.9             # m ~ 0.9*V + 0.1*B    (Canon monochrome mag)

CLEAN=2               # Clean up all the cruft?

VERB=0                # Verbosity level

eval $@

FISHEXTRA=""
if [[ $FIXCTR == 1 ]] ; then FISHEXTRA="$FISHEXTRA -fixctr" ; fi
if [[ $FIXDIST == 1 ]] ; then FISHEXTRA="$FISHEXTRA -fixdist" ; fi

# If local sidereal time not specified, get it from the UTC in the header
if [[ -z $LST ]] ; then
  utcobs=`fitshdr $obs.fits | grep "MJD-OBS =" | awk '{print $3}'`
  gmst=`echo $utcobs | awk '{gmst=(280.460618+360.985647366*($1-51544.5))%360.0; gmst=gmst<0?gmst+360:gmst; printf "%.6f", gmst}'`
  LST=`echo $gmst $LNG | awk '{lst=$1+$2; lst=lst>0?lst:lst+360; printf "%.6f\n", lst}'`
fi

echo "I AM HERE:" $LST

# Convert Hipparcos RA to HA at this LST
if [[ -e $CATDIR/hip.dat ]] ; then
  awk -v lst=$LST -v mlim=$MLIM 'NR>1 && $10<mlim{ha=lst-$2; ha=ha<=180?ha:ha-360; ha=ha>=-180?ha:ha+360; printf "%9.4f %9.4f %9.4f %6.2f %5.2f %5.2f %5.2f %5.2f\n", ha,$3,$2,$6,$8,$10,$11,$13}' $CATDIR/hip.dat > $obs.hhip
elif [[ -e $CATDIR/hip.dat.bz2 ]] ; then
  bunzip2 -c $CATDIR/hip.dat.bz2 | awk -v lst=$LST -v mlim=$MLIM 'NR>1 && $10<mlim{ha=lst-$2; ha=ha<=180?ha:ha-360; ha=ha>=-180?ha:ha+360; printf "%9.4f %9.4f %9.4f %6.2f %5.2f %5.2f %5.2f %5.2f\n", ha,$3,$2,$6,$8,$10,$11,$13}' > $obs.hhip
else
  echo "Error: fishiter cannot find '$CATDIR/hip.dat or hip.dat.bz2'"
  exit 1
fi

# Convert these to first cut at xy coordinates for this obs
fisheye -ha $HA -dec $DEC -az $AZ -scale $SCALE -quad $QUAD -cube $CUBE -cx $CX -cy $CY -lat $LAT -out - -sky2xy < $obs.hhip | awk '{printf "%7.2f %7.2f\n", $1,$2}' > $obs.xyhip

# Write a Hipparcos file for this obs: $obs.hip
echo "   RA         Dec      HA      xcalc   ycalc     BT    VT    V   (B-V) (V-I)" > $obs.hip
colmerge 0 $obs.xyhip 0 $obs.hhip -nobar | awk '{printf "%9.4f %9.4f %9.4f %7.2f %7.2f %6.2f %5.2f %5.2f %5.2f %5.2f\n", $5,$4,$3,$1,$2,$6,$7,$8,$9,$10}' >> $obs.hip

# Match up these predictions with the brightest NSTAR stars in the tphot file
sort -g -r -k 8 $obs.tph | head -$NSTAR | colmerge 4,5 $obs.hip 1,2 - -tol $MCHTOL > $obs.mch0

# How many matches?
nmch=`wc $obs.mch0 | awk '{print $1}'`
if [[ $nmch -lt 3 ]] ; then
  echo WHOA: only $nmch stars matched, something is wrong
  exit 0
fi

# Extract just x,y HA,Dec: $obs.xyrd
awk '{printf "%7.2f %7.2f %9.4f %9.4f\n",$12,$13,$3,$2}' $obs.mch0 > $obs.xyrd0

if [[ $VERB -gt 1 ]] ; then
  echo Initial params: $HA $DEC $AZ $SCALE $QUAD $CUBE $CX $CY
fi

# Fit this new match up: $obs.fish0
fisheye -ha $HA -dec $DEC -az $AZ -scale $SCALE -quad $QUAD -cube $CUBE -cx $CX -cy $CY -lat $LAT -in $obs.xyrd0 -out $obs.fish0 $FISHEXTRA

# Prune bad points update x,y, HA,Dec: $obs.xyrd1
RMED=`awk 'NR>2{print $17}' $obs.fish0 | median verb | awk '{print $1}'`
awk -v rmed=$RMED -v rsig=$RSIG 'NR>3 && $17<rmed*rsig{printf "%7.2f %7.2f %9.4f %9.4f\n",$1,$2,$3,$4}' $obs.fish0 > $obs.xyrd1

# Refit cleaned list
fisheye -ha $HA -dec $DEC -az $AZ -scale $SCALE -quad $QUAD -cube $CUBE -cx $CX -cy $CY -lat $LAT -in $obs.xyrd1 -out $obs.fish $FISHEXTRA

# Get parameters of revised fit
parm=(`grep "HA=" $obs.fish`)
HA=${parm[2]} DEC=${parm[4]} AZ=${parm[6]} SCALE=${parm[8]}
QUAD=${parm[10]} CUBE=${parm[12]} CX=${parm[14]} CY=${parm[16]}

if [[ $VERB -gt 1 ]] ; then
  echo Revised params: $HA $DEC $AZ $SCALE $QUAD $CUBE $CX $CY
fi

# Dump out all we care from the tphot output: x,y,m,dm,sky
awk 'NR>1 && $8>0{printf "%7.2f %7.2f %6.2f %5.2f %5.1f\n",$1,$2,-2.5*log($8)/log(10),$9/$8,$4}' $obs.tph > $obs.minst

# Match astrometry with the photometry: $obs.mch
printf "    RA        Dec      HA      xcalc  ycalc    BT     VT     V    (B-V)  (V-I)     m"  > $obs.mch
printf "      x       y    HAcalc Decalc  alt    rad   minst   dm   sky  dvig\n" >> $obs.mch
colmerge 3,2 $obs.hip 3,4 $obs.fish -tol 0.01 | colmerge 12,13 - 1,2 $obs.minst -tol 0.1  | awk -v f=$VFRAC -v frat=$FRAT '{skyflux=$34; az=$19; alt=$20; jac=$21; skysec=skyflux/jac; mu=skysec>0?-2.5/log(10)*log(skysec):0.0; mcat=f*$8+(1-f)*($8+$9); x=$30; y=$31; r=$24; minst=$32; dm=$33; dvig=frat>0?1-24*exp(-3*log(frat))*(r/1000)*(r/1000):1.0; dvig = 2.5/log(10)*log(dvig); printf "%9.4f %9.4f %9.4f %6.1f %6.1f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %7.2f %7.2f %7.2f %6.2f %5.2f %6.1f %6.2f %5.2f %5.2f %5.2f\n", $1,$2,$3,$22,$23,$6,$7,$8,$9,$10,mcat,x,y,$16,$17,alt,r,minst+dvig,dm,mu+dvig,dvig}' >> $obs.mch

if [[ $CLEAN > 0 ]] ; then
  rm $obs.hhip $obs.xyhip $obs.hip $obs.xyrd0 $obs.xyrd1 $obs.minst
  rm $obs.mch0 $obs.fish0 
fi
