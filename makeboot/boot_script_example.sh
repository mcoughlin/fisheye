Identify 8 stars by hand with stellarium and monsta, join their x,y to RA,Dec, write boot.ra:
{{{
1039.9  318.3  187.79141 -57.11319    1.63  1.59  M3.5III   Gam_Cru     Gacrux (top S cross)
 856.6 1040.5  188.59683 -23.39672    2.65  0.89  G5II      9Bet_Crv    Beta-crv (east)
 569.1 1259.8  201.29829 -11.16124    0.98 -0.23  B1III     67Alp_Vir   Spica
2121.7  122.2   95.98789 -52.69565   -0.72  0.15  F0II      Alp_Car     Canopus (under Orion)
2442.2 1475.2  114.82523   5.22496    0.38  0.42  F5IV-V    10Alp_CMi   Procyon
2376.5  641.9  104.65646 -28.97213    1.50 -0.21  B2II      21Eps_CMa   Adhara (r foot)
1650.2 1799.9  152.09291  11.96719    1.35 -0.11  B7V       32Alp_Leo   Regulus
 317.4 1889.6  213.91533  19.18219   -0.04  1.23  K1.5III   16Alp_Boo   Arcturus
}}}

CTIO is LAT=-30.67833, LNG=-70.73669, and the MJD is 56784.015602, so LST=159.655400:
{{{
 LAT=-30.67833 LNG=-70.73669
 utcobs=56784.015602
 gmst=`echo $utcobs | awk '{gmst=(280.460618+360.985647366*($1-51544.5))%360.0; gmst=gmst<0?gmst+360:gmst; printf "%.6f", gmst}'`
 LST=`echo $gmst $LNG | awk '{lst=$1+$2; lst=lst>0?lst:lst+360; printf "%.6f\n", lst}'`
}}}

Convert boot.ra to boot.ha and feed it to fisheye:
{{{
 awk -v lst=$LST '{printf "%7.1f %7.1f %9.4f %9.4f %9.4f  %s\n", $1+0.5,$2+0.5,lst-$3,$4,$3,$9}' boot.ra > boot.ha

 fisheye < boot.ha
HA=0.015 DEC=-27.092 AZ=95.082 SCALE=0.04705 PARITY=1 MED=1.054
}}}

Emulate the imstats_fish.sh's use of fishiter.sh:
{{{
 bunzip2 -c /local/catalogs/hip.dat.bz2 > /tmp/hip.dat

 CATDIR=/tmp

 CLEAN=0 VERB=2
 MCHTOL1=30           # First match tolerance [pix] at m<3 and N<100
 MCHTOL=5             # Subsequence match tolerances [pix]
 RSIG=3               # Final pruning sigma

 LAT=-30.67833 LNG=-70.73669
 eval $(fisheye < boot.ha)

 obs=ut050614.0200.long.M
 frat=`fitshdr $obs.fits | grep APERTURE | awk '{printf "%.6f", $2}'`
# NOT!
# frat=2.8 will blow up fishiter.sh, fishiter.sh is for the Canon 10-24mm f/4 lens
 frat=4

 fishiter.sh $obs CATDIR=$CATDIR MLIM=3 NSTAR=100 MCHTOL=$MCHTOL1 HA=$HA DEC=$DEC AZ=$AZ SCALE=$SCALE LNG=$LNG LAT=$LAT FRAT=$frat CLEAN=$CLEAN VERB=$VERB

 parm=(`grep "HA=" $obs.fish`)
 HA=${parm[2]} DEC=${parm[4]} AZ=${parm[6]} SCALE=${parm[8]} QUAD=${parm[10]} CUBE=${parm[12]} CX=${parm[14]} CY=${parm[16]}
 fishiter.sh $obs CATDIR=$CATDIR MLIM=5 NSTAR=1000 MCHTOL=$MCHTOL HA=$HA DEC=$DEC AZ=$AZ SCALE=$SCALE CX=$CX CY=$CY QUAD=$QUAD CUBE=$CUBE LNG=$LNG LAT=$LAT FRAT=$frat CLEAN=$CLEAN VERB=$VERB

 parm=(`grep "HA=" $obs.fish`)
 HA=${parm[2]} DEC=${parm[4]} AZ=${parm[6]} SCALE=${parm[8]} QUAD=${parm[10]} CUBE=${parm[12]} CX=${parm[14]} CY=${parm[16]}
 fishiter.sh $obs CATDIR=$CATDIR MLIM=6 NSTAR=2000 MCHTOL=5 HA=$HA DEC=$DEC AZ=$AZ SCALE=$SCALE CX=$CX CY=$CY QUAD=$QUAD CUBE=$CUBE RSIG=$RSIG LNG=$LNG LAT=$LAT FRAT=$frat CLEAN=$CLEAN VERB=$VERB
}}}

