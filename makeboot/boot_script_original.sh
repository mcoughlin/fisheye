 LAT=-30.67833 LNG=-70.73669
 utcobs=56784.015602
 gmst=`echo $utcobs | awk '{gmst=(280.460618+360.985647366*($1-51544.5))%360.0; gmst=gmst<0?gmst+360:gmst; printf "%.6f", gmst}'`
 LST=`echo $gmst $LNG | awk '{lst=$1+$2; lst=lst>0?lst:lst+360; printf "%.6f\n", lst}'`

 awk -v lst=$LST '{printf "%7.1f %7.1f %9.4f %9.4f %9.4f  %s\n", $1+0.5,$2+0.5,lst-$3,$4,$3,$9}' boot.ra.original > boot.ha.original

 fisheye -verb < boot.ha.original
#HA=0.015 DEC=-27.092 AZ=95.082 SCALE=0.04705 PARITY=1 MED=1.054

