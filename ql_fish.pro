! Report some statistics of a fisheye image

! Estimate of saturation level
set SAT=12000
! Radius where masking starts
set RAD=1000
set MAXRAD=1150

rd 1 {arg2}.fits silent

! Count saturated pixels
cop 20 1
clip 20 min=SAT
di 20 20
abx 20 all total=nsat silent

set bias={1:BIAS}
sc 1 bias

set cx={1:NAXIS1}/2 cy={1:NAXIS2}/2
cop 2 1 
surface 2 params=(cx,cy,0.01,0,0,1,1,0) silent
sqrt 2
clip 2 max=RAD
di 2 2 
mi 2 1

! Exposure time and aperture relative to f/4.0
set etime={1:EXPTIME}
set ap={1:APERTURE}/4.0

! Sky level and variation
abx 2 all median=sky medrms=rms silent

! Sky brightness per unit time
set bright=sky*ap*ap/etime

! Tell us the stats...
printf '%f9.2 %f9.2 %f15.3 %i10' sky rms bright nsat

! Write a jpeg stretched to really show the stars
box 1 sx=cx-MAXRAD ex=cx+MAXRAD sy=0 ny={2:NAXIS2}
cop 3 1 box=1
smooth 3 fw=2
cop 3 3 bin=2
set thresh=sky/2 imgsat=10*rms+sky
if imgsat>SAT
  set imgsat=SAT
end_if
tv 3 cf=bw thresh=thresh sat=imgsat jpeg={arg2}_bw.jpg

end
