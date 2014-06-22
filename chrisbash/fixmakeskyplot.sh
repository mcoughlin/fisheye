#!/bin/bash
cd $fixdirpath/$fixdirname
rm $fixdirname.sky.eps
grep -v '#' $fixdirname.short.obslog | awk '{print $4, $5, $6, $7, $8, $9}' > temp.dat
sm << EOF
data temp.dat
read MJD 1
read Msky 2
read Bsky 3
read Gsky 4
read Rsky 5
read nstars 6
device postencap Gsky.eps
set Gskymag=2.5*lg(Gsky)
set Rskymag=2.5*lg(Rsky)
set Bskymag=2.5*lg(Bsky)
limits MJD Rskymag
box
ctype green
points MJD Gskymag
ctype blue
points MJD Bskymag
ctype red
points MJD Rskymag
ctype default
xlabel MJD
ylabel magnitudes
toplabel $fixdirname
hardcopy
EOF
mv Gsky.eps $fixdirname.sky.eps
rm temp.dat
