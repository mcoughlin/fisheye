set xhsize 1050
set xvsize 900
terminal

! Run on fish file: astrometry ut030114.0650
! Astrometric residuals
define astrometry
data &1.fish
line 3 100000
xc 1
yc 2
uc 11
vc 12
set d 17 colu
limits
set u u - x 
set v v - y
expand 1
lweight 1
color 3 1 0 0
vfield 1 100
color 1
points
color 2 0 1 0
grid
color 4 0 1 1
relocate 1444 968
expand -1100
ptype 100 0
dot
expand 1
ptype 40
color 1
expand 1.5
lweight 2
box
xlabel x
ylabel y
id
expand 1
lweight 1
end

define cubicdropoff
set x r / 1000
set y m - n
limits 0 0.4 \0 \1
poly 0 0
set \2 \12      ! Constant
set y \2 - y
set y y ln
set y y * 0.3333
set y y exp
set y y / x
limits
limits 0.8 x2 0.5 1
poly 0 0
set \3 \12     ! poly
set x r
limits
set x 0 to x2
set y x / 1000
set y y * \3
set y y ln
set y y * 3
set y y exp
set y y * -1
set y y + \2
limits 0 x2 \0 \1
connect
end

! Run on mch file: photometry ut030114.0650
define photometry
data &1.mch
lines 3 100000
set m 11 column
set n 18 column
set e 19 column
set r 17 column
set a 18 column

set x r
set y m - n
limits
poly 0 0
set \0 \12 - 2    ! lower zp limit
set \1 \12 + 2    ! upper zp limit

! cubicdropoff

limits 0 x2 \0 \1
set x r
set y m - n
color 3 0 1 0
grid
color 2 1 0 0
points
color 1
expand 1.5
lweight 2
box
ylabel m - minst
xlabel radius [pix]
id
expand 1
lweight 1
end
