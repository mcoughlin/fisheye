! Make a picture of the sky from a fisheye .mch file
! Syntax: 
!
! jpg outname.fig (or tube)
! zpmap mchfilename (or skymap mchfilename)
! hardcopy (if jpg)

set \20 1444   ! x center of the fisheye
set \21  968   ! y center of the fisheye
! set \22 1176   ! radius of the fisheye
set \22 1150   ! radius of the fisheye
set \23    6   ! expansion factor
set \24   13   ! crude zeropoint
set \25    9.5 ! sky faint end
set \26   15   ! sky color factor

define rainbow
color   1    0.00 0.00 0.00
color   2    0.53 0.00 0.76
color   3    0.37 0.00 0.82
color   4    0.04 0.00 0.85
color   5    0.04 0.22 0.78
color   6    0.04 0.39 0.78
color   7    0.03 0.53 0.80
color   8    0.03 0.66 0.83
color   9    0.02 0.73 0.77
color  10    0.00 0.73 0.62
color  11    0.00 0.72 0.45
color  12    0.00 0.71 0.25
color  13    0.28 0.72 0.00
color  14    0.51 0.73 0.00
color  15    0.69 0.75 0.00
color  16    0.85 0.76 0.00
color  17    0.96 0.75 0.00
color  18    1.00 0.66 0.00
color  19    1.00 0.58 0.00
color  20    1.00 0.49 0.00
color  21    1.00 0.40 0.00
color  22    1.00 0.31 0.00
color  23    1.00 0.19 0.00
color  24    1.00 0.00 0.00
color  25    1.00 0.22 0.23
color  26    1.00 0.33 0.34
color  27    1.00 0.42 0.44
color  28    1.00 0.51 0.54
color  29    1.00 0.61 0.64
color  30    1.00 0.75 0.80
color  31    1.00 1.00 1.00
color  1
end

define jpg
print 7 /tmp/foo.fig
location 200 9500 200 9500
set \23 6
end

define tube
set xhsize 1000
set xvsize 1000
terminal 7
location 300 3900 300 3900
set \23 3
end

define getem
set \0 \20 - \22
set \1 \20 + \22
set \2 \21 - \22
set \3 \21 + \22
limits \0 \1 \2 \3
color 31 1 1 0
relocate \20 \21
ptype 100 0
set \0 \22 * -1
expand \0
dot
relocate x1 y1
draw     x2 y1
draw     x2 y2
draw     x1 y2
draw     x1 y1
color 1
expand 1
rainbow

data &1.mch
line 2 10000

xcolumn 11        ! x fisheye
ycolumn 12        ! y fisheye
set n 17 column   ! minst
set m 10 column   ! m
set z m - n       ! zeropoint
set s 19 column   ! inst skymag
ptype 103
pcolumn 0 0 1

set d 7 - m       ! expansion ~ magnitude
set d d / \23
set p d expand
end

define zpmap
getem &1 
set c z - \24      ! color ~ zp
set c c + 2
set c c * 10
set p c color
points
end

define skymap
getem &1 
set c \25 - s      ! color ~ sky
set c c sqrt
set c c * \26
set c c min 31
set p c color
points
end

! hardcopy
! fig2dev -L png -g black /tmp/foo.fig /tmp/foo.png
