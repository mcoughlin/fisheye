#!/bin/bash
# Make a movie of all of the processed fisheye movies in a directory
# Syntax: fishmovie.sh [color=sky | color=zp]
# (hardly robust, but still useful)

fishdir=/svn/src/atlas/trunk/user/jt/fisheye
nite=`ls *.mch | head -1 | tr . ' ' | awk '{print $1}'`
color=zp

CLEAN=1

eval $@

n=0
for f in *.mch ; do
cat > /tmp/foo.pl <<EOF
input $fishdir/fishsky.pl
jpg /tmp/foo.fig
${color}map `basename $f .mch`
hardcopy wait
end
EOF
  mongo /tmp/foo.pl
  fig2dev -L png -g black /tmp/foo.fig /tmp/foo.png
  convert /tmp/foo.png `printf "fishy-%03d.jpg" $n`
  let n++
#  sleep 1     # Slow down mongo a bit?  Shouldn't be needed with "wait"
done

ffmpeg -f image2 -i fishy-%03d.jpg -sameq $color$nite.mp4

if [[ $CLEAN == 1 ]] ; then rm fishy-???.jpg /tmp/foo.{pl,fig,png} ; fi
