linux=`uname -a | grep -i linux`
osx=`uname -a | grep -i darwin`

# Known OS?
if [ "$linux" == "" -a "$osx" == "" ] ; then
  echo "Unsupported operating system?"
  uname -a
  exit 1
fi

# Convert generic "Tue Aug 6 05:38:41 2013" date to YYMMDD
# Syntax: datetoymd "date" zonehr
# E.g:    datetoymd "Tue Aug  6 19:00:00 HST 2013" -10
if [ "$linux" != "" ] ; then
  function datetoymd() { date -d "$1 $((-$2)) hour" +%y%m%d ; }
elif [ "$osx" != "" ] ; then
  function datetoymd() { 
    date -v `printf "%+dH" $2` -j -f "%a %b %d %T %Y" "$1""+%y%m%d"
  }
fi
