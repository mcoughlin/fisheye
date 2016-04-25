#!/bin/sh

# $1 = month
# $2 = year


# ls -d /lsst/all-sky/ut01*16 | sed 's/.*\///' | xargs -I'{}' ~/gitrepos/fisheye/fast_hp/reduce_night.py '{}'

lsst
setup sims_skybrightness

ls -d /lsst/all-sky/ut$1*$2 | sed 's/.*\///' | xargs -I'{}' ~/gitrepos/fisheye/fast_hp/reduce_night.\
py '{}'
