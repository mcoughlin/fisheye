
#python2.7 run_fisheye_single.py --doMakeFits --doFisheye --doClouds --doCopyFilesSingle --file /Volumes/ALLSKY1/ut042914/ut042914.0500.long.cr2.gz --doGetLatest

#python2.7 run_fisheye_single.py --doMakeFits --doFisheye --doClouds --doCopyFilesSingle --file /Volumes/2TB/ut070314/ut070314.0500.long.cr2.gz

source setup_paths_single.sh
python2.7 run_fisheye_single.py --doMakeFits --doFisheye --doClouds --doCopyFilesSingle --doGetLatest

