!#/bin/sh

# to run:
# source '/home/yoachim/anaconda2/bin/eups-setups.sh'
# setup sims_skybrightness
# nohup ./run_parallel.sh &
cat run_many_months.txt | xargs -I '{}' -n 1 -P 3 sh -c '{}'
