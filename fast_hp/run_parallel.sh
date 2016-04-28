!#/bin/bash

# to run:
# nohup ./run_parallel.sh 
cat run_many_months.txt | xargs -I '{}' -n 1 -P 3 run_month.sh '{}'
