#!/bin/bash
source global_env.sh
set -ex

pair=$(sed -n "$1 p" $2)
sname=$(basename $pair | awk -F".fq.gz" '{print $1}')

#cd ${home_project}
treads=$(rapidgzip -P 4 --count-lines ${home_project}/$pair | awk '{print $1/4}')
echo -e "$sname\t$treads" >> ~/sandro/phage_approach/supp_files/summary_samples.tsv
