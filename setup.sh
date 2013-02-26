#!/bin/bash

# Setup the cython proxies, find input ntuple files, and compute luminosity.

source jobid.sh
## export jobid=$jobid7
## export datasrc=$(ls -d /scratch/*/data/$jobid | awk -F$jobid '{print $1}')
## 
## for dir in $datasrc; do
##     rake "meta:getinputs[$jobid, $datasrc]"
##     rake "meta:getmeta[inputs/$jobid, mm/metaInfo, 7]"
## done
## 
export jobid=$jobid8
export datasrc=$(ls -d /scratch/*/data/$jobid | awk -F$jobid '{print $1}')
./make_proxies.sh
for dir in $datasrc; do
    rake "meta:getinputs[$jobid, $datasrc]"
    rake "meta:getmeta[inputs/$jobid, mm/metaInfo, 8]"
done
# Use the 7TeV WH samples for 8TeV
#pushd inputs/$jobid/
# Symlink the list of input files and the counts of the number of events.
# For the effectively lumis, we have to recompute using the 8 TeV x-section.
#ls ../../inputs/$jobid7/WH_*HWW* | grep -v lumicalc | xargs -n 1 ln -s 
#popd
