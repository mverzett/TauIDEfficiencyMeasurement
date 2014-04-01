#! /usr/bin/env bash

source jobid.sh

publichtml=$HOME/public_html/tauIdEffMeas
mkdir -p $publichtml
export jobid=$jobid8
for file in $(find results/$jobid/ -type f -name '*.png' ); do
    dest=$publichtml/$(echo $file | awk -F$jobid '{print $2}')
    mkdir -p $(dirname $dest)
    cp -u $file $dest
done

for file in $(find results/$jobid/ -type f -name '*.raw_txt' ); do
    dest=$publichtml/$(echo $file | awk -F$jobid '{print $2}')
    mkdir -p $(dirname $dest)
    cp -u $file $dest
done


$DOTBIN/web.py
