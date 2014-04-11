#!/bin/bash 

echo "Setting up CMSSW runtime environment"
eval `scramv1 ru -sh`
source $CMSSW_BASE/src/FinalStateAnalysis/environment.sh

export MEGAPATH=$hdfs
source jobid.sh
export jobid=$jobid8
export pub=$HOME/public_html/tauIdEffMeas/$jobid/
export results=$CMSSW_BASE/src/TauIDEfficiencyMeasurement/results/$jobid/
export plots=$results/plots/