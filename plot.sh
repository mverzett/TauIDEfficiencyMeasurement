#!/bin/bash
# Plot all of the analysis

set -o nounset
set -o errexit

source jobid.sh
export jobid=$jobid7

export jobid=$jobid8
python TauEffPlotter.py