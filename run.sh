#!/bin/bash
# Run all of the analysis

set -o nounset
set -o errexit

source jobid.sh

export jobid=$jobid8
rake mm
rake mt
