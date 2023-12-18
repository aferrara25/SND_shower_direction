#!/bin/bash

SNDLHC_soft=/afs/cern.ch/user/f/fmei/private/sndas
export ALIBUILD_WORK_DIR=$SNDLHC_soft/sw #for alienv
source /cvmfs/sndlhc.cern.ch/SNDLHC-2023/Aug30/setUp.sh
eval `alienv load --no-refresh sndsw/latest`
cd /afs/cern.ch/work/f/fmei/private/sndlhc_bo_tbanalysis

runN=$1 # read from the input_args.txt
nFiles=$2 # read from the input_args.txt
isTB=$3 # read from the input_args.txt

echo "runAnalysis($runN, $nFiles, $isTB)"

root -l -b -x <<EOF
.L SciFiPlaneView.cpp
.L USPlaneView.cpp
.L ConvertedDataAnalyser.cpp
runAnalysis($runN, $nFiles, $isTB)
.q
EOF

echo “Job is done”