#!/bin/bash

# Check if the number of arguments is as expected
if [ "$#" -ne 1 ]; then
  echo "Usage: source neutrino_check.sh <isTB>"
  exit 1
fi

isTB=$1

# Define an array of commands
if $isTB; then
    root -l -b -x <<EOF
    .L SciFiPlaneView.cpp
    .L USPlaneView.cpp
    .L ConvertedDataAnalyser.cpp
    runAnalysis(100639, 10, true, true, 789)
    runAnalysis(100647, 11, true, true, 997)
    runAnalysis(100672, 8, true, true, 824)
    runAnalysis(100673, 9, true, true, 671)
    runAnalysis(100631, 9, true, true, 871)
    .q
EOF
else
    root -l -b -x <<EOF
    .L SciFiPlaneView.cpp
    .L USPlaneView.cpp
    .L ConvertedDataAnalyser.cpp
    runAnalysis(4752, 51, false, true, 384014)
    runAnalysis(4815, 3, false, true, 511466)
    runAnalysis(4819, 49, false, true, 14045)
    runAnalysis(4992, 10, false, true, 724635)
    runAnalysis(5013, 42, false, true, 493532)
    runAnalysis(5056, 101, false, true, 849600)
    runAnalysis(5099, 75, false, true, 339702)
    runAnalysis(5120, 99, false, true, 612370)
    runAnalysis(4809, 6, false, true, 626367)
    runAnalysis(4976, 64, false, true, 882391)
    runAnalysis(5132, 49, false, true, 831153)
    runAnalysis(5152, 40, false, true, 781267)
    runAnalysis(5171, 18, false, true, 810277)
    runAnalysis(5180, 30, false, true, 578821)
    runAnalysis(5239, 174, false, true, 814223)
    runAnalysis(5389, 50, false, true, 929436)
    runAnalysis(5981, 18, false, true, 71685)
    runAnalysis(6050, 177, false, true, 763312)
    runAnalysis(6069, 73, false, true, 227538)
    runAnalysis(6250, 41, false, true, 301634)
    runAnalysis(6252, 142, false, true, 468058)
    runAnalysis(6268, 38, false, true, 637101)
    runAnalysis(6279, 63, false, true, 954103)
    runAnalysis(6286, 178, false, true, 922363)
    runAnalysis(6290, 116, false, true, 495522)
    runAnalysis(6295, 119, false, true, 394005)
    runAnalysis(6296, 4, false, true, 544667)
    runAnalysis(6296, 9, false, true, 764224)
    runAnalysis(6568, 256, false, true, 40759)
    runAnalysis(6590, 49, false, true, 84534)
    runAnalysis(6596, 192, false, true, 65981)
    runAnalysis(6610, 59, false, true, 826348)
    runAnalysis(6640, 72, false, true, 255465)
    .q
EOF
fi


echo "All commands have finished."
