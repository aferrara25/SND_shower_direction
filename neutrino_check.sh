#!/bin/bash

# Define an array of commands
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
.q
EOF


echo "All commands have finished."
