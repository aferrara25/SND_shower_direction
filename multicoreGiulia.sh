#!/bin/bash

# Check if the number of arguments is as expected
if [ "$#" -ne 3 ]; then
  echo "Usage: source multicoreGiulia.sh <run> <nFiles> <isTB>"
  exit 1
fi

run=$1
nFiles=$2
isTB=$3


# Define an array of commands
commands=()
for ((i=0; i<$nFiles; i++)); do
  commands+=("root -l -b -x <<EOF
.L SciFiPlaneView.cpp
.L USPlaneView.cpp
.L ConvertedDataAnalyser.cpp
runAnalysis($run, $i, $isTB, true)
.q
EOF")
done


# Create temporary script files for each command and make them executable
for ((i=0; i<${#commands[@]}; i++)); do
    temp_script="/tmp/temp_script_$i.sh"
    echo "${commands[i]}" > "$temp_script"
    chmod +x "$temp_script"
    "$temp_script" &
    temp_pids[i]=$!
done

# Wait for all background processes to finish
for pid in "${temp_pids[@]}"; do
    wait "$pid"
done

# Clean up temporary script files
for temp_script in "${temp_scripts[@]}"; do
    rm "$temp_script"
done


echo "Merging files..."

root -l -b -q "sumHisto.C($run, $nFiles, $isTB)"

# Additional commands to execute after parallel tasks are done
echo "All commands have finished."
