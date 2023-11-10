#!/bin/bash

# Check if the number of arguments is as expected
if [ "$#" -ne 3 ]; then
  echo "Usage: ${0} <INPUT_FILE> <OUTPUT_FOLDER> <isTB>"
  exit 1
fi

INPUT_FILE=$1
OUTPUT_FOLDER=$2
isTB=$3


# Define an array of commands
root -l -x -b -q "skim.cpp(\"$INPUT_FILE\", \"$OUTPUT_FOLDER\", $isTB);"
