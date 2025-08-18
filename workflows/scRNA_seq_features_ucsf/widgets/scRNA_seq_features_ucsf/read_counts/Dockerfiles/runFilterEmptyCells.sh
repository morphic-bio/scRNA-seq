#!/bin/bash
dir=$1
pattern=$2
if [ -z "$dir" ]; then
  echo "Usage: $0 <dir>"
  exit 1
fi
if [ ! -d "$dir" ]; then
  echo "Error: $dir is not a directory"
  exit 1
fi
if [ -z "$pattern" ]; then
  echo "Usage: $0 <pattern>"
  exit 1
fi
#find all files in the directory that match the pattern
for file in $(find $dir -type f -name "$pattern"); do
  echo "working on $file"
  #find parent directory of the file and make a logfile name
  logfile=$(dirname $file)/filter_empty_cells.log 
  echo "Logging to $logfile"
  #tee both stderr and stdout of the R script to the logfile
  filterEmptyCells.R $file 2>&1 | tee $logfile
done