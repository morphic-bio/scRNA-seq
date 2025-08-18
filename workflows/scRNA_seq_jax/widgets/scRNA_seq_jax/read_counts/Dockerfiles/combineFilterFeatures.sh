#!/bin/bash

base_input_dir=$1
pattern=$2

[ -z "$base_input_dir" ] && echo "Usage: $0 <base_input_dir> <pattern>" && exit 1
[ -z "$pattern" ] && pattern="unfiltered_counts.h5ad"

for file in $(find $base_input_dir -type f -name "$pattern"); do
    dir=$(dirname $file)
    #replace raw with filtered in the output filename
    output_file="$dir/filtered_hq_features.h5ad"
    #if [ -f "$output_file" ]; then
    #    echo "output exists skipping $file"
    #    continue
    #fi
    echo "python3 combineFiltersforFeatures.py --input_file $file "
    python3 combineFiltersforFeatures.py --input_file $file
    
done