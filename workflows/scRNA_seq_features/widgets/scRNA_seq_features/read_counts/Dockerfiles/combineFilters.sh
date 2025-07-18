#!/bin/bash

base_input_dir=$1
pattern=$2

[ -z "$base_input_dir" ] && echo "Usage: $0 <base_input_dir> <pattern>" && exit 1
[ -z "$pattern" ] && pattern="counts.h5ad"

for file in $(find $base_input_dir -type f -name "$pattern"); do
    dir=$(dirname $file)
    
    if [ -f "$dir/filtered_counts.h5ad" ]; then
        echo "output exists skipping $file"
        continue
    fi
    echo "python3 combineFilters.py --input_file $file"
    python3 combineFilters.py --input_file $file
    
done
