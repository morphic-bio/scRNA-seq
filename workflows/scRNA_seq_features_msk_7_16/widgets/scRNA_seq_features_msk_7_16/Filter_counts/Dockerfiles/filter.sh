#!/bin/bash
parse_string_into_array() {
    local string=$1
    local -n array=$2
    # Remove square brackets and initial and final double quotes
    local arrayString=$(echo "$string" | sed 's/[][]//g; s/\"//g')
    # Split on comma
    IFS=',' read -r -a array <<< "$arrayString"
    # Make IFS default
    IFS=' '
}
filter_counts() {
    local counts_file=$1
    #echo "$counts_file"
    outputfile="$(dirname $counts_file)/filtered_counts.h5ad"
    [ -f "$outputfile" ] && [ -z "$overwrite" ] && echo "File exists: $outputfile" && return 0
    echo "will create $outputfile"
    logfile=$(dirname $counts_file)/valid_cells.log 
    echo "Logging to $logfile"
    #tee both stderr and stdout of the R script to the logfile
    #filterEmptyCells.R $counts_file
    findValidCells.R $counts_file
    
    echo "combineFilters.py --input_file $counts_file --adaptive_filter --n_mad $n_mad"
    eval "combineFilters.py --input_file $counts_file --adaptive_filter --n_mad $n_mad"
}

run_with_limit() {
    while [ $(jobs | wc -l) -ge $nthreads ]; do
        # Wait for any background process to finish if we hit the thread limit
        echo $(jobs | wc -l) "jobs running, waiting for a slot"
        wait -n
    done
}

[ -z "$pattern" ] && pattern="counts.h5ad"
[ -z "$nthreads" ] && nthreads=1
[ -z "$n_mad" ] && n_mad=3
parse_string_into_array "$alignsDir" dirs
counts_files=()

for dir in "${dirs[@]}"; do
    counts_files+=($(find $dir -name $pattern -type f | tr '\n' ' '))
done

#run the command in loop with nthreads
for file in "${counts_files[@]}"; do
    run_with_limit
    filter_counts $file &
done
wait