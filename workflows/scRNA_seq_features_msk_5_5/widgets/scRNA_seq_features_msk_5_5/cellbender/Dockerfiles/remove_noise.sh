#!/bin/bash
export PYTORCH_CUDA_ALLOC_CONF=garbage_collection_threshold:0.6,max_split_size_mb:128
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
clear_gpu_memory() {
    echo "Clearing GPU memory"
    # Kill processes using GPU memory, change to `grep python` if only targeting Python processes.
    for pid in $(nvidia-smi | grep -E 'python|cellbender' | awk '{print $5}'); do
        echo "Killing process $pid to free GPU memory"
        kill -9 $pid
    done
    # Alternatively, if using PyTorch, you could clear cache in Python like this:
    python3 -c "import torch; torch.cuda.empty_cache()"
}
repeat_command() {
    local cmd=$1
    local cpu_cmd=$2
    echo "repeating command $cmd"
    if [ -n "$usecpu" ]; then
        eval "$cmd"
        return $?
    fi
    clear_gpu_memory
    eval "$cmd" && return $?
    echo "Command failed, retrying without cuda"
    eval "$cpu_cmd"
    return $? 
}

filter_counts() {
    local counts_file=$1
    local full_cb_file
    local logfile
    local cmd
    local outputfile
    local full_cb_subdir
    full_cb_subdir="$(dirname $counts_file)/$cb_subdir" 
    full_cb_file="$full_cb_subdir/$cb_file"
    outputfile="$(dirname $counts_file)/$cb_subdir/$output_pattern"
    if [ -n "$overwrite_cellbender" ]; then
        [ -d "$full_cb_subdir" ] && rm -rf "$full_cb_subdir" && echo "Overwriting previous results - removed $full_cb_subdir"
    fi
    if [ -f "$full_cb_file" ] && [ -z "$overwrite_cellbender" ]; then
        echo "File $full_cb_file already exists, skipping cellbender"
        addCounts.py -i $file -c $full_cb_file -o $outputfile -l $layername
        return
    fi
    [ -f "$full_cb_file" ] && echo "will create $full_cb_file" || echo "will overwrite $full_cb_file"
    logfile="$full_cb_subdir/cellbender.log" 
    echo "Logging to $logfile"
    [ -d "$full_cb_subdir" ] || mkdir -p "$full_cb_subdir"
    cmd="cellbender remove-background $runMode --input $counts_file --output $full_cb_file $additional_flags --cpu-threads $cpu_cores"
    echo "$cmd"
    eval "$cmd" || repeat_command "$cmd" "cellbender remove-background --force-use-checkpoint --input $counts_file --output $full_cb_file $additional_flags --cpu-threads $cpu_cores"
    #check if the command was successful
    if [ $? -ne 0 ]; then
        echo "retries failed giving up and going to next file"
    else
        echo "Retry if cellbender completed successfully"
        addCounts.py -i $file -c $full_cb_file -o $outputfile -l $layername
    fi
}

run_with_limit() {
    while [ $(jobs | wc -l) -ge $nThreads ]; do
        # Wait for any background process to finish if we hit the thread limit
        echo $(jobs | wc -l) "jobs running, waiting for a slot"
        wait -n
    done
}

[ -z "$usecpu" ] && runMode="--cuda"
[ -z "$max_retries" ] && max_retries=3
[ -z "$input_pattern" ] && input_pattern="filtered_counts.h5ad"
[ -z "$output_pattern" ] && output_pattern="final_counts.h5ad"
[ -z "$nThreads" ] && nThreads=1
[ -z "$cpu_cores" ] && cpu_cores=1
parse_string_into_array "$alignsDir" dirs
counts_files=()

for dir in "${dirs[@]}"; do
    counts_files+=($(find $dir -name $input_pattern -type f | tr '\n' ' '))
done

#run the command in loop with nthreads
for file in "${counts_files[@]}"; do
    run_with_limit
    echo "Processing $file"
    filter_counts $file &
done
wait
