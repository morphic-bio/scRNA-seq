#!/bin/bash
setup_alevin() {
    [ -n "$alevin_setup" ] && return
    export AF_SAMPLE_DIR=$indices_dir/alevin
    export ALEVIN_FRY_HOME="$AF_SAMPLE_DIR/af_home"
    mkdir -p $AF_SAMPLE_DIR
    mkdir -p $ALEVIN_FRY_HOME
    simpleaf set-paths
    alevin_setup=1
}

run_star() {
    for fastq_dir in "${fastq_array[@]}"; do
        echo "find -L ${fastq_dir} -name '*$R1pattern*.$fastq_suffix' -type f | sort | tr '\n' ',' | sed 's/,$//'"
        reads1=$(find -L ${fastq_dir} -name "*$R1pattern*.$fastq_suffix" -type f | sort | tr '\n' ',' | sed 's/,$//')
        reads2=$(find -L ${fastq_dir} -name "*$R2pattern*"$fastq_suffi -type f | sort | tr '\n' ',' | sed 's/,$//')
        echo "$reads1"
        echo "$reads2"
    done

}

run_kallisto() {
    local nac_out_dir="$indices_dir/kallisto/nac_offlist_1"
    [ -z "$overwrite" ] && [ -f "$nac_out_dir/c2" ] && return
    rm -rf /data/kallisto-temp
    local cmd="kb ref --tmp /data/kallisto-temp --d-list=$genome_file -t $n_threads -i $nac_out_dir/index.idx  --workflow nac --overwrite -f1 $nac_out_dir/f1 -f2 $nac_out_dir/f2 -c1 $nac_out_dir/c1 -c2 $nac_out_dir/c2 -g $nac_out_dir/g $genome_file $gtf_file "    
    
    mkdir -p "$nac_out_dir"
    echo "$cmd" 
    eval $cmd
}
run_salmon() {
    local out_dir="$indices_dir/salmon"
    [ -z "$overwrite_index" ] && [ -d "$out_dir/splicei/index" ] && return
    setup_alevin
    mkdir -p $indices_dir/salmon/splicei
    local cmd="simpleaf set-paths && simpleaf index -f $genome_file -g $gtf_file -o $out_dir/splicei -t $n_threads --ref-type splici --no-piscem --overwrite"
    echo "$cmd"
    eval $cmd
}

run_piscem() {
    local out_dir="$indices_dir/salmon"
    [ -z "$overwrite_index" ] && [ -d "$out_dir/spliceu/index" ] && return
    setup_alevin
    mkdir -p $indices_dir/salmon/spliceu
    local cmd="simpleaf set-paths && simpleaf index -f $genome_file -g $gtf_file -o $out_dir/spliceu -t $n_threads --ref-type spliceu --no-piscem --overwrite"
    echo "$cmd"
    eval $cmd
}
parse_string_into_array(){
    local string=$1
    local -n array=$2
    # Remove square brackets and initial and final double quotes
    array=($(echo "$string" | sed 's/[][]//g; s/^"//; s/"$//'))
}

[ -z "$useStar" ] && [ -z "$useKallisto" ] && [ -z "$useSalmon" ] && [ -z "$usePiscem" ] && echo "Must choose one aligner" && exit 1
if [ -n "$usePiscem" ] || [ -n "$useSalmon" ]; then
    [ -z "$useSplicei" ] && [ -z "$useSpliceu" ] && echo "Must choose one reference for salmon/piscem" && exit 1
fi
# check if fastq_suffix ends in .gz
[[ "$fastq_suffix" != *.gz ]] && [ -n "$gzipped_fastq" ] && fastq_suffix="$fastq_suffix.gz"
# parse fastqdirs into an array
parse_string_into_array "$fastqdirs" fastq_array

[ -n "$useStar" ] && run_star
exit 0
[ -n "$useKallisto" ] && run_kallisto
[ -n "$useSalmon" ] && run_salmon
[ -n "$usePiscem" ] && run_piscem   
