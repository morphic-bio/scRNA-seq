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
create_kallisto_index() {
    local nac_out_dir="$indices_dir/kallisto/nac_offlist_1"
    [ -z "$overwrite" ] && [ -f "$nac_out_dir/c2" ] && return
    rm -rf /data/kallisto-temp
    local cmd="kb ref --tmp /data/kallisto-temp --d-list=$genome_file -t $n_threads -i $nac_out_dir/index.idx  --workflow nac --overwrite -f1 $nac_out_dir/f1 -f2 $nac_out_dir/f2 -c1 $nac_out_dir/c1 -c2 $nac_out_dir/c2 -g $nac_out_dir/g $genome_file $gtf_file "    
    
    mkdir -p "$nac_out_dir"
    echo "$cmd" 
    eval $cmd
}
create_star_index() {
    local out_dir="$indices_dir/star"
    [ -z "$overwrite" ] && [ -f "$out_dir/SAindex" ] && return
    mkdir -p "$out_dir"
    local cmd="STAR --runThreadN $n_threads --runMode genomeGenerate --genomeDir $out_dir --genomeFastaFiles $genome_file --sjdbGTFfile $gtf_file --sjdbOverhang 100"
    echo "$cmd"
    eval $cmd
}

create_salmon_index_splicei() {
    local out_dir="$indices_dir/salmon"
    [ -z "$overwrite" ] && [ -d "$out_dir/splicei/index" ] && return
    setup_alevin
    mkdir -p $indices_dir/salmon/splicei
    local cmd="simpleaf set-paths && simpleaf index -f $genome_file -g $gtf_file -o $out_dir/splicei -t $n_threads --ref-type splici --no-piscem --overwrite"
    echo "$cmd"
    eval $cmd
}

create_salmon_index_spliceu() {
    local out_dir="$indices_dir/salmon"
    [ -z "$overwrite" ] && [ -d "$out_dir/spliceu/index" ] && return
    setup_alevin
    mkdir -p $indices_dir/salmon/spliceu
    local cmd="simpleaf set-paths && simpleaf index -f $genome_file -g $gtf_file -o $out_dir/spliceu -t $n_threads --ref-type spliceu --no-piscem --overwrite"
    echo "$cmd"
    eval $cmd
}

create_piscem_index_splicei() {
    local out_dir="$indices_dir/piscem"
    [ -z "$overwrite" ] && [ -d "$out_dir/splicei/index" ] && return
    setup_alevin
    mkdir -p $indices_dir/piscem/splicei
    local cmd="simpleaf set-paths && simpleaf index -f $genome_file -g $gtf_file -o $out_dir/splicei -t $n_threads --ref-type splici --use-piscem --overwrite"
    echo "$cmd"
    eval $cmd
}

create_piscem_index_spliceu() {
    local out_dir="$indices_dir/piscem"
    [ -z "$overwrite" ] && [ -d "$out_dir/spliceu/index" ] && return
    setup_alevin
    mkdir -p $indices_dir/piscem/spliceu
    local cmd="simpleaf set-paths && simpleaf index -f $genome_file -g $gtf_file -o $out_dir/spliceu -t $n_threads --ref-type spliceu --use-piscem --overwrite"
    echo "$cmd"
    eval $cmd
}

if [ -n "$useStar" ]; then
    create_star_index
fi
if [ -n "$useKallisto" ]; then
    create_kallisto_index
fi
if [ -n "$useSalmon" ] || [ -n "$usePiscem" ]; then
    [ -n "$useSalmon" ] && [ -n "$useSplicei" ] && create_salmon_index_splicei
    [ -n "$useSalmon" ] && [ -n "$useSpliceu" ] && create_salmon_index_spliceu
    [ -n "$usePiscem" ] && [ -n "$useSplicei" ] && create_piscem_index_splicei
    [ -n "$usePiscem" ] && [ -n "$useSpliceu" ] && create_piscem_index_spliceu
fi
