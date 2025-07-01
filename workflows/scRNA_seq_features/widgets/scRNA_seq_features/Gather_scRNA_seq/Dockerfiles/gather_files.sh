#!/bin/bash

process_directory_list() {
    #return an array of directories
    local string=$1
    #remove the first and last element of the string if they are double quotes
    string=$(echo $string | sed 's/^"//;s/"$//')
    #remove the first and last element of the string if they are [ ]
    string=$(echo $string | sed 's/^\[//;s/\]$//')
    #split the string by comma
    IFS=',' read -r -a array <<< "$string"
    array=("${array[@]//\"}")
    echo ${array[@]}
    #restore the IFS
    IFS=$'\n'
}
#check that the output_dir environment variable is set
if [ -z "$output_dir" ]; then
    echo "output_dir environment variable is not set"
    exit 1
fi
mkdir -p $output_dir
#check that the counts_dir environment variable is set
if [ -z "$counts_dir" ]; then
    echo "counts_dir environment variable is not set"
    exit 1
fi
directory_list=$(process_directory_list $counts_dir)
#loop through the directory_list
for directory in ${directory_list[@]}; do
    #check that the directory is a directory
    if [ ! -d "$directory" ]; then
        echo "$directory is not a directory"
        exit 1
    fi
done

for directory in ${directory_list[@]}; do
   expression_files=($(find $directory -name "final_counts.h5ad" -type f))
   for expression_file in ${expression_files[@]}; do
        #get the name of the expression file
        expression_file_name=$(basename $expression_file)
        #get the directory of the expression file - should be sample/aligner/cellbender
        sample_dir=$(dirname $(dirname $(dirname $expression_file)))
        #find the aligners which are the directories in the sample_dir
        aligners=($(find $sample_dir -mindepth 1 -maxdepth 1 -type d))
        sample_name=$(basename $sample_dir)
        #if there is more than one aligner, then the expression file is in the cellbender directory
        if [ ${#aligners[@]} -gt 1 ]; then
            #get the name of the cellbender directory
            for aligner in ${aligners[@]}; do       
                #copy the expression file to the output directory
                echo "mkdir -p $output_dir/$sample_name/$aligner"
                echo "cp $expression_file $output_dir/$sample_name/$aligner/$expression_file_name"
            done
        else
         #don't bother making a subdirectory for the aligner
         #just copy the expression file to the output directory
         echo "mkdir -p $output_dir/$sample_name"
         echo "cp $expression_file $output_dir/$sample_name/$expression_file_name"
        fi
    done
done
#check that the counts_dir is a directory