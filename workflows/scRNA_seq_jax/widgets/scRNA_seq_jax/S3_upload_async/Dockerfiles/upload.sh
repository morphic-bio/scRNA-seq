#!/bin/bash
#will do a directory based recursive upload - the file level parallelism will not work well with low bandwidths esp on servers that are not EC2
awsdir=$1
s3bucket=$2
s3prefix=$3

#check that s3bucket starts with s3:// if not then add it
if [[ $s3bucket != s3://* ]]; then
	s3bucket="s3://$s3bucket"
fi	
function check_sync() {
    local source_dir=$1
	#make sure that source_dir ends with / or else add it
	if [ "${source_dir: -1}" != "/" ]; then	
		source_dir="$source_dir/"
	fi
    local dest_url=$2
	#need to remove all trailing / from dest_url
	dest_url=${dest_url%/}
	local basedir=$(basename "$source_dir")
	local cmd="aws s3 sync $source_dir $dest_url/$basedir"
	# Check if there is a profile
	[ -n "$aws_profile" ] && cmd="$cmd --profile $aws_profile"
	#run cmd with --dryrun and grep that files are being copied
	local output=$($cmd --dryrun) 
	if echo $output | grep -q "(dryrun) upload"; then
		return 1
	fi
	return 0
}
function sync_dir() {
    local source_dir=$1
	#make sure that source_dir ends with / or else add it
	if [ "${source_dir: -1}" != "/" ]; then	
		source_dir="$source_dir/"
	fi
    local dest_url=$2
	#need to remove all trailing / from dest_url
	dest_url=${dest_url%/}
    local basedir=$(basename $source_dir)
    local cmd="aws s3 sync $source_dir $dest_url/$basedir"
    # Check if there is a profile
    [ -n "$aws_profile" ] && cmd="$cmd --profile $aws_profile"
    #run cmd with --dryrun and grep that files are being copied

    local output=$($cmd --dryrun) 

    if echo $output | grep -q "(dryrun) upload"; then
        echo "$cmd"  
        eval  $cmd
    fi
    return 1
}
# Parse the string of directories into an array
parse_string_into_array() {
    local string=$1
    local -n array=$2
    # Remove square brackets and initial and double quotes
    local array_string=$(echo "$string" | sed 's/[][]//g; s/\"//g')
	#split on comma
	IFS=',' read -r -a array <<< "$array_string"
	#make IFS default
	IFS=' '
}
#check that s3bucket, s3prefix and localdirs and done directory are set
if [ -z "$s3bucket" ] || [ -z "$s3prefix" ] || [ -z "$localdirs" ] || [ -z "$doneDir" ]; then
    echo "Error: Missing required environment variables."
	echo "s3bucket=$s3bucket"
	echo "s3prefix=$s3prefix"
	echo "localdirs=$localdirs"
	echo "donedir=$doneDir"
    exit 1
fi
# Check if aws credentials are present
if [ -f "$HOME/.aws/credentials" ]; then
    echo "aws credentials present"
else
    mkdir -p -m 0600 $HOME/.aws	
    cp -r $awsdir/* $HOME/.aws/
fi
# parse localdirs into an array
parse_string_into_array "$localdirs" DIRS
echo "DIRS=${DIRS[@]}"
#default minimum time between checks is 60 seconds
[ -z "$MIN_TIME_BETWEEN_CHECKS" ] && MIN_TIME_BETWEEN_CHECKS=60
dest_url="$s3bucket/$s3prefix"
while true; do
	if [ -z "$doneDir" ] || [ -d "$doneDir" ]; then
		#last round of uploads
		echo "last round of uploads"
		for dir in "${DIRS[@]}"; do
			echo "checking $dir"
			if [ -d "$dir" ]; then
				#sudo chown -R lhhung:lhhung $dir
				echo "syncing $dir"	
				sync_dir $dir $dest_url
			fi
		done
		#check if all files are uploaded

		for dir in "${DIRS[@]}"; do
			echo "checking $dir"
			if [ -d "$dir" ]; then
				echo "syncing $dir"	
			   check_sync "$dir" "$dest_url" && exit 1
			fi
		done		
		exit 0
	else
		for dir in "${DIRS[@]}"; do
			if [ -d "$dir" ]; then
				#sudo chown -R lhhung:lhhung $dir
				echo "checking $dir"
				sync_dir $dir $dest_url
			fi
		done	
	fi
    sleep "$MIN_TIME_BETWEEN_CHECKS"
done




