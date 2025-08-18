#!/bin/bash

awsdir=$1
bucket=$2
outputDir=$3

error=0
mkdir -p "$outputDir"
[ -d /root/.aws/ ] || cp -r "$awsdir" /root/.aws/

copy_wildcard(){
 echo "parsing wildcards in $1"
 local command attempts max_attempts my_str my_glob
 my_str=$1
 max_attempts=4
 #if there is no / then we search from the base bucket
 my_glob=""
 wildcard=$my_str
 if [[ $glob == */* ]]; then
	#split into a glob (directory) and wildcard string
	no_wc="${my_str%%['!'@#\$%^\&*()+]*}"
	my_glob="${no_wc%/*}"/
	wildcard="${my_str#$my_glob}"
 fi
 if [ -n "$overwrite" ]; then 
   command=(nice aws s3 cp --exclude "*" --include="$wildcard" --recursive s3://$bucket/$my_glob $outputDir)
 else
   cd $outputDir
   command=(nice aws s3 sync --exclude "*" --include="$wildcard" s3://$bucket/$my_glob .)
 fi
 for attempts in {1..$max_attempts}; do
	echo "${command[@]}"
 	if "${command[@]}" ; then
   	return
 	fi
 done
 echo "error in ${command[@]}"
 exit 1

}

copy_directory(){
 echo "copying directory object $1"
 local attempts command max_attempts
 max_attempts=2
 if [ -n "$overwrite" ]; then
    command=(nice aws s3 cp --recursive s3://$bucket/$1 $outputDir) 
 else
    cd "$outputDir"
    command=(nice aws s3 sync s3://$bucket/$1 .) 
 fi
 for attempts in {1..$max_attempts}; do
	echo "${command[@]}"
 	if "${command[@]}" ; then
   	return
 	fi
 done
 echo "error in ${command[@]}"
 exit 1

}

copy_file(){
 echo "copying file object $1"
 local destination attempts command max_attempts
 max_attempts=2
 destination=$(basename "$1")
 if [ -n "$overwrite" ]; then
    command=(nice aws s3 cp s3://$bucket/$1 $outputDir/$destination) 
 else
    cd "$outputDir"
    command=(nice aws s3 sync s3://$bucket/$1 .) 
 fi
 for attempts in {1..$max_attempts}; do
	echo "${command[@]}"
 	if "${command[@]}" ; then
   	return
 	fi
 done
 echo "error in ${command[@]}"
 exit 1
}

copy(){
   local my_glob=$1
   echo "$my_glob"
   if [[ $my_glob == *['!'@#\$%^\&*()+]* ]]; then
     copy_wildcard $my_glob || error=1
   elif [ "${my_glob: -1}" == "/" ]; then
     copy_directory $my_glob || error=1
   else
    copy_file $my_glob || error=1
   fi	
}

multiCopy(){
 lasti=$((${#globs[@]} - 1))
 for i in $(seq 0 ${lasti}); do
  if ( mkdir $lockDir/lock$i 2> /dev/null ); then
   glob=${globs[i]}
   echo "thread $1 copying $glob"
   copy $glob
  fi
 done
}
if [ -z $DIRS ] || [ "$DIRS" == "[]" ]; then
    echo "no bucket object given to download"
	exit 1
fi
globs=( $(echo $DIRS | jq -r '.[]') )

if [ -z $nThreads ] || (( $nThreads == 1 )) || (( $nThreads == 0 )); then
	#use single thread
	echo "Using single thread"
	for glob in "${globs[@]}"; do
		copy $glob
	done
else
	lockDir=/tmp/locks.$$
	mkdir -p $lockDir
	for i in $(seq 2 $nThreads); do
	  multiCopy $i &
	done
	multiCopy 1 &
	wait
	rm -rf $lockDir
fi
exit $error
