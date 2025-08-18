#!/bin/bash
#will do a directory based recursive upload - the file level parallelism will not work well with low bandwidths esp on servers that are not EC2
#sync is now working as fast as the recursive copy so use sync
awsdir=$1
DIR=$2
bucket=$3
glob=$4
#do not nuke existing credentials
[ -d /root/.aws/ ] || cp -r "$awsdir" /root/.aws/

#with one thread just copy the directory
if [ -z "${NTHREADS}" ]; then
 basedir=$(basename "$DIR")
 aws s3 sync "$DIR" "s3://$bucket/$glob/$basedir" 
 exit $?
fi

runJob(){
 lasti=$((${#dirs[@]} - 1))
 for i in $(seq 0 ${lasti}); do
  if (mkdir $lockDir/lock$i 2> /dev/null ); then
   dir=${dirs[$i]}
   echo thread $1 working on $dir
   echo "cd $DIR && nice aws s3 sync $dir s3://$bucket/$glob/$dir"
   cd $DIR && nice aws s3 $copy $dir s3://$bucket/$glob/$dir --recursive
  fi
 done
 exit
}
dirs=( $(cd "$DIR" && find * -maxdepth 0 -type d) )
#if there are directories
if (( ${#dirs[@]} )); then
	lockDir=/tmp/locks.$$
	mkdir -p $lockDir
    for i in $(seq 2 $NTHREADS); do
	  runJob "$i" &
	done
	runJob 1 &
	wait
	rm -rf $lockDir
else
	 basedir=$(basename "$DIR")
 	aws s3 sync "$DIR" "s3://$bucket/$glob/$basedir" 
	exit $?
fi




