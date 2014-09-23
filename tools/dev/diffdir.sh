#!/bin/bash

if [ -z $1 ]
then
    echo "diffdir"
    echo "Recursive file difference comparison between files in two directories."
    echo "Usage: $0 dir1 dir2"
    echo "dir1 dir2 - directories to compare"
    exit 1
fi
dir1=$1
dir2=$2
for file in $(ls $dir1 | grep -v CVS)
  do 
  if [ -f $dir2/$file ]; then
      echo "Differences from $file in $dir1 to $dir2"
      diff $dir1/$file $dir2/$file
  else
      echo "No $file in $dir2"
  fi
done



