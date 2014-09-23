#!/bin/bash
if [ -z $1 ]; then
  echo "usage: $0 filename.patch"
  exit -1
fi

svn diff . plugins/* > $1

