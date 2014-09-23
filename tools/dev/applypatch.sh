#!/bin/bash
if [ -z $1 ]; then
  echo "usage: $0 filename.path"
  exit -1
fi

patch -p0 < $1

