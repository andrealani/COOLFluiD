#!/bin/bash
for dd in $(find plugins -maxdepth 1 -type d | grep -v \.svn )
do 
  if [ "$(svn st $dd)" != "" ] ; then 
    echo $dd 
  fi
done
