#!/bin/bash
# author Tiago Quintino
# 
# $1, first parameter is the file name from where ti remove the data
# $2, second is the number of lines to skip at the beginning, usefull to skip
# format headers or ignore initial data
# $3, the number of the column to actual print, the X axis is always the first
# column

if [ -z $3 ]; then
  echo "Usage: convergenc-plot.sh file nbskip nbcolumn"
  exit 0
fi

TMPCONVFILE=.tmp.$(date +%H%M%S%d%m%Y)
tail -n $(($(cat $1 | wc -l)-$2)) $1 | awk '{print $1, $'$3'}' > $TMPCONVFILE
echo "p '$TMPCONVFILE' w l" | gnuplot -persist
rm $TMPCONVFILE

