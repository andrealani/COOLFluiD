#!/bin/bash
svn st . | egrep "^C" | cut -b 8-
for dd in $(find plugins -maxdepth 1 -type d | grep -v \.svn )
do 
   svn st $dd | egrep "^C" | cut -b 8-
done
