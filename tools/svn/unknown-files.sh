#!/bin/bash
svn st . | egrep "^\?" | cut -b 8-
for dd in $(find plugins -type d | grep -v \.svn )
do 
   svn st $dd | egrep "^\?" | cut -b 8- | sort | uniq -u
done
