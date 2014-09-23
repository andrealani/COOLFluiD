#!/bin/bash
#
# Script for searching the source and replacing a string for another
#

if [ -z $2 ]
then
    echo "search-replace"
    echo "Script for searching the source code for a string and replacing for another"
    echo "Ignores upper or lower case characters."
    echo "Usage: $0 string1 string2 [dir]"
    echo "string1 - original string to search in the source files"
    echo "string2 - string to replace into the source files"
    echo "dir - top directory to start replacing"
    echo "      If no directory is given, the current dir is searched"
    exit 1
fi

# make th default dir search the local one
D=""
if [ -z $3 ]
then
	D="./"
fi

# take the first two parameters and the rest are directories
FF=$1
RR=$2
shift 2

files=$( find $D $* \
	        -iname "CMakeLists.txt" -or -iname "*.cmake" \
                -or -iname "*.h" -or -iname "*.hh" \
                -or -iname "*.c"  -or -iname "*.cxx" -or -iname "*.cpp" -or -iname "*.ci" \
		-or -iname "*.f" -or -iname "*.f77" -or -iname "*.f90" -or -iname "*.for" \
		-or -iname "*.CFcase" -or -iname "*.in" \
		-or -iname "*.pl" -or -iname "*.sh" \
       )
for filename in $files
do
  perl -pi -e "s#$FF#$RR#gi" $filename
done


