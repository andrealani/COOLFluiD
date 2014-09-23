#!/bin/bash
#
# Script for searching the source code for some regular expression
#

if [ -z $1 ]
then
    echo "search-source"
    echo "Script for searching the source code for some regular expression"
    echo "Usage: $0 regexp [dir] [dir] ... "
    echo "dir - directory where to search"
    echo "regexp - regular expression to search in the source files"
    echo
    echo "If no directory is given, the current dir is searched"
    exit 1
fi

REG="$1"
shift 1


D=""
if [ -z $1 ]
then
	D="."
fi

find $D $* \(   \
	        -iname "CMakeLists.txt" -or -iname "*.cmake" \
                -or -iname "*.h" -or -iname "*.hh" \
                -or -iname "*.c"  -or -iname "*.cxx" -or -iname "*.cpp" -or -iname "*.ci" \
		-or -iname "*.f" -or -iname "*.f77" -or -iname "*.f90" -or -iname "*.for" \
		-or -iname "*.CFcase" -or -iname "*.in" \
		-or -iname "*.pl" -or -iname "*.sh" \
           \) -print -follow | grep -v "\.svn/" |  sed "s/ /\\\ /g" | \
xargs egrep -H -n -e "$REG"
exit
