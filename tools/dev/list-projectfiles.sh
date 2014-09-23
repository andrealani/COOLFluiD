#!/bin/bash
#
# script for the source code files
#

if [ -z $1 ]
then
    echo "list-projectfiles.sh : script for listing the source code files"
    echo
    echo "Usage: $0 [dir] ... "
    echo "dir - directory where to search"
    echo
    exit 1
fi

find $* \( -iname "*.py" -or -iname "*.sh" -or -iname "*.CFcase" -or -iname "CMakeLists.txt" -or -name "*.cmake" -or -iname "*.pl" -or -iname "*.in" -or -iname "*.c" -or -iname "*.hh" -or -iname "*.h" -or -iname "*.cxx" -or -iname "*.ci" -or -iname "*.f90" -or -iname "*.f" -or -iname "*.cpp" \) -print -follow | grep -v "\.svn/"

exit
