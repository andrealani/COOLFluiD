#!/bin/bash
#
# script for the source code files
#
export kdevfile="COOLFluiD.kdevelop.filelist"

tools/dev/list-projectfiles.sh src plugins testcases cmake tools config doc > $kdevfile

echo "prepare.pl" >> $kdevfile
echo "CMakeLists.txt" >> $kdevfile
echo "CTestConfig.cmake" >> $kdevfile

exit
