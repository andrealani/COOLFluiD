#!/bin/bash
#
# Script for searching the source and replacing a string for another
#
echo "tagremove.sh"
echo "Script for removing CVS tags from the source files "
echo "Ignores upper or lower case characters."

files=$(find . -name "*.hh" -or -name "*.cxx" -or -name "*.ci")

for filename in $files
do
  perl -pi -e 's/\ *?\*\ *?\$Revision.*?\$/\ */gi' $filename
  perl -pi -e 's/\ *?\*\ *?\$Source.*?\$/\ */gi' $filename
  perl -pi -e 's/\ *?\*\ *?\$Date.*?\$/\ */gi' $filename
done


