: Ms Windows prepare script
: Creates a directory with the coolfluid build

: no repeated output
@echo off

: check for the build name
if [%1]==[] ECHO Please provide a build name. Usage 'prepare.bat buildname'
if [%1]==[] EXIT /B

: defined the build name
set buildname=%1
set builddir=builds\%buildname%

echo Build Type %buildname%

: ensure build dir exists
md %builddir%

: go to build dir
pushd %builddir%

: execute cmake
cmake ..\.. -G"Visual Studio 9 2008" -DCF_SKIP_FORTRAN:BOOL=1 -DCMAKE_BUILD_TYPE=DEBUG

: go to original directory
popd