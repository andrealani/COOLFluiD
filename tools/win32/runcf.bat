: Ms Windows script to run coolfluid
: Must copy and adapt:
:  - the CF_BUILDDIR to your case
:  - the CF_BUILDTYPE to your case

@echo off

echo PATH [%PATH%]

: ADAPT THIS ==========================

: your build directory
set CF_BUILDDIR=C:\Work\coolfluid\dev\ufem\builds\vs2008x86
: your build type as selected inside Visual Studio
set CF_BUILDTYPE=Debug

: =====================================

: check for the cfcase file nae
if [%1]==[] ECHO Missing CFcase file name. Usage 'runcf.bat CFcaseFile'
if [%1]==[] EXIT /B

: set some variables
set PATH=%CF_BUILDDIR%\dso;%PATH%
set CF_SOLVER=%CF_BUILDDIR%\src\Solver\%CF_BUILDTYPE%\coolfluid-solver.exe
set CF_SOLVERXML=%CF_BUILDDIR%\src\Solver\coolfluid-solver.xml
set CF_CASE=%1

: execute it
%CF_SOLVER% --conf %CF_SOLVERXML% --scase %CF_CASE%
