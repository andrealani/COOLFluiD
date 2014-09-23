: Ms Windows update script
: Updates to the HEAD version the kernel and plugin sources 
: Note that you must install the svn command line tool
: See the coolfluid webpage for more instructions

@echo off

: update kernel sources
svn up

: go to plugins dir
pushd plugins

: update plugin sources
for /f "Tokens=*" %%G in ('dir /B') do svn up "%%G"

: back to base dir
popd