#1/bin/bash

# print help if needed
if [ -z "$1" ] ; then
  echo "Usage: $0 TAG [test]"
  exit 0
fi

# enable dry run
if [ "$2" = "test" ] ; then
  export RUN="echo"
else
  export RUN=""
fi

echo ""
echo "+++ Svn::Copy Kernel/trunk -> Kernel/tags/$1"
$RUN svn copy Kernel/trunk Kernel/tags/$1
for i in $( svn list https://ardisksrv1.private.vki.eu/svn/coolfluid/Sources/Plugins );
do

  echo ""
  # plugin exists in this revision
  if [ -d Plugins/$i ]; then 

	# tags dir exists if not create it
  	if [ ! -d Plugins/$i/tags ] ; then
		svn mkdir Plugins/$i/tags
  	fi
	# tags dir exists if not create it
  	if [ ! -d Plugins/$i/branches ] ; then
		svn mkdir Plugins/$i/branches
  	fi
  	echo "+++ Svn::Copy Plugins/$i/trunk -> Plugins/$i/tags/$1"
  	$RUN svn copy Plugins/$i/trunk Plugins/$i/tags/$1
  else
	echo "=== Plugin $i does not exist on this revision" 
  fi
done
