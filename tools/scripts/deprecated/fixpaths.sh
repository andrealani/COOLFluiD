#!/bin/bash
svn switch https://coolfluidsrv.vki.ac.be/svn/coolfluid/Sources/Kernel/trunk
cd plugins
for i in $( find ./ -maxdepth 1 -type d | grep -v .svn | sed -e 's/\.\///g' ) ; do
  cd $i
  svn switch https://coolfluidsrv.vki.ac.be/svn/coolfluid/Sources/Plugins/$i/trunk
  cd ..
done
