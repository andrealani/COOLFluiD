#!/bin/bash

if [ -z "$1" ] ; then
    echo 'install.sh needs the full path to the COOLFluiD base directory as full path'
    echo '"Usage: install.sh baseDir"'
    exit 1
fi

export WORKDIR=$1
export DEPSDIR=$1/DEPS

cd $WORKDIR

# set the environment variables for the COOLFluiD dependencies
cp tools/conf/coolfluid.conf.radiation coolfluid.conf.0
sed -e "s|BASEDIR|$WORKDIR|" coolfluid.conf.0 > coolfluid.conf.1
sed -e "s|DEPSDIR|$DEPSDIR|" coolfluid.conf.1 > coolfluid.conf.radiation
rm -fr coolfluid.conf.0 coolfluid.conf.1

cd tools/scripts

cp run.sh.bkp run.sh.0
sed -e "s|BASEDIR|$WORKDIR|" run.sh.0 > run.sh.1
sed -e "s|DEPSDIR|$DEPSDIR|" run.sh.1 > run.sh 
rm -fr run.sh.0 run.sh.1

# installation of dependencies in $DEPSDIR
rm -fr $DEPSDIR/tmp/
./install-coolfluid-deps.pl --install=openmpi,parmetis,boost --install-dir=$DEPSDIR --install-mpi-dir=$DEPSDIR --tmp-dir=$DEPSDIR/tmp

# install mutation++
cd $DEPSDIR/mutation++
mkdir build ; cd build ; cmake .. ; make -j4 ; make install

export PATH=$DEPSDIR/bin:$DEPSDIR/mutation++/install/bin:$PATH
echo 'PATH is $PATH'
export LD_LIBRARY_PATH=$DEPSDIR/lib:$DEPSDIR/mutation++/install/lib:$LD_LIBRARY_PATH
echo 'LD_LIBRARY_PATH is $LD_LIBRARY_PATH'
export MPP_DATA_DIRECTORY=$DEPSDIR/mutation++/data
echo 'MPP_DATA_DIRECTORY is $MPP_DATA_DIRECTORY'

#install PARADEv3.2.1
export COMPILER=GNU
cd $DEPSDIR/parade-3.2.1
make clean ; make

# install COOLFluiD
cd $WORKDIR
./prepare.pl --build=optim --config-file=coolfluid.conf.radiation
cd OPENMPI/optim ; make -j4
# unpack data for HSNB
#cd $WORKDIR/plugins/RadiativeTransfer/RadiationLibrary/Models/HSNB/
#if [ -f "data.tgz" ]; then tar xvfz data.tgz; fi

