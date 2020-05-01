# Updated procedure 2018-07-03
# The same CoolFluid_VSC.conf should be used

if [ -z "$1" ] ; then
    echo '"Please choose and run one of the following options:"'
    echo '"./install_COOLFluiD.sh $1"'
    echo '"where $1 is either --download=0 (update) or --download=1 (update,configure) or --download=2 (download,update,configure)"' 
    exit 1 
fi

module load CMake/3.10.2-GCCcore-6.4.0
module load PETSc/3.9.0-foss-2018a-cpu
#module load Subversion/1.8.14-foss-2015a
module load Boost/1.66.0-foss-2018a 
#module load Boost/1.70.0-foss-2018a
module load ParMETIS/4.0.3-foss-2018a

export TOP_DIR="${VSC_DATA}"
export COOLFLUID_TOP_DIR="${TOP_DIR}/COOLFluiD_VSC"
#download COOLFluiD
if [ "$1" == "--download=2" ] ; then
svn co https://github.com/andrealani/COOLFluiD/trunk ${COOLFLUID_TOP_DIR}
elif [ "$1" == "--download=0" ] || [ "$1" == "--download=1" ] ; then
#update COOLFluiD
cd ${COOLFLUID_TOP_DIR} 
svn up .
fi

export COOLFLUID_BASEBUILD_DIR="${COOLFLUID_TOP_DIR}/OPENMPI"
export BUILD_MODE=optim
export COOLFLUID_CONF_FILE="${TOP_DIR}/COOLFluid_VSC.conf"
export COOLFLUID_INSTALL_DIR="${COOLFLUID_BASEBUILD_DIR}/INSTALL"
export ALL_ACTIVE=1

if [ "$1" == "--download=2" ] ; then
cp ${COOLFLUID_TOP_DIR}/tools/conf/VSC/COOLFluid_VSC.conf ${TOP_DIR}
cd ${COOLFLUID_TOP_DIR}
./prepare.pl --config-file=${COOLFLUID_CONF_FILE} --build=${BUILD_MODE}
elif [ "$1" == "--download=1" ] ; then
# clean up old object files and libraries
rm -fr ${COOLFLUID_BASEBUILD_DIR}/${BUILD_MODE} 
cd ${COOLFLUID_TOP_DIR}
./prepare.pl --config-file=${COOLFLUID_CONF_FILE} --build=${BUILD_MODE}
fi

cd ${COOLFLUID_BASEBUILD_DIR}/${BUILD_MODE}
make -j 4
make install
