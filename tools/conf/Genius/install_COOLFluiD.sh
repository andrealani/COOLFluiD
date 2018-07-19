# Updated procedure 2018-07-03
# The same CoolFluid_VSC.conf should be used

if [ -z "$1" ] ; then
    echo '"Please choose and run one of the two following options:"'
    echo '"./install_COOLFluiD.sh DEBUG  (with debugging, development version)"'
    echo '"./install_COOLFluiD.sh OPTIM  (w/o debugging, production version)"'
    exit 1 
fi

module load CMake/3.10.2-GCCcore-6.4.0
module load paralution/1.1.0-foss-2018a
module load Boost/1.66.0-foss-2018a
module load PETSc/3.9.0-foss-2018a
module load ParMETIS/4.0.3-foss-2018a

export TOP_DIR="${VSC_DATA}"
export COOLFLUID_TOP_DIR="${TOP_DIR}/COOLFluiD_Genius"
#download COOLFluiD
svn co https://github.com/andrealani/COOLFluiD/trunk ${COOLFLUID_TOP_DIR}

if [ "$1" == "DEBUG" ] ; then
# with debugging
export BUILD_MODE=geniuscuda
elif [ "$1" == "OPTIM" ] ; then
# w/o debugging (production mode)
export BUILD_MODE=geniuscudafast
fi

export COOLFLUID_BASEBUILD_DIR="${COOLFLUID_TOP_DIR}/OPENMPI/${BUILD_MODE}"

cp ${COOLFLUID_TOP_DIR}/tools/conf/Genius/COOLFluid_Genius.conf ${TOP_DIR}
export COOLFLUID_CONF_FILE="${TOP_DIR}/COOLFluid_Genius.conf"
export COOLFLUID_INSTALL_DIR="${COOLFLUID_BASEBUILD_DIR}/INSTALL"
export ALL_ACTIVE=1

cd ${COOLFLUID_TOP_DIR}
./prepare.pl --config-file=${COOLFLUID_CONF_FILE} --build=${BUILD_MODE}

cd ${COOLFLUID_BASEBUILD_DIR}/${BUILD_MODE}
make -j 10
make install
