# Updated procedure 2018-07-03
# The same CoolFluid_VSC.conf should be used

if [ -z "$1" ] ; then
    echo '"Please choose and run one of the two following options:"'
    echo '"./install_COOLFluiD.sh DEBUG_CUDA   (with debugging, development version with CUDA)"'
    echo '"./install_COOLFluiD.sh OPTIM_CUDA   (w/o debugging, production version with CUDA)"'
    echo '"./install_COOLFluiD.sh DEBUG_NOCUDA (with debugging, development version w/o CUDA)"'
    echo '"./install_COOLFluiD.sh OPTIM_NOCUDA (w/o debugging, production version w/o CUDA)"'
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

if [ "$1" == "DEBUG_CUDA" ] ; then
# with debugging
export BUILD_MODE=geniuscuda
export CONF_FILE="COOLFluid_Genius.conf"
elif [ "$1" == "OPTIM_CUDA" ] ; then
# w/o debugging (production mode)
export BUILD_MODE=geniuscudafast
export CONF_FILE="COOLFluid_Genius.conf"
elif [ "$1" == "DEBUG_NOCUDA" ] ; then
# w/o debugging (production mode)
export BUILD_MODE=optim
export CONF_FILE="COOLFluid_Genius_nocuda.conf"
elif [ "$1" == "OPTIM_NOCUDA" ] ; then
# w/o debugging (production mode)
export BUILD_MODE=release
export CONF_FILE="COOLFluid_Genius_nocuda.conf"
fi

export COOLFLUID_BASEBUILD_DIR="${COOLFLUID_TOP_DIR}/OPENMPI"

cp ${COOLFLUID_TOP_DIR}/tools/conf/Genius/${CONF_FILE} ${COOLFLUID_TOP_DIR}
export COOLFLUID_CONF_FILE="${COOLFLUID_TOP_DIR}/${CONF_FILE}"
export COOLFLUID_INSTALL_DIR="${COOLFLUID_BASEBUILD_DIR}/${BUILD_MODE}/INSTALL"
export ALL_ACTIVE=1

cd ${COOLFLUID_TOP_DIR}
./prepare.pl --config-file=${COOLFLUID_CONF_FILE} --build=${BUILD_MODE}

cd ${COOLFLUID_BASEBUILD_DIR}/${BUILD_MODE}
make -j 10
make install
