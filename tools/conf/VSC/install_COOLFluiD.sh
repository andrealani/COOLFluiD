# Updated procedure 2018-07-03
# The same CoolFluid_VSC.conf should be used

source switch_to_2015a
module load CMake/3.1.0-foss-2015a

module load OpenMPI/1.8.4-GCC-4.9.2
#module load cURL/7.43.0-foss-2015a
module load Boost/1.66.0-foss-2015a-Python-2.7.9
module load PETSc/3.6.3-foss-2015a-Python-2.7.9
module load ParMETIS/4.0.3-foss-2015a

# extra modules
module load GSL/2.1-foss-2015a
module load METIS/5.1.0-foss-2015a

export TOP_DIR="/data/leuven/304/vsc30484/temp"
export COOLFLUID_TOP_DIR="${TOP_DIR}/YOUR_COOLFluiD"
export COOLFLUID_BASEBUILD_DIR="${TOP_DIR}/YOUR_COOLFluiD/OPENMPI"
export BUILD_MODE=optim
export COOLFLUID_CONF_FILE="${TOP_DIR}/COOLFluid_VSC.conf"
export COOLFLUID_INSTALL_DIR="${TOP_DIR}/COOLFluid_Install_Dir"
export ALL_ACTIVE=1

cd ${COOLFLUID_TOP_DIR}
./prepare.pl --config-file=${COOLFLUID_CONF_FILE} --build=${BUILD_MODE}

cd ${COOLFLUID_BASEBUILD_DIR}/${BUILD_MODE}
make -j 10
make install
