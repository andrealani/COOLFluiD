module purge
module load cluster/genius/batch
module load CMake/3.20.1-GCCcore-10.3.0
module load Boost/1.76.0-GCC-10.3.0
module load ParMETIS/4.0.3-gompi-2021a
module load PETSc/3.15.1-foss-2021a


export TOP_DIR="${VSC_DATA}"
export COOLFLUID_TOP_DIR="${TOP_DIR}/COOLFluiD_Genius_Mine"

export BUILD_MODE=optim
export CONF_FILE="COOLFluid_Genius_nocuda.conf"
module load PETSc/3.15.1-foss-2021a


export COOLFLUID_BASEBUILD_DIR="${COOLFLUID_TOP_DIR}/OPENMPI"
export COOLFLUID_CONF_FILE="${COOLFLUID_TOP_DIR}/${CONF_FILE}"
export COOLFLUID_INSTALL_DIR="${COOLFLUID_BASEBUILD_DIR}/${BUILD_MODE}/INSTALL"
export ALL_ACTIVE=1

#./prepare.pl --config-file=${COOLFLUID_CONF_FILE} --build=${BUILD_MODE}

echo "Changing: "
echo ${COOLFLUID_BASEBUILD_DIR}/${BUILD_MODE}
cd ${COOLFLUID_BASEBUILD_DIR}/${BUILD_MODE}
make -j 4
make install

