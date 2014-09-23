cd ./vkiles
rm -r ./lib
mkdir ./lib
cd ./lib
echo "CMake used: $(which cmake)"
cmake ..
make
echo "***************************************************************************"
echo "* WARNING: Don't forget to add the following line to your coolfluid.conf: *"
echo "* lesmodels_dir = /---path-to-coolfluid---/plugins/UFEM/les               *" 
echo "***************************************************************************"