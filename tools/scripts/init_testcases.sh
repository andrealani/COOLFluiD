#!/bin/bash
# unzip large files for COCONUT regression testcases
if [ -e plugins/MHD/testcases/COCONUT/Dipole_WTD/30Rs.CFmesh.xz ]; then
   unxz plugins/MHD/testcases/COCONUT/Dipole_WTD/30Rs.CFmesh.xz 
fi
if [ -e plugins/MHD/testcases/COCONUT/Dipole_Shifted/2008_lmax25.dat.xz ]; then
   unxz plugins/MHD/testcases/COCONUT/Dipole_Shifted/2008_lmax25.dat.xz 
fi
if [ -e plugins/MHD/testcases/COCONUT/Dipole_Shifted/2Rslvl6.CFmesh.xz ]; then
   unxz plugins/MHD/testcases/COCONUT/Dipole_Shifted/2Rslvl6.CFmesh.xz
fi
cd plugins/MHD/testcases/COCONUT/Eclipse/Mesh
if compgen -G "*.xz" > /dev/null; then
## file was split using:
## split -b 40m corona_720.CFmesh corona_720.CFmesh.  
## this reassembles the original file
   unxz corona_720.CFmesh.*
   cat corona_720.CFmesh.?? > corona_720.CFmesh
fi
cd -
cd plugins/MHD/testcases/COCONUT/Eclipse/MapData      
if compgen -G "*.xz" > /dev/null; then
   unxz map_gong_lmax25*
fi
cd -
if [ -e plugins/RadiativeTransfer/testcases/SolarCorona/corona_fullMHD.CFmesh.xz ]; then
   unxz plugins/RadiativeTransfer/testcases/SolarCorona/corona_fullMHD.CFmesh.xz 
fi
echo "######## init_testcases done! ########"
