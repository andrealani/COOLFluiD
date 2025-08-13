#!/bin/bash
## file was split using:
#split -b 40m corona_720.CFmesh corona_720.CFmesh.  
## this reassembles the original file
cd Mesh
unxz corona_720.CFmesh.*
cat corona_720.CFmesh.?? > corona_720.CFmesh
cd ..
cd MapData	
unxz map_gong_lmax25*
cd ..
