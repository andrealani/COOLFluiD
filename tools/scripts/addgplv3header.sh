#!/bin/bash

EXTENSIONS_ARRAY=(.hh .ci .cu .cxx)
#DIRS_ARRAY=(./AnalyticalModel ./BackwardEuler ./Burgers ./CFmesh2THOR ./CFmeshCellSplitter ./CFmeshExtruder ./CFmeshFileReader ./CFmeshFileWriter ./CGNS2CFmesh ./ConvertStructMesh ./Dpl2CFmesh ./EmptySpaceMethod ./FAST2CFmesh ./ForwardEuler ./Gambit2CFmesh ./Gmsh2CFmesh ./LinearAdv ./LinearAdvSys ./MeshGenerator1D ./MeshTools ./NewtonMethod ./NonLinearAdv ./ParaViewWriter ./Pardiso ./Petsc ./PhysicalModelDummy ./RotationAdvSys ./SAMGLSS ./SimpleGlobalMeshAdapter ./TecplotWriter ./Tecplot2CFmesh ./Trilinos ./THOR2CFmesh ./RungeKutta ./RungeKutta2 ./RungeKuttaLS ./XCFcaseConverter)
DIRS_ARRAY=(./RungeKuttaLS)

# loop on extensions
for item in ${EXTENSIONS_ARRAY[*]} ; do
  echo "Processing extension: $item"
for i in $(tree --noreport -fi ${DIRS_ARRAY[*]} | grep $item) ; do 
    echo "Processing: $i"
    mv $i $i.change
    echo "// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium" > $i
    echo "//" >> $i
    echo "// This software is distributed under the terms of the" >> $i
    echo "// GNU Lesser General Public License version 3 (LGPLv3)." >> $i
    echo "// See doc/lgpl.txt and doc/gpl.txt for the license text." >> $i
    echo "" >> $i
    cat $i.change >> $i
    rm $i.change
  done
done
