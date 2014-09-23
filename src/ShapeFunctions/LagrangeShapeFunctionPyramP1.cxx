// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "ShapeFunctions/LagrangeShapeFunctionPyramP1.hh"
#include "Common/NotImplementedException.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



   namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////

CFuint LagrangeShapeFunctionPyramP1::_interpolatorID = 0;
RealVector LagrangeShapeFunctionPyramP1::_vec1 = RealVector(DIM_3D);
RealVector LagrangeShapeFunctionPyramP1::_vec2 = RealVector(DIM_3D);
RealVector LagrangeShapeFunctionPyramP1::_vec3 = RealVector(DIM_3D);
std::vector<RealVector> LagrangeShapeFunctionPyramP1::_normals = std::vector<RealVector>(5);

//////////////////////////////////////////////////////////////////////////////

void LagrangeShapeFunctionPyramP1::computeFaceJacobianDeterminant(
          const std::vector<RealVector>& mappedCoord,
          const std::vector<Framework::Node*>& nodes,
          const Framework::IntegratorPattern& pattern,
                std::vector<RealVector>& faceJacobian)
{
  throw Common::NotImplementedException (FromHere(),"LagrangeShapeFunctionPyramP1::computeFaceJacobianDeterminant()");
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace ShapeFunctions



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
