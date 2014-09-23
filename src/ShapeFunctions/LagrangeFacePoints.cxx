// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/GeometricEntityProvider.hh"

#include "ShapeFunctions/ShapeFunctions.hh"

#include "ShapeFunctions/LagrangeShapeFunctionPointP1.hh"
#include "ShapeFunctions/LagrangeShapeFunctionPointP0.hh"

#include "Framework/Face.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////
// POINT
//////////////////////////////////////////////////////////////////////////////

  /// Lagrange Point Face with P1 geometry and P1 solution.
  GeometricEntityProvider<Face,
    LagrangeShapeFunctionPointP1,
    LagrangeShapeFunctionPointP1,
                          ShapeFunctionsLib>
  FacePointLagrangeP1LagrangeP1("FacePointLagrangeP1LagrangeP1");

  /// Lagrange Line Face with P1 geometry and P0 solution.
  GeometricEntityProvider<Face,
    LagrangeShapeFunctionPointP1,
    LagrangeShapeFunctionPointP0,
                          ShapeFunctionsLib>
  FacePointLagrangeP1LagrangeP0("FacePointLagrangeP1LagrangeP0");

//////////////////////////////////////////////////////////////////////////////

    } // namespace ShapeFunctions

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
