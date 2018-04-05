// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/GeometricEntityProvider.hh"

#include "ShapeFunctions/ShapeFunctions.hh"

#include "ShapeFunctions/LagrangeShapeFunctionQuadP2.hh"
#include "ShapeFunctions/LagrangeShapeFunctionQuadP1.hh"
#include "ShapeFunctions/LagrangeShapeFunctionQuadP0.hh"

#include "Framework/Face.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////
// QUAD
//////////////////////////////////////////////////////////////////////////////

  /// Lagrange Quadrilateral Face with P1 geometry and P2 solution.
  GeometricEntityProvider<Face,
                          LagrangeShapeFunctionQuadP1,
                          LagrangeShapeFunctionQuadP2,
                          ShapeFunctionsLib>
  FaceQuadLagrangeP1LagrangeP2("FaceQuadLagrangeP1LagrangeP2");

  /// Lagrange Quadrilateral Face with P1 geometry and P1 solution.
  GeometricEntityProvider<Face,
    LagrangeShapeFunctionQuadP1,
    LagrangeShapeFunctionQuadP1,
                          ShapeFunctionsLib>
  FaceQuadLagrangeP1LagrangeP1("FaceQuadLagrangeP1LagrangeP1");

  /// Lagrange Quadrilateral Face with P1 geometry and P0 solution.
  GeometricEntityProvider<Face,
          LagrangeShapeFunctionQuadP1,
          LagrangeShapeFunctionQuadP0,
                          ShapeFunctionsLib>
  FaceQuadLagrangeP1LagrangeP0("FaceQuadLagrangeP1LagrangeP0");

//////////////////////////////////////////////////////////////////////////////

    } // namespace ShapeFunctions
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
