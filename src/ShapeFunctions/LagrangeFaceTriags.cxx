// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/GeometricEntityProvider.hh"

#include "ShapeFunctions/ShapeFunctions.hh"

#include "ShapeFunctions/LagrangeShapeFunctionTriagP0.hh"
#include "ShapeFunctions/LagrangeShapeFunctionTriagP1.hh"
#include "ShapeFunctions/LagrangeShapeFunctionTriagP2.hh"
#include "ShapeFunctions/LagrangeShapeFunctionTriagP3.hh"

#include "Framework/Face.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////
// TRIAG
//////////////////////////////////////////////////////////////////////////////

  /// Lagrange Triangular Face with P1 geometry and P3 solution.
  GeometricEntityProvider<Face,
                          LagrangeShapeFunctionTriagP1,
                          LagrangeShapeFunctionTriagP3,
                          ShapeFunctionsLib>
  FaceTriagLagrangeP1LagrangeP3("FaceTriagLagrangeP1LagrangeP3");

  /// Lagrange Triangular Face with P1 geometry and P2 solution.
  GeometricEntityProvider<Face,
                          LagrangeShapeFunctionTriagP1,
                          LagrangeShapeFunctionTriagP2,
                          ShapeFunctionsLib>
  FaceTriagLagrangeP1LagrangeP2("FaceTriagLagrangeP1LagrangeP2"); 

    GeometricEntityProvider<Face,
                          LagrangeShapeFunctionTriagP2,
                          LagrangeShapeFunctionTriagP2,
                          ShapeFunctionsLib>
  FaceTriagLagrangeP2LagrangeP2("FaceTriagLagrangeP2LagrangeP2");

  /// Lagrange Triangular Face with P1 geometry and P1 solution.
  GeometricEntityProvider<Face,
                                   LagrangeShapeFunctionTriagP1,
                                   LagrangeShapeFunctionTriagP1,
                          ShapeFunctionsLib>
  FaceTriagLagrangeP1LagrangeP1("FaceTriagLagrangeP1LagrangeP1");

  /// Lagrange Triangular Face with P1 geometry and P0 solution.
  GeometricEntityProvider<Face,
          LagrangeShapeFunctionTriagP1,
          LagrangeShapeFunctionTriagP0,
                          ShapeFunctionsLib>
          FaceTriagLagrangeP1LagrangeP0("FaceTriagLagrangeP1LagrangeP0");

//////////////////////////////////////////////////////////////////////////////

    } // namespace ShapeFunctions

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
