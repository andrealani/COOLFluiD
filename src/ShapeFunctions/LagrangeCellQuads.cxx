// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/GeometricEntityProvider.hh"

#include "ShapeFunctions/ShapeFunctions.hh"

#include "ShapeFunctions/LagrangeShapeFunctionQuadP0.hh"
#include "ShapeFunctions/LagrangeShapeFunctionQuadP1.hh"
#include "ShapeFunctions/LagrangeShapeFunctionQuadP2.hh"

#include "Framework/Cell.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////
// QUAD
//////////////////////////////////////////////////////////////////////////////

  /// Lagrange Quadrilateral Cell with P1 geometry and P0 solution.
  GeometricEntityProvider<Cell,
            LagrangeShapeFunctionQuadP1,
            LagrangeShapeFunctionQuadP0,
                          ShapeFunctionsLib>
    CellQuadLagrangeP1LagrangeP0("CellQuadLagrangeP1LagrangeP0");

  /// Lagrange Quadrilateral Cell with P1 geometry and P1 solution.
  GeometricEntityProvider<Cell,
                          LagrangeShapeFunctionQuadP1,
                          LagrangeShapeFunctionQuadP1,
                          ShapeFunctionsLib>
  CellQuadLagrangeP1LagrangeP1("CellQuadLagrangeP1LagrangeP1");

  /// Lagrange Quad Cell with P1 geometry and P2 solution.
  GeometricEntityProvider<Cell,
                          LagrangeShapeFunctionQuadP1,
                          LagrangeShapeFunctionQuadP2,
                          ShapeFunctionsLib>
  CellQuadLagrangeP1LagrangeP2("CellQuadLagrangeP1LagrangeP2");

//////////////////////////////////////////////////////////////////////////////

    } // namespace ShapeFunctions

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
