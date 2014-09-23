// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/GeometricEntityProvider.hh"

#include "ShapeFunctions/ShapeFunctions.hh"

#include "ShapeFunctions/LagrangeShapeFunctionTetraP0.hh"
#include "ShapeFunctions/LagrangeShapeFunctionTetraP1.hh"
#include "ShapeFunctions/LagrangeShapeFunctionTetraP2.hh"

#include "Framework/Cell.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////
// TETRA
//////////////////////////////////////////////////////////////////////////////

  /// Lagrange Tetrahedral Cell with P1 geometry and P0 solution.
  GeometricEntityProvider<Cell,
            LagrangeShapeFunctionTetraP1,
            LagrangeShapeFunctionTetraP0,
                          ShapeFunctionsLib>
  CellTetraLagrangeP1LagrangeP0("CellTetraLagrangeP1LagrangeP0");

  /// Lagrange Tetrahedral Cell with P1 geometry and P1 solution.
  GeometricEntityProvider<Cell,
                          LagrangeShapeFunctionTetraP1,
                          LagrangeShapeFunctionTetraP1,
                          ShapeFunctionsLib>
  CellTetraLagrangeP1LagrangeP1("CellTetraLagrangeP1LagrangeP1");

  /// Lagrange Tetrahedral Cell with P1 geometry and P2 solution.
  GeometricEntityProvider<Cell,
                          LagrangeShapeFunctionTetraP1,
                          LagrangeShapeFunctionTetraP2,
                          ShapeFunctionsLib>
  CellTetraLagrangeP1LagrangeP2("CellTetraLagrangeP1LagrangeP2");

    GeometricEntityProvider<Cell,
			    LagrangeShapeFunctionTetraP2,
			    LagrangeShapeFunctionTetraP2,
			    ShapeFunctionsLib>
    CellTetraLagrangeP2LagrangeP2("CellTetraLagrangeP2LagrangeP2");

//////////////////////////////////////////////////////////////////////////////

    } // namespace ShapeFunctions
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
