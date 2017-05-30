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

#include "Framework/Cell.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////
// TRIANGLE
//////////////////////////////////////////////////////////////////////////////

  /// Lagrange Triangular Cell with P1 geometry and P0 solution.
  GeometricEntityProvider<Cell,
            LagrangeShapeFunctionTriagP1,
            LagrangeShapeFunctionTriagP0,
            ShapeFunctionsLib>
    CellTriagLagrangeP1LagrangeP0("CellTriagLagrangeP1LagrangeP0");

  /// Lagrange Triangular Cell with P1 geometry and P1 solution.
  GeometricEntityProvider<Cell,
                          LagrangeShapeFunctionTriagP1,
                          LagrangeShapeFunctionTriagP1,
                          ShapeFunctionsLib>
  CellTriagLagrangeP1LagrangeP1("CellTriagLagrangeP1LagrangeP1");

  /// Lagrange Triangular Cell with P1 geometry and P2 solution.
  GeometricEntityProvider<Cell,
                          LagrangeShapeFunctionTriagP1,
                          LagrangeShapeFunctionTriagP2,
                          ShapeFunctionsLib>
  CellTriagLagrangeP1LagrangeP2("CellTriagLagrangeP1LagrangeP2");

  /// Lagrange Triangle Cell with P1 geometry and P3 solution.
  GeometricEntityProvider<Cell,
                          LagrangeShapeFunctionTriagP1,
                          LagrangeShapeFunctionTriagP3,
                          ShapeFunctionsLib>
  CellTriagLagrangeP1LagrangeP3("CellTriagLagrangeP1LagrangeP3");

  /// Lagrange Triangular Cell with P2 geometry and P2 solution.
  GeometricEntityProvider<Cell,
  		 LagrangeShapeFunctionTriagP2,
  		 LagrangeShapeFunctionTriagP2,
                          ShapeFunctionsLib>
  CellTriagLagrangeP2LagrangeP2("CellTriagLagrangeP2LagrangeP2");

  /// Lagrange Triangular Cell with P2 geometry and P1 solution.
  GeometricEntityProvider<Cell,
  		 LagrangeShapeFunctionTriagP2,
  		 LagrangeShapeFunctionTriagP1,
                          ShapeFunctionsLib>
  CellTriagLagrangeP2LagrangeP1("CellTriagLagrangeP2LagrangeP1");

  /// Lagrange Triangular Cell with P3 geometry and P3 solution.
  GeometricEntityProvider<Cell,
       LagrangeShapeFunctionTriagP3,
       LagrangeShapeFunctionTriagP3,
                          ShapeFunctionsLib>
  CellTriagLagrangeP3LagrangeP3("CellTriagLagrangeP3LagrangeP3");

//////////////////////////////////////////////////////////////////////////////

    } // namespace ShapeFunctions
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
