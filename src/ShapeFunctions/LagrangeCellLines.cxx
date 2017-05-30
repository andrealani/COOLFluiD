// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/GeometricEntityProvider.hh"

#include "ShapeFunctions/ShapeFunctions.hh"

#include "ShapeFunctions/LagrangeShapeFunctionLineP0.hh"
#include "ShapeFunctions/LagrangeShapeFunctionLineP1.hh"
#include "ShapeFunctions/LagrangeShapeFunctionLineP2.hh"
#include "ShapeFunctions/LagrangeShapeFunctionLineP3.hh"

#include "Framework/Cell.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////
// LINE
//////////////////////////////////////////////////////////////////////////////

  /// Lagrange Line Cell with P1 geometry and P0 solution.
  GeometricEntityProvider<Cell,
            LagrangeShapeFunctionLineP1,
            LagrangeShapeFunctionLineP0,
                          ShapeFunctionsLib>
    CellLineLagrangeP1LagrangeP0("CellLineLagrangeP1LagrangeP0");

    /// Lagrange Line Cell with P1 geometry and P1 solution.
  GeometricEntityProvider<Cell,
                          LagrangeShapeFunctionLineP1,
                          LagrangeShapeFunctionLineP1,
                          ShapeFunctionsLib>
  CellLineLagrangeP1LagrangeP1("CellLineLagrangeP1LagrangeP1");

  /// Lagrange Line Cell with P1 geometry and P2 solution.
  GeometricEntityProvider<Cell,
                          LagrangeShapeFunctionLineP1,
                          LagrangeShapeFunctionLineP2,
                          ShapeFunctionsLib>
  CellLineLagrangeP1LagrangeP2("CellLineLagrangeP1LagrangeP2");

//////////////////////////////////////////////////////////////////////////////

    } // namespace ShapeFunctions

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
