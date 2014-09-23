// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/GeometricEntityProvider.hh"

#include "ShapeFunctions/ShapeFunctions.hh"

#include "ShapeFunctions/LagrangeShapeFunctionPrismP0.hh"
#include "ShapeFunctions/LagrangeShapeFunctionPrismP1.hh"

#include "Framework/Cell.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////

  /// Lagrange Prism Cell with P1 geometry and P0 solution.
  GeometricEntityProvider<Cell,
            LagrangeShapeFunctionPrismP1,
            LagrangeShapeFunctionPrismP0,
                          ShapeFunctionsLib>
  CellPrismLagrangeP1LagrangeP0("CellPrismLagrangeP1LagrangeP0");

  /// Lagrange Prism Cell with P1 geometry and P1 solution.
  GeometricEntityProvider<Cell,
                          LagrangeShapeFunctionPrismP1,
                          LagrangeShapeFunctionPrismP1,
                          ShapeFunctionsLib>
  CellPrismLagrangeP1LagrangeP1("CellPrismLagrangeP1LagrangeP1");

//////////////////////////////////////////////////////////////////////////////

    } // namespace ShapeFunctions
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
