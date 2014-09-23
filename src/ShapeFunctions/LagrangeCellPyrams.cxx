// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/GeometricEntityProvider.hh"

#include "ShapeFunctions/ShapeFunctions.hh"

#include "ShapeFunctions/LagrangeShapeFunctionPyramP1.hh"
#include "ShapeFunctions/LagrangeShapeFunctionPyramP0.hh"

#include "Framework/Cell.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////

  /// Lagrange Pyramid Cell with P1 geometry and P1 solution.
  GeometricEntityProvider<Cell,
                          LagrangeShapeFunctionPyramP1,
                          LagrangeShapeFunctionPyramP1,
                          ShapeFunctionsLib>
  CellPyramLagrangeP1LagrangeP1("CellPyramLagrangeP1LagrangeP1");

  /// Lagrange Pyramid Cell with P1 geometry and P0 solution.
  GeometricEntityProvider<Cell,
            LagrangeShapeFunctionPyramP1,
            LagrangeShapeFunctionPyramP0,
                          ShapeFunctionsLib>
  CellPyramLagrangeP1LagrangeP0("CellPyramLagrangeP1LagrangeP0");

//////////////////////////////////////////////////////////////////////////////

    } // namespace ShapeFunctions
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
