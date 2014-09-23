// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/GeometricEntityProvider.hh"

#include "ShapeFunctions/ShapeFunctions.hh"

#include "ShapeFunctions/LagrangeShapeFunctionHexaP0.hh"
#include "ShapeFunctions/LagrangeShapeFunctionHexaP1.hh"
#include "ShapeFunctions/LagrangeShapeFunctionHexaP2.hh"

#include "Framework/Cell.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////

  /// Lagrange Hexadron Cell with P1 geometry and P0 solution.
  GeometricEntityProvider<Cell,
            LagrangeShapeFunctionHexaP1,
            LagrangeShapeFunctionHexaP0,
                          ShapeFunctionsLib>
  CellHexaLagrangeP1LagrangeP0("CellHexaLagrangeP1LagrangeP0");

  /// Lagrange Hexahedron Cell with P1 geometry and P1 solution.
  GeometricEntityProvider<Cell,
                          LagrangeShapeFunctionHexaP1,
                          LagrangeShapeFunctionHexaP1,
                          ShapeFunctionsLib>
  CellHexaLagrangeP1LagrangeP1("CellHexaLagrangeP1LagrangeP1");

  /// Lagrange Hexahedron Cell with P1 geometry and P2 solution.
  GeometricEntityProvider<Cell,
                          LagrangeShapeFunctionHexaP1,
                          LagrangeShapeFunctionHexaP2,
                          ShapeFunctionsLib>
  CellHexaLagrangeP1LagrangeP2("CellHexaLagrangeP1LagrangeP2");

//////////////////////////////////////////////////////////////////////////////

    } // namespace ShapeFunctions
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
