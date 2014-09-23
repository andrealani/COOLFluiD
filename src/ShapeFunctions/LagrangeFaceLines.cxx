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

#include "Framework/Face.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace ShapeFunctions {


//////////////////////////////////////////////////////////////////////////////
// LINE
//////////////////////////////////////////////////////////////////////////////

  /// Lagrange Line Face with P1 geometry and P3 solution.
  GeometricEntityProvider<Face,
                          LagrangeShapeFunctionLineP1,
                          LagrangeShapeFunctionLineP3,
                          ShapeFunctionsLib>
  FaceLineLagrangeP1LagrangeP3("FaceLineLagrangeP1LagrangeP3");

  /// Lagrange Line Face with P1 geometry and P2 solution.
  GeometricEntityProvider<Face,
                          LagrangeShapeFunctionLineP1,
                          LagrangeShapeFunctionLineP2,
                          ShapeFunctionsLib>
  FaceLineLagrangeP1LagrangeP2("FaceLineLagrangeP1LagrangeP2");

  /// Lagrange Line Face with P2 geometry and P1 solution.
  GeometricEntityProvider<Face,
                          LagrangeShapeFunctionLineP2,
                          LagrangeShapeFunctionLineP1,
                          ShapeFunctionsLib>
  FaceLineLagrangeP2LagrangeP1("FaceLineLagrangeP2LagrangeP1");

  /// Lagrange Line Face with P2 geometry and P2 solution.
  GeometricEntityProvider<Face,
                          LagrangeShapeFunctionLineP2,
                          LagrangeShapeFunctionLineP2,
                          ShapeFunctionsLib>
  FaceLineLagrangeP2LagrangeP2("FaceLineLagrangeP2LagrangeP2");

  /// Lagrange Line Face with P3 geometry and P3 solution.
  GeometricEntityProvider<Face,
                          LagrangeShapeFunctionLineP3,
                          LagrangeShapeFunctionLineP3,
                          ShapeFunctionsLib>
  FaceLineLagrangeP3LagrangeP3("FaceLineLagrangeP3LagrangeP3");

  /// Lagrange Line Face with P1 geometry and P1 solution.
  GeometricEntityProvider<Face,
                          LagrangeShapeFunctionLineP1,
                          LagrangeShapeFunctionLineP1,
                          ShapeFunctionsLib>
  FaceLineLagrangeP1LagrangeP1("FaceLineLagrangeP1LagrangeP1");
  
  /// Lagrange Line Face with P1 geometry and P0 solution.
  GeometricEntityProvider<Face,
          LagrangeShapeFunctionLineP1,
          LagrangeShapeFunctionLineP0,
                          ShapeFunctionsLib>
  FaceLineLagrangeP1LagrangeP0("FaceLineLagrangeP1LagrangeP0");
  
//////////////////////////////////////////////////////////////////////////////

    } // namespace ShapeFunctions
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
