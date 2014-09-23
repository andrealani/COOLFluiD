// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_ShapeFunctions_SetHexaLagrangeP1LagrangeP2StateCoord_hh
#define COOLFluiD_ShapeFunctions_SetHexaLagrangeP1LagrangeP2StateCoord_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/SetElementStateCoord.hh"
#include "ShapeFunctions/ShapeFunctionsAPI.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework { class Node;  }
  namespace Framework { class State; }

  namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////

/// This class is a functor and offers an abstract interface
/// for setting the corrisponding space coordinates (Framework::Node) in
/// the State's in a hexahedron with P1 geometrical and P2 solution
/// interpolation.
/// @author Thomas Wuilbaut
class ShapeFunctions_API SetHexaLagrangeP1LagrangeP2StateCoord : public Framework::SetElementStateCoord {

public:

  /// Constructor
  SetHexaLagrangeP1LagrangeP2StateCoord() :
  Framework::SetElementStateCoord(),
  _tempMappedCoords(20),
  _tempShapeFunctions(8)
  {
    const CFuint nbMappedCoordinates = 3;
    for(CFuint i=0;i<_tempMappedCoords.size(); ++i)
    {
      _tempMappedCoords[i].resize(nbMappedCoordinates);
    }
  }


  /// Destructor
  ~SetHexaLagrangeP1LagrangeP2StateCoord()
  {
  }

  /// Overloading of the operator () to make this class act as a
  /// functor
  /// @param nodes   list of the nodes in the current element
  /// @param states  list of the states in the current element
  void operator() (const std::vector<Framework::Node*>& nodes,
                   std::vector<Framework::State*>& states);

  /// Function allowing to update the StateCoord
  /// @param nodes   list of the nodes in the current element
  /// @param states  list of the states in the current element
  void update(const std::vector<Framework::Node*>& nodes,
                            std::vector<Framework::State*>& states);

private:

  ///temp vector for storing the mapped coordinates
  std::vector<RealVector> _tempMappedCoords;

  ///temp vector for storing the shape functions
  RealVector _tempShapeFunctions;

}; // end of class SetHexaLagrangeP1LagrangeP2StateCoord

//////////////////////////////////////////////////////////////////////////////

  } // namespace ShapeFunctions

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_ShapeFunctions_SetHexaLagrangeP1LagrangeP2StateCoord_hh
