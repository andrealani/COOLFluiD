// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_SetElementStateCoord_hh
#define COOLFluiD_Framework_SetElementStateCoord_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/StringOps.hh"
#include "Common/OwnedObject.hh"
#include "MathTools/RealVector.hh"
#include "Environment/ConcreteProvider.hh"
#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

    class Node;
    class State;

//////////////////////////////////////////////////////////////////////////////

/// This class is a functor and offers an abstract interface
/// for setting the corresponding space coordinates (Node) in
/// the State's in an element.
/// @author Andrea Lani
class Framework_API SetElementStateCoord : public Common::OwnedObject {

public:

  typedef Environment::ConcreteProvider<SetElementStateCoord> PROVIDER;

  /// Constructor
  SetElementStateCoord();

  /// Virtual destructor
  virtual ~SetElementStateCoord()
  {
  }

  /// Overloading of the operator () to make this class act as a
  /// functor
  /// @param nodes   list of the nodes in the current element
  /// @param states  list of the states in the current element
  virtual void operator() (const std::vector<Node*>& nodes,
                           std::vector<State*>& states) = 0;

  /// Function allowing to update the StateCoord
  /// @param nodes   list of the nodes in the current element
  /// @param states  list of the states in the current element
  virtual void update(const std::vector<Node*>& nodes,
                            std::vector<State*>& states) = 0;

  /// Gets the Class name
  static std::string getClassName()
  {
    return "SetElementStateCoord";
  }

protected:

  /// Set the coordinates in the states in isoparametric elements
  /// @param nodes   list of the nodes in the current element
  /// @param states  list of the states in the current element
  /// @pre states.size() == nodes.size()
  void setIsoParamStateCoord(const std::vector<Node*>& nodes,
                             std::vector<State*>& states);

  /// Set the coordinates in the states in P1 (geometry) - P0
  /// (solution) elements
  /// @param nodes   list of the nodes in the current element
  /// @param states  list of the states in the current element
  /// @pre states.size() == 0
  void setLagrangeP1LagrangeP0StateCoord(const std::vector<Node*>& nodes,
                         std::vector<State*>& states);

  /// Update the coordinates in the states in P1 (geometry) - P0
  /// (solution) elements
  /// @param nodes   list of the nodes in the current element
  /// @param states  list of the states in the current element
  /// @pre states.size() == 0
  void updateLagrangeP1LagrangeP0StateCoord(const std::vector<Node*>& nodes,
                         std::vector<State*>& states);

protected:

  /// vector to store temporary coordinates
  RealVector _tempCoord;

}; // end of class SetElementStateCoord

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

Framework_Factory(SetElementStateCoord)

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_SetElementStateCoord_hh
