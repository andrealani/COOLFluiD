// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_JacobianLinearizer_hh
#define COOLFluiD_Framework_JacobianLinearizer_hh

//////////////////////////////////////////////////////////////////////////////



#include "MathTools/RealVector.hh"
#include "Common/NullableObject.hh"
#include "Environment/ConcreteProvider.hh"
#include "PhysicalModel.hh"
#include "Common/OwnedObject.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class State;

//////////////////////////////////////////////////////////////////////////////

/// This class performs the linearization of the jacobians for the
/// corresponding physical model
/// @author Andrea Lani
/// @author Tiago Quintino
/// @author Thomas Wuilbaut
class Framework_API JacobianLinearizer : public Common::OwnedObject,
  		                     public Common::NullableObject {
public: // functions

  typedef Environment::ConcreteProvider<JacobianLinearizer,1> PROVIDER;
  typedef Common::SafePtr<Framework::PhysicalModel> ARG1;

  /// Default constructor without arguments
  JacobianLinearizer(Common::SafePtr<Framework::PhysicalModel> model);

  /// Default destructor
  virtual ~JacobianLinearizer();

  /// Linearize the states
  /// @param states states to be linearized, passed as pointers to State
  /// @post the result is stored in linearStates
  virtual void linearize(const std::vector<State*>& states) = 0;

  /// Set the variations of extra values
  void setExtraValuesDelta(RealVector *const delta)
  {
    _extraValuesDelta = delta;
  }

  /// Set the update values (they can be useful in case of not exact linearization)
  void setUpdateStates(std::vector<State*> *const upStates)
  {
    _upStates = upStates;
  }

  /// Set the maximum number of states
  void setMaxNbStates(const CFuint& maxNbStates)
  {
    _maxNbStates = maxNbStates;
  }

  /// Get the maximum number of states
  CFuint getMaxNbStates()
  {
    return _maxNbStates;
  }

  /// Gets the Class name
  static std::string getClassName()
  {
    return "JacobianLinearizer";
  }

protected: // data
  
  /// vector to store the temporary sum of linearizing variables
  RealVector        _sumZ;
  
  /// vector to store the temporary average of linearizing variables
  RealVector        _avZ;

  /// max number of states to use for the linearization
  CFuint _maxNbStates;

  /// variations of extra values
  RealVector* _extraValuesDelta;

  /// update values
  std::vector<State*>* _upStates;

  /// cache the physical model pointer
  Common::SafePtr<Framework::PhysicalModel> m_physmodel;

}; // end of class JacobianLinearizer

//////////////////////////////////////////////////////////////////////////////

    } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

Framework_Factory(JacobianLinearizer) // define the factoty instance

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_JacobianLinearizer_hh
