// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_IdentityVarSetTransformer_hh
#define COOLFluiD_Framework_IdentityVarSetTransformer_hh

//////////////////////////////////////////////////////////////////////////////

#include "VarSetTransformer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a transformer of variables from conservative
/// to conservative variables
/// @author Andrea Lani
class Framework_API IdentityVarSetTransformer : public VarSetTransformer {
public:

  /// Default constructor without arguments
  IdentityVarSetTransformer(Common::SafePtr<Framework::PhysicalModelImpl> model);

  /// Default destructor
  ~IdentityVarSetTransformer();

  /// Checks if this object is a IdentityVarSet object.
  /// Since this is a IdentityVarSetTransformer it returns true.
  bool isIdentity() const
  {
    return true;
  }

  /// Transform a set of state vectors into another one
  std::vector<State*>* transform(std::vector<State*> *const  states)
  {
    return states;
  }

  /// Transform a state into another one
  State* transform(State* const state)
  {
    return state;
  }

  /// Transform an array of variables into another one of potentially different length
  void transform(const RealVector& state, RealVector& result)
  {
    copy(state, result);
  }
  
  /// Transform a set of state vectors into another one from
  /// reference precomputed values (physical data) associated to the
  /// given states
  std::vector<State*>* transformFromRef(std::vector<State*> *const  states)
  {
    return states;
  }

  /// Transform a state into another one from reference precomputed
  /// values (physical data)associated to the given state
  State* transformFromRef(State* const state)
  {
    return state;
  }

  /// Transform a state into another one
  void transform(const Framework::State& state, Framework::State& result)
  {
  }

  /// Transform a state into another one from reference precomputed
  /// values (physical data)associated to the given state
  void transformFromRef(const RealVector& data, Framework::State& result)
  {
  }

};

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_IdentityVarSetTransformer_hh
