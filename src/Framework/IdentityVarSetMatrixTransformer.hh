// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_IdentityVarSetMatrixTransformer_hh
#define COOLFluiD_Framework_IdentityVarSetMatrixTransformer_hh

//////////////////////////////////////////////////////////////////////////////

#include "VarSetMatrixTransformer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a transformer of variables from conservative
/// to conservative variables
/// @author Andrea Lani
class Framework_API IdentityVarSetMatrixTransformer : public VarSetMatrixTransformer {
public:

  /// Default constructor without arguments
  IdentityVarSetMatrixTransformer(Common::SafePtr<Framework::PhysicalModelImpl> model);

  /// Default destructor
  ~IdentityVarSetMatrixTransformer();

  /// Checks if this object is a IdentityVarSet object.
  /// Since this is a IdentityVarSetMatrixTransformer it returns true.
  bool getIsIdentityTransformation() const
  {
    return true;
  }

  /// Transform a set of state vectors into another one
  std::vector<Framework::State*>* transform
  (std::vector<Framework::State*> *const  states)
  {
    return states;
  }

  /// Transform a state into another one
  Framework::State* transform(Framework::State* const state)
  {
    return state;
  }

  /// Transform a set of state vectors into another one from
  /// reference precomputed values (physical data) associated to the
  /// given states
  std::vector<Framework::State*>* transformFromRef
  (std::vector<Framework::State*> *const states)
  {
    return states;
  }

  /// Transform a state into another one from reference precomputed
  /// values (physical data)associated to the given state
  Framework::State* transformFromRef(Framework::State* const state)
  {
    return state;
  }

};

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_IdentityVarSetMatrixTransformer_hh
