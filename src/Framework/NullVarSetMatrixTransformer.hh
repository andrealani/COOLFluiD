// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_NullVarSetMatrixTransformer_hh
#define COOLFluiD_Framework_NullVarSetMatrixTransformer_hh

//////////////////////////////////////////////////////////////////////////////

#include "VarSetMatrixTransformer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a transformer of variables from conservative
/// to conservative variables
/// @author Tiago Quintino
class Framework_API NullVarSetMatrixTransformer : public VarSetMatrixTransformer {
public:

  /// Default constructor without arguments
  NullVarSetMatrixTransformer(Common::SafePtr<Framework::PhysicalModelImpl> model);

  /// Default destructor
  ~NullVarSetMatrixTransformer();

  /// Checks if this object is a NullVarSetMatrix object.
  /// Since this is a NullVarSetMatrixTransformer it returns true.
  bool isNull() const
  {
    return true;
  }

  /// Transform a set of state vectors into another one
  std::vector<State*>* transform(std::vector<State*> *const  states);

  /// Transform a state into another one
  virtual State* transform(State* const state);

private:

  /// Set the flag telling if the transformation is an identity one
  /// @pre this method must be called during set up
  bool getIsIdentityTransformation() const;

}; // end of class NullVarSetMatrixTransformer

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_NullVarSetMatrixTransformer_hh
