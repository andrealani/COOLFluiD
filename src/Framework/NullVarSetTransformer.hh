// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_NullVarSetTransformer_hh
#define COOLFluiD_Framework_NullVarSetTransformer_hh

//////////////////////////////////////////////////////////////////////////////

#include "VarSetTransformer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a transformer of variables from conservative
/// to conservative variables
/// @author Tiago Quintino
/// @author Andrea Lani
class Framework_API NullVarSetTransformer : public VarSetTransformer {
public:

  /// Default constructor without arguments
  NullVarSetTransformer(Common::SafePtr<Framework::PhysicalModelImpl> model);

  /// Default destructor
  ~NullVarSetTransformer();

  /// Checks if this object is a NullVarSet object.
  /// Since this is a NullVarSetTransformer it returns true.
  bool isNull() const { return true; }

  /// Transform a state into another one
  void transform(const Framework::State& state, Framework::State& result);

  /// Transform a state into another one from reference precomputed
  /// values (physical data)associated to the given state
  void transformFromRef(const RealVector& pdata, Framework::State& result);

}; // end of class NullVarSetTransformer

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_NullVarSetTransformer_hh
