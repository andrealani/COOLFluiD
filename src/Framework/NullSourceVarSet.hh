// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_NullSourceVarSet_hh
#define COOLFluiD_Framework_NullSourceVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "SourceVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

    class State;

//////////////////////////////////////////////////////////////////////////////

/// This class represents the basic interface for a Source
/// variable set, which provides physical model dependent data
/// and methods associated to a choice of variables
/// @author Tiago Quintino
/// @author Andrea Lani
class Framework_API NullSourceVarSet : public SourceVarSet {
public:

  /// Constructor
  NullSourceVarSet(const std::string& name);

  /// Default destructor
  ~NullSourceVarSet();

  /// Checks if the source term has a part independent from the solution
  bool hasIndepCoef() const
  {
    return false;
  }

  /// Checks if the source term has a linear part
  bool hasLinearCoef() const
  {
    return false;
  }

  /// Gets the Linear Source Coeficients
  void getLinearSourceCoefs(const Framework::State& state, const RealVector& normals, RealMatrix& coef);

  /// Gets the IndepSource Coeficients
  void getIndepSourceCoefs(const Framework::State& state, const RealVector& normals, RealVector& coef);

  /// Gets the Linear Source Coeficients
  void getLinearSourceCoefs(const Framework::State& state, RealMatrix& coef);

  /// Gets the IndepSource Coeficients
  void getIndepSourceCoefs(const Framework::State& state, RealVector& coef);

private: // methods

  /// Private Copy Constructor
  NullSourceVarSet(const NullSourceVarSet& v);

  /// Private Assignement operator
  const NullSourceVarSet& operator=(const NullSourceVarSet& v);

}; // end of class NullSourceVarSet

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_NullSourceVarSet_hh
