// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_NullInertiaVarSet_hh
#define COOLFluiD_Framework_NullInertiaVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "InertiaVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

    class State;

//////////////////////////////////////////////////////////////////////////////

/// This class represents the basic interface for a Inertia
/// variable set, which provides physical model dependent data
/// and methods associated to a choice of variables
/// @author Andrea Lani
class Framework_API NullInertiaVarSet : public InertiaVarSet {
public:

  /// Constructor
  NullInertiaVarSet(const std::string& name);

  /// Default destructor
  ~NullInertiaVarSet();

/// Gets the Inertia Coeficients
  void getWeakMassMat(const CFreal& W,
                      const CFreal& N,
                      const Framework::State& state,
                      RealMatrix& result)
  {
  }

//////////////////////////////////////////////////////////////////////////////


private: // methods

  /// Private Copy Constructor
  NullInertiaVarSet(const NullInertiaVarSet& v);

  /// Private Assignement operator
  const NullInertiaVarSet& operator=(const NullInertiaVarSet& v);

}; // end of class NullInertiaVarSet

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_NullInertiaVarSet_hh
