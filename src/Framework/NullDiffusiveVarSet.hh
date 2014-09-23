// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_NullDiffusiveVarSet_hh
#define COOLFluiD_Framework_NullDiffusiveVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "DiffusiveVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

    class State;

//////////////////////////////////////////////////////////////////////////////

/// This class represents the basic interface for a diffusive
/// variable set, which provides physical model dependent data
/// and methods associated to a choice of variables
/// @author Andrea Lani
/// @author Tiago Quintino
class Framework_API NullDiffusiveVarSet : public DiffusiveVarSet {
public:

  /// Constructor
  NullDiffusiveVarSet(const std::string& name, Common::SafePtr<Framework::PhysicalModelImpl> model);

  /// Default destructor
  ~NullDiffusiveVarSet();

  /// Get the diffusive flux
  RealVector& getFlux(const RealVector& values,
                      const std::vector<RealVector*>& gradients,
                      const RealVector& normal,
                      const CFreal& radius);

  /// Get the diffusive flux vector
  RealMatrix& getFlux(const RealVector& values,
                      const std::vector<RealVector*>& gradients,
                      const CFreal& radius);

private: // methods

  /// Private Copy Constructor
  NullDiffusiveVarSet(const NullDiffusiveVarSet& v);

  /// Private Assignement operator
  const NullDiffusiveVarSet& operator=(const NullDiffusiveVarSet& v);

}; // end of class NullDiffusiveVarSet

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_NullDiffusiveVarSet_hh
