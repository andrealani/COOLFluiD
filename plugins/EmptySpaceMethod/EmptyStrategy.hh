// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_EmptySpaceMethod_EmptyStrategy_hh
#define COOLFluiD_EmptySpaceMethod_EmptyStrategy_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"
#include "EmptySpaceMethod/EmptySolverData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace EmptySpaceMethod {

//////////////////////////////////////////////////////////////////////////////

/// This class represent an empty strategy
/// @author Pedro Maciel
class EmptyStrategy : public EmptySolverStrategy {

public:  // types

  typedef Framework::BaseMethodStrategyProvider<
      EmptySolverData,EmptyStrategy > PROVIDER;

public:  // methods

  /// Constructor
  EmptyStrategy(const std::string& name);

  /// Destructor
  ~EmptyStrategy();

  /// Add compute the term to add in the jacobian
  void compute();

  /// Gets the Class name
  static std::string getClassName()
  {
    return "EmptyStrategy";
  }

  /// Set up private data and data
  virtual void setup();

private: // data

}; // class EmptyStrategy

//////////////////////////////////////////////////////////////////////////////

  }  // namespace EmptySpaceMethod

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_EmptySpaceMethod_EmptyStrategy_hh

