// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_NewtonMethod_BDF2_hh
#define COOLFluiD_Numerics_NewtonMethod_BDF2_hh

//////////////////////////////////////////////////////////////////////////////

#include "NewtonIterator.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class NumericalCommand;
  }

  namespace Numerics {

    namespace NewtonMethod {

//////////////////////////////////////////////////////////////////////////////

/// This class defines a ConvergenceMethod that implements the BDF2
/// method for solving non linear systems. Backward Euler is
/// used for the first step.
/// @author Thomas Wuilbaut

class BDF2 : public NewtonIterator {
public:

  /// Default constructor without arguments
  /// @param name missing documentation
  explicit BDF2(const std::string& name);

  /// Default destructor
  ~BDF2();

protected: // abstract interface implementations

  /// Take one timestep
  /// @see ConvergenceMethod::takeStep()
  virtual void takeStepImpl();
  
  /// Perform the prepare phase before any iteration
  virtual void prepare ();
  
}; // class BDF2

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_BDF2_hh
