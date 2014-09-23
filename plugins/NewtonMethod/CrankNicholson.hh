// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_NewtonMethod_CrankNicholson_hh
#define COOLFluiD_Numerics_NewtonMethod_CrankNicholson_hh

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

/// This class defines a ConvergenceMethod that implements the CrankNicholson
/// method for solving non linear systems
/// @author Andrea Lani
class CrankNicholson : public NewtonIterator {
public:

  /// Default constructor without arguments
  /// @param name missing documentation
  explicit CrankNicholson(const std::string& name);

  /// Default destructor
  ~CrankNicholson();

protected: // abstract interface implementations

  /// Take one timestep
  /// @see ConvergenceMethod::takeStep()
  virtual void takeStepImpl();

}; // class CrankNicholson

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_CrankNicholson_hh
