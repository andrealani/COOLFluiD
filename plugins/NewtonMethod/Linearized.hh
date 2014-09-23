// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_NewtonMethod_Linearized_hh
#define COOLFluiD_Numerics_NewtonMethod_Linearized_hh

//////////////////////////////////////////////////////////////////////////////

#include "CrankNicholson.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class NumericalCommand;
  }

  namespace Numerics {

    namespace NewtonMethod {

//////////////////////////////////////////////////////////////////////////////

/// This class defines a ConvergenceMethod that implements the Linearized
/// method for solving non linear systems
/// @author Thomas Wuilbaut

class Linearized : public CrankNicholson {
public:

  /// Default constructor without arguments
  /// @param name missing documentation
  explicit Linearized(const std::string& name);

  /// Default destructor
  ~Linearized();

}; // class Linearized

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_Linearized_hh
