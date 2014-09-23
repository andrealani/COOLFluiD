// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_NewtonMethod_LinearizedBDF2_hh
#define COOLFluiD_Numerics_NewtonMethod_LinearizedBDF2_hh

//////////////////////////////////////////////////////////////////////////////

#include "BDF2.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class NumericalCommand;
  }

  namespace Numerics {

    namespace NewtonMethod {

//////////////////////////////////////////////////////////////////////////////

/// This class defines a ConvergenceMethod that implements the LinearizedBDF2
/// method for solving non linear systems
/// @author Thomas Wuilbaut

class LinearizedBDF2 : public BDF2 {
public:

  /// Default constructor without arguments
  /// @param name missing documentation
  explicit LinearizedBDF2(const std::string& name);

  /// Default destructor
  ~LinearizedBDF2();


}; // class LinearizedBDF2

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_LinearizedBDF2_hh
