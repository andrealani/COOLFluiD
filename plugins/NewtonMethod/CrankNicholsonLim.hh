// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_NewtonMethod_CrankNicholsonLim_hh
#define COOLFluiD_Numerics_NewtonMethod_CrankNicholsonLim_hh

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

/// This class defines a ConvergenceMethod that implements the CrankNicholsonLim
/// method for solving non linear systems
/// @author Thomas Wuilbaut
class CrankNicholsonLim : public CrankNicholson {
public:

  /// Default constructor without arguments
  /// @param name missing documentation
  explicit CrankNicholsonLim(const std::string& name);

  /// Default destructor
  ~CrankNicholsonLim();

}; // class CrankNicholsonLim

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_CrankNicholsonLim_hh
