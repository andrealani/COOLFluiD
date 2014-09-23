// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_NewtonMethod_CrankNichUnSetup_hh
#define COOLFluiD_Numerics_NewtonMethod_CrankNichUnSetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "StdUnSetup.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace NewtonMethod {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a NumericalCommand action to be
/// to be executed in order to deallocate data specific to the
/// method to which it belongs.
/// @author Andrea Lani
class CrankNichUnSetup : public StdUnSetup {
public:

  /// Constructor.
  explicit CrankNichUnSetup(std::string name) :
    StdUnSetup(name)
  {
  }

  /// Destructor.
  ~CrankNichUnSetup()
  {
  }

  /// Execute Processing actions
  void execute();

}; // class Setup

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_NewtonMethod_CrankNichUnSetup_hh

