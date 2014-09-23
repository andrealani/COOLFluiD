// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_RungeKutta_StdUnSetup_hh
#define COOLFluiD_Numerics_RungeKutta_StdUnSetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "RKData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace RungeKutta {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a NumericalCommand action to be
/// to be executed in order to deallocate data specific to the
/// method to which it belongs.
/// @author Tiago Quintino
class StdUnSetup : public RKCom {
public:

  /// Constructor.
  explicit StdUnSetup(const std::string& name) : RKCom(name)
  {
  }

  /// Destructor.
  ~StdUnSetup()
  {
  }

  /// Execute Processing actions
  void execute();

}; // class UnSetup

//////////////////////////////////////////////////////////////////////////////

    } // namespace RungeKutta

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_RungeKutta_StdUnSetup_hh

