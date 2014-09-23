// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Physics_LinearAdvSys_AdvectionDiffusionSysData_hh
#define COOLFluiD_Physics_LinearAdvSys_AdvectionDiffusionSysData_hh

//////////////////////////////////////////////////////////////////////////////

#include "MathTools/RealVector.hh"
#include "Framework/PhysicalModel.hh"

#include "LinearAdvSys/LinearAdvSys.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
namespace LinearAdvSys {

//////////////////////////////////////////////////////////////////////////////

/// This class represents the data for the NavierStokes physical model
/// @author Andrea Lani
class LinearAdvSys_API AdvectionDiffusionSysData {

public:

  /// Overloading of the assignment operator
  const AdvectionDiffusionSysData& operator=(const AdvectionDiffusionSysData& data)
  {
    nu  = data.nu; // dynamic viscosity
    return *this;
  }

  /// Set up private data
  void setup()  {}

  CFreal nu; // dynamic viscosity

};

//////////////////////////////////////////////////////////////////////////////

} // namespace LinearAdv
  }
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_LinearAdv_AdvectionDiffusion_hh
