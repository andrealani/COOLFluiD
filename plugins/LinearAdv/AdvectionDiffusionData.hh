// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Physics_LinearAdv_AdvectionDiffusionData_hh
#define COOLFluiD_Physics_LinearAdv_AdvectionDiffusionData_hh

//////////////////////////////////////////////////////////////////////////////

#include "MathTools/RealVector.hh"
#include "Framework/PhysicalModel.hh"

#include "LinearAdv/LinearAdv.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace Physics {
namespace LinearAdv {

//////////////////////////////////////////////////////////////////////////////

/// This class represents the data for the NavierStokes physical model
/// @author Andrea Lani
class LinearAdv_API AdvectionDiffusionData {

public:

  /// Overloading of the assignment operator
  const AdvectionDiffusionData& operator=(const AdvectionDiffusionData& data)
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
