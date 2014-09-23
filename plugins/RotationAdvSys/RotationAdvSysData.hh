// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Physics_RotationAdvSys_RotationAdvSysData_hh
#define COOLFluiD_Physics_RotationAdvSys_RotationAdvSysData_hh

//////////////////////////////////////////////////////////////////////////////

#include "MathTools/RealVector.hh"
#include "Framework/PhysicalModel.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
namespace RotationAdvSys {

//////////////////////////////////////////////////////////////////////////////

/// This class represents the data for the LinearAdv physical model
/// @author Andrea Lani
class RotationAdvSys_API RotationAdvSysData {

public:

  /// This method has to be called before using the
  /// PhysicalData
  void setup()
  {
    avVel.resize(Framework::PhysicalModelStack::getActive()->getDim());
  }

  RealVector avVel; //average velocity

};

//////////////////////////////////////////////////////////////////////////////

} // namespace LinearAdv
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_LinearAdv_LinearAdvData_hh
