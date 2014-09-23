// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "LinearAdvSys/LinearAdvSys.hh"
#include "AdvectionDiffusionSysPhysicalModel.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
namespace LinearAdvSys {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<AdvectionDiffusionSysPhysicalModel<DIM_2D>, PhysicalModelImpl, LinearAdvSysModule, 1>
advectiondiffusionSys2DProvider("AdvectionDiffusionSys2D");

Environment::ObjectProvider<AdvectionDiffusionSysPhysicalModel<DIM_3D>, PhysicalModelImpl, LinearAdvSysModule, 1>
advectiondiffusionSys3DProvider("AdvectionDiffusionSys3D");

//////////////////////////////////////////////////////////////////////////////

} // namespace LinearAdv
  }
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

