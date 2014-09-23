// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "LinearAdv/LinearAdv.hh"
#include "AdvectionDiffusionPhysicalModel.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace Physics {
namespace LinearAdv {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<AdvectionDiffusionPhysicalModel<DIM_2D>, PhysicalModelImpl, LinearAdvModule, 1>
advectiondiffusion2DProvider("AdvectionDiffusion2D");

Environment::ObjectProvider<AdvectionDiffusionPhysicalModel<DIM_3D>, PhysicalModelImpl, LinearAdvModule, 1>
advectiondiffusion3DProvider("AdvectionDiffusion3D");

//////////////////////////////////////////////////////////////////////////////

} // namespace LinearAdv
  }
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

