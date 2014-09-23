// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "LinearAdv/LinearAdv.hh"
#include "LinearAdvPhysicalModel.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {


  namespace Physics {
namespace LinearAdv {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<LinearAdvPhysicalModel<DIM_1D>, PhysicalModelImpl, LinearAdvModule, 1>
linearAdv1DModelProvider("LinearAdv1D");

Environment::ObjectProvider<LinearAdvPhysicalModel<DIM_2D>, PhysicalModelImpl, LinearAdvModule, 1>
linearAdv2DModelProvider("LinearAdv2D");

Environment::ObjectProvider<LinearAdvPhysicalModel<DIM_3D>, PhysicalModelImpl, LinearAdvModule, 1>
linearAdv3DModelProvider("LinearAdv3D");

//////////////////////////////////////////////////////////////////////////////

} // namespace LinearAdv
  }
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
