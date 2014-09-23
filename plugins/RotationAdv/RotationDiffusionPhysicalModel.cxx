// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "RotationAdv/RotationAdv.hh"
#include "RotationDiffusionPhysicalModel.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace RotationAdv {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<RotationDiffusionPhysicalModel<DIM_2D>, PhysicalModelImpl, RotationAdvModule, 1>
RotationDiffusion2DProvider("RotationDiffusion2D");

Environment::ObjectProvider<RotationDiffusionPhysicalModel<DIM_3D>, PhysicalModelImpl, RotationAdvModule, 1>
RotationDiffusion3DProvider("RotationDiffusion3D");

//////////////////////////////////////////////////////////////////////////////

    } // namespace RotationAdv

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

