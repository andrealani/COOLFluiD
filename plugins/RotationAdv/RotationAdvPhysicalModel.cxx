// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Environment/ObjectProvider.hh"
#include "Framework/PhysicalModel.hh"

#include "RotationAdv/RotationAdv.hh"
#include "RotationAdv/RotationAdvPhysicalModel.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace RotationAdv {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<RotationAdvPhysicalModel<DIM_2D>,
	       PhysicalModelImpl,
 	       RotationAdvModule,
	       1>
rotationAdv2DProvider("RotationAdv2D");

Environment::ObjectProvider<RotationAdvPhysicalModel<DIM_3D>,
	       PhysicalModelImpl,
 	       RotationAdvModule,
	       1>
rotationAdv3DProvider("RotationAdv3D");

//////////////////////////////////////////////////////////////////////////////

    } // namespace RotationAdv

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
