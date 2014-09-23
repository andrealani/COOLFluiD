// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Environment/ObjectProvider.hh"
#include "Framework/PhysicalModel.hh"

#include "RotationAdvSys/RotationAdvSys.hh"
#include "RotationAdvSys/RotationAdvSysPhysicalModel.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace RotationAdvSys {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<RotationAdvSysPhysicalModel<DIM_2D>,
	       PhysicalModelImpl,
 	       RotationAdvSysModule,
	       1>
rotationAdvSys2DProvider("RotationAdvSys2D");

Environment::ObjectProvider<RotationAdvSysPhysicalModel<DIM_3D>,
	       PhysicalModelImpl,
 	       RotationAdvSysModule,
	       1>
rotationAdvSys3DProvider("RotationAdvSys3D");

//////////////////////////////////////////////////////////////////////////////

    } // namespace RotationAdv

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
