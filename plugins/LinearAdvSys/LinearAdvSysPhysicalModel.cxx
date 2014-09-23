// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "LinearAdvSys/LinearAdvSys.hh"
#include "LinearAdvSys/LinearAdvSysPhysicalModel.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace LinearAdvSys {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<LinearAdvSysPhysicalModel<DIM_2D>, PhysicalModelImpl,LinearAdvSysModule, 1>
lineuler2DModelProvider("LinearAdvSys2D");

Environment::ObjectProvider<LinearAdvSysPhysicalModel<DIM_3D>, PhysicalModelImpl,LinearAdvSysModule, 1>
lineuler3DModelProvider("LinearAdvSys3D");

//////////////////////////////////////////////////////////////////////////////

    } // namespace LinearAdvSys

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

