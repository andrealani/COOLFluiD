// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Burgers/Burgers.hh"
#include "BurgersPhysicalModel.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace Burgers {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<BurgersPhysicalModel<DIM_2D>,
	       PhysicalModelImpl,
	       BurgersModule,
	       1>
burgers2DProvider("Burgers2D");

//////////////////////////////////////////////////////////////////////////////

    } // namespace Burgers

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
