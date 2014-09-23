// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Environment/ObjectProvider.hh"
#include "Framework/PhysicalModel.hh"

#include "NonLinearAdv/NonLinearAdv.hh"
#include "NonLinearAdv/NonLinearAdvPhysicalModel.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NonLinearAdv {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<NonLinearAdvPhysicalModel<DIM_2D>,
	       PhysicalModelImpl,
 	       NonLinearAdvModule,
	       1>
NonLinearAdv2DProvider("NonLinearAdv2D");


//////////////////////////////////////////////////////////////////////////////

    } // namespace NonLinearAdv

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
