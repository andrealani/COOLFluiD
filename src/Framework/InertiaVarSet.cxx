// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "InertiaVarSet.hh"
#include "PhysicalModel.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

InertiaVarSet::InertiaVarSet(const std::string& name) :
  Common::OwnedObject(),
  Config::ConfigObject(name),
  Common::NullableObject(),
  _varNames()
{
}

//////////////////////////////////////////////////////////////////////////////

InertiaVarSet::~InertiaVarSet()
{
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

