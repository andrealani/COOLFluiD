// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/NullInertiaVarSet.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<NullInertiaVarSet,
               InertiaVarSet,
               FrameworkLib,
               1>
nullInertiaVarSetProvider("Null");

//////////////////////////////////////////////////////////////////////////////

NullInertiaVarSet::NullInertiaVarSet(const std::string& name) :
InertiaVarSet(name)
{
}

//////////////////////////////////////////////////////////////////////////////

NullInertiaVarSet::~NullInertiaVarSet()
{
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
