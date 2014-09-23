// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "SourceVarSet.hh"
#include "PhysicalModel.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

SourceVarSet::SourceVarSet(const std::string& name) :
  Common::OwnedObject(),
  Config::ConfigObject(name),
  Common::NullableObject(),
  _varNames()
{
}

//////////////////////////////////////////////////////////////////////////////

SourceVarSet::~SourceVarSet()
{
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

