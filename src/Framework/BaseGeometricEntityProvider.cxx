// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/BaseGeometricEntityProvider.hh"
#include "Environment/Factory.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

BaseGeometricEntityProvider::BaseGeometricEntityProvider(const std::string& name) 
  : Environment::Provider<GeometricEntity>(name)
{
#ifndef CF_HAVE_SINGLE_EXEC
  Environment::Factory<GeometricEntity>::getInstance().regist(this);
#endif
}

//////////////////////////////////////////////////////////////////////////////

BaseGeometricEntityProvider::~BaseGeometricEntityProvider() 
{
}

//////////////////////////////////////////////////////////////////////////////

std::string BaseGeometricEntityProvider::getProviderType() const
{
  return "GeometricEntity";
}

//////////////////////////////////////////////////////////////////////////////

std::string BaseGeometricEntityProvider::getProviderName() const
{
  return getName();
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework
} // namespace COOLFluiD
