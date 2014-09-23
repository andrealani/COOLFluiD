// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/PhysicalPropertyLibrary.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

void PhysicalPropertyLibrary::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::string >("path","Library path.");
}

//////////////////////////////////////////////////////////////////////////////

PhysicalPropertyLibrary::PhysicalPropertyLibrary(const std::string& name) :
    Common::NonCopyable<PhysicalPropertyLibrary>(),
    Common::OwnedObject(),
    Common::SetupObject(),
    Common::NullableObject(),
    Config::ConfigObject(name)
{
  addConfigOptionsTo(this);
   
  m_libPath = "";
  setParameter("path",&m_libPath);
}

//////////////////////////////////////////////////////////////////////////////

PhysicalPropertyLibrary::~PhysicalPropertyLibrary() {}

//////////////////////////////////////////////////////////////////////////////

void PhysicalPropertyLibrary::configure ( Config::ConfigArgs& args )
{
  Config::ConfigObject::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
