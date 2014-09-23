// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/MethodData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

void MethodData::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::vector<std::string> >("CollaboratorNames","Names of the collaborating Methods.");
}

//////////////////////////////////////////////////////////////////////////////

MethodData::MethodData(Common::SafePtr<Method> owner) :
  OwnedObject(),
  ConfigObject("Data"),
  CollaboratorAccess(owner),
  m_strategies()
{
  addConfigOptionsTo(this);

  setParameter("CollaboratorNames",&m_CollaboratorNames);
}

//////////////////////////////////////////////////////////////////////////////

MethodData::~MethodData()
{
}

//////////////////////////////////////////////////////////////////////////////

void MethodData::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  ConfigObject::configure(args);
  CollaboratorAccess::process();
}

//////////////////////////////////////////////////////////////////////////////

void MethodData::setup()
{
  CFAUTOTRACE;

  SetupObject::setup();
  CollaboratorAccess::setup();
}

//////////////////////////////////////////////////////////////////////////////

void MethodData::unsetup()
{
  CFAUTOTRACE;

  CollaboratorAccess::setup(); 
  SetupObject::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD
