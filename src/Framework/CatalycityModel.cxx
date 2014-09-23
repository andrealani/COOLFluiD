// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/CatalycityModel.hh"
#include "Framework/PhysicalChemicalLibrary.hh"
#include "Framework/PhysicalModel.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    
//////////////////////////////////////////////////////////////////////////////

void CatalycityModel::defineConfigOptions(Config::OptionList& options)
{
  //   options.addConfigOption< bool >("freezeChemistry","Flag to freeze the chemistry.");
}
    
//////////////////////////////////////////////////////////////////////////////

CatalycityModel::CatalycityModel(const std::string& name)
  : PhysicalPropertyLibrary(name),
    m_library(CFNULL)
{ 
  addConfigOptionsTo(this);
  
  // _freezeChemistry = false;
  //   setParameter("freezeChemistry",&_freezeChemistry);
}
    
//////////////////////////////////////////////////////////////////////////////

CatalycityModel::~CatalycityModel()
{
}

//////////////////////////////////////////////////////////////////////////////

void CatalycityModel::configure ( Config::ConfigArgs& args )
{
  PhysicalPropertyLibrary::configure(args);
}

//////////////////////////////////////////////////////////////////////////////
 
void CatalycityModel::setup()
{
  PhysicalPropertyLibrary::setup();
  
  m_library = PhysicalModelStack::getActive()->
    getImplementor()->getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
}

//////////////////////////////////////////////////////////////////////////////
  
} // namespace Framework
  
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
