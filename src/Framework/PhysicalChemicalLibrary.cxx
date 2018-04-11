// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "PhysicalChemicalLibrary.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    
//////////////////////////////////////////////////////////////////////////////

void PhysicalChemicalLibrary::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFint >("electrEnergyID",
				   "ID of the variables that gets the electronic energy.");
  options.addConfigOption< bool >("freezeChemistry","Flag to freeze the chemistry.");
  options.addConfigOption< bool >("ShiftH0","Shift the formation enthalpy to have H(T=0K)=0.");
  options.addConfigOption< CFdouble, Config::DynamicOption<> >
    ("MaxTe","Maximum value for the electron temperature."); 
}
    
//////////////////////////////////////////////////////////////////////////////

PhysicalChemicalLibrary::PhysicalChemicalLibrary(const std::string& name)
  : PhysicalPropertyLibrary(name),
    _NS(0),
    _NC(0),
    _nbTvib(0),
    _nbTe(0),
    _hasElectrons(false),
    _Rgas(0.),
    _electronPress(0.),
    m_H0(0.),
    m_vecH0(),
    _extraData(),
    _atomicityCoeff(),
    _molecule2EqIDs()
{ 
  addConfigOptionsTo(this);

  _electrEnergyID = -1;
  setParameter("electrEnergyID",&_electrEnergyID);
  
  _freezeChemistry = false;
  setParameter("freezeChemistry",&_freezeChemistry);
  
  m_shiftHO = false;
  setParameter("ShiftH0",&m_shiftHO);
  
  _maxTe = 200000.0;
  setParameter("MaxTe",&_maxTe);
}
    
//////////////////////////////////////////////////////////////////////////////

PhysicalChemicalLibrary::~PhysicalChemicalLibrary()
{
}

//////////////////////////////////////////////////////////////////////////////

void PhysicalChemicalLibrary::configure ( Config::ConfigArgs& args )
{
  PhysicalPropertyLibrary::configure(args);
}

//////////////////////////////////////////////////////////////////////////////
  
} // namespace Framework
  
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
