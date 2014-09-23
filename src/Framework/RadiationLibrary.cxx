// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/RadiationLibrary.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Framework {
    
//////////////////////////////////////////////////////////////////////////////

void RadiationLibrary::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >
    ("WavelengthMin","Minimum wavelenght to consider.");
  options.addConfigOption< CFreal >
    ("WavelengthMax","Maximum wavelenght to consider.");
  options.addConfigOption< CFuint > 
    ("WavelengthStride","Max nb of spectral points for which radiative properties are computed at once.");
}
 
//////////////////////////////////////////////////////////////////////////////
   
RadiationLibrary::RadiationLibrary(const std::string& name) : 
  Framework::PhysicalPropertyLibrary(name)
{
  addConfigOptionsTo(this);
  
  m_wavMin = -1.0; // default value is on purpose <0
  setParameter("WavelengthMin", &m_wavMin);
  
  m_wavMax = -1.0; // default value is on purpose <0
  setParameter("WavelengthMax", &m_wavMax);
  
  m_wavStride = 1; // default value allows code to run without an actual radiation library
  setParameter("WavelengthStride",&m_wavStride);
}

//////////////////////////////////////////////////////////////////////////////

RadiationLibrary::~RadiationLibrary()
{
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
