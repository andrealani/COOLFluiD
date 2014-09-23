// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Environment/ObjectProvider.hh"

#include "Framework/Framework.hh"
#include "Framework/MaxFilterState.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<MaxFilterState,
			    FilterState,
			    FrameworkLib,1>
maxFilterStateProvider("Max");

//////////////////////////////////////////////////////////////////////////////

MaxFilterState::MaxFilterState(const std::string& name) : 
  FilterState(name)
{
  addConfigOptionsTo(this);
  
  m_minValues = vector<CFreal>();
  Config::ConfigObject::setParameter("minValues",&m_minValues);
}

//////////////////////////////////////////////////////////////////////////////
 
MaxFilterState::~MaxFilterState()
{
}
  
//////////////////////////////////////////////////////////////////////////////

void MaxFilterState::filter (RealVector& state) const
{
  cf_assert( m_maskIDs.size() == m_minValues.size() );
  
  for (CFuint iEq = 0; iEq < state.size(); ++iEq)
  {
    // only apply to the equations selected by the mask
    if ( m_maskIDs[iEq] )
      state[iEq] = max (m_minValues[iEq], state[iEq]);
  }
}

//////////////////////////////////////////////////////////////////////////////

void MaxFilterState::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< vector<CFreal> >
    ("minValues","Minimal values for the variables to filter.");
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
