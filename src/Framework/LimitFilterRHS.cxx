// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Environment/ObjectProvider.hh"

#include "Framework/Framework.hh"
#include "Framework/LimitFilterRHS.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<LimitFilterRHS,
			    FilterRHS,
			    FrameworkLib,1>
limitFilterRHSProvider("LimitRHS");

//////////////////////////////////////////////////////////////////////////////

LimitFilterRHS::LimitFilterRHS(const std::string& name) : 
  FilterRHS(name)
{
  addConfigOptionsTo(this);
  
  m_maxDelta = vector<CFreal>();
  Config::ConfigObject::setParameter("maxDelta",&m_maxDelta);
}

//////////////////////////////////////////////////////////////////////////////
 
LimitFilterRHS::~LimitFilterRHS()
{
}
  
//////////////////////////////////////////////////////////////////////////////

void LimitFilterRHS::filter(CFuint iVar, CFreal& dU) 
{
  cf_assert( m_maskIDs.size() == m_maxDelta.size() );
  
  if ( m_maskIDs[iVar] ) {
    dU = (std::abs(dU) < m_maxDelta[iVar]) ? dU : m_maxDelta[iVar];
  }
}

//////////////////////////////////////////////////////////////////////////////

void LimitFilterRHS::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< vector<CFreal> >
    ("maxDelta","Maximum allowable value to impose to the RHS component.");
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
