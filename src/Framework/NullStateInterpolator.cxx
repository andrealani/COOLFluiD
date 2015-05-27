// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/Framework.hh"
#include "Framework/NullStateInterpolator.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Environment;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

ObjectProvider<NullStateInterpolator, StateInterpolator, FrameworkLib, 1>
nullInterpolatorProvider("Null");

//////////////////////////////////////////////////////////////////////////////

NullStateInterpolator::NullStateInterpolator(const std::string& name) :
  StateInterpolator(name)
{
  addConfigOptionsTo(this);
}


//////////////////////////////////////////////////////////////////////////////

NullStateInterpolator::~NullStateInterpolator()
{
}

//////////////////////////////////////////////////////////////////////////////

void NullStateInterpolator::setup()
{
  CFAUTOTRACE;
  
  StateInterpolator::setup(); 
  
  CFLog(VERBOSE, "NullStateInterpolator::setup()\n");
}

//////////////////////////////////////////////////////////////////////////////

void NullStateInterpolator::unsetup()
{
  CFAUTOTRACE;
  
  StateInterpolator::unsetup(); 
  
  CFLog(VERBOSE, "NullStateInterpolator::unsetup()\n");
}
    
//////////////////////////////////////////////////////////////////////////////

void NullStateInterpolator::interpolate(const CFuint varID, 
				   const CFreal& in, CFreal& out)
{
  CFLog(VERBOSE, "NullStateInterpolator::interpolate()\n");
}
    
//////////////////////////////////////////////////////////////////////////////

void NullStateInterpolator::interpolate(const RealVector& in, RealVector& out)
{
  CFLog(VERBOSE, "NullStateInterpolator::interpolate()\n");
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
