// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/ConvectiveVarSet.hh"
#include "Framework/PhysicalModel.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

ConvectiveVarSet::ConvectiveVarSet(Common::SafePtr<Framework::BaseTerm> term) :
  Common::OwnedObject(),
  Common::NullableObject(),
  _nbEqs(0),
  _entityID(0),
  _iEqSubSys(0),
  _extraData(false),
  _varNames(),
  _maskArray(),
  _hasSourceTerm(false),
  _jacobDissip(0.0),
  _eValuesP(),
  _eValuesM(),
  _fluxArray(),
  _physFlux(),
  _flux(new Flux(*this)),
  _delta(0.0)
{
}

//////////////////////////////////////////////////////////////////////////////

ConvectiveVarSet::~ConvectiveVarSet()
{
}

//////////////////////////////////////////////////////////////////////////////

void ConvectiveVarSet::setup()
{
  cf_assert(_flux.get() != CFNULL);
  
  _nbEqs = PhysicalModelStack::getActive()->getNbEq();
  
  _maskArray.resize(_nbEqs, true);
  _flux->setup();
  
  _eValuesP.resize(_nbEqs);
  _eValuesM.resize(_nbEqs);
}
    
//////////////////////////////////////////////////////////////////////////////

void ConvectiveVarSet::unsetup()
{
  cf_assert(_flux.get() != CFNULL);
  _flux->unsetup();

  _eValuesP.resize(0);
  _eValuesM.resize(0);
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
