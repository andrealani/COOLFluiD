// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "VarSetTransformer.hh"
#include "PhysicalModel.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

VarSetTransformer::VarSetTransformer(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  Common::OwnedObject(),
  Common::NullableObject(),
  _iState(0),
  _localStateID(0),
  _extraValues(CFNULL),
  _transState(CFNULL),
  _transStateVec()
{
}

//////////////////////////////////////////////////////////////////////////////

VarSetTransformer::~VarSetTransformer()
{
  CFLog(VERBOSE, "VarSetTransformer::~VarSetTransformer() => START\n");
  
  deletePtr(_transState);
  
  CFLog(VERBOSE, "VarSetTransformer::~VarSetTransformer() => size = " << _transStateVec.size()  << "\n");
  
  for (CFuint i = 0; i < _transStateVec.size(); ++i) {
    CFLog(VERBOSE, "VarSetTransformer::~VarSetTransformer() => B[" << i<< "] IN\n");
    deletePtr(_transStateVec[i]); 
    CFLog(VERBOSE, "VarSetTransformer::~VarSetTransformer() => B[" << i<< "] OUT\n");
  }
  CFLog(VERBOSE, "VarSetTransformer::~VarSetTransformer() => END\n");
}
    
//////////////////////////////////////////////////////////////////////////////

void VarSetTransformer::setup(const CFuint maxNbTransStates)
{
  if (_transState == CFNULL) {
    _transState = new State();
  }

  // if the resize function is called on an already setup
  // variable set, clean up the old array of State's
  if (maxNbTransStates > _transStateVec.size()) {
    if (_transStateVec.size() > 0) {
      for (CFuint i = 0; i < _transStateVec.size(); ++i) {
	deletePtr(_transStateVec[i]);
      }
    }

    _transStateVec.resize(maxNbTransStates);
    for (CFuint i = 0; i < _transStateVec.size(); ++i) {
      _transStateVec[i] = new State();
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

std::string VarSetTransformer::getProviderName
(const std::string& physicalModelName,
 const std::string& first,
 const std::string& second)
{
  if (first == second) return "Identity";
  return (physicalModelName + first + "To" + second);
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
