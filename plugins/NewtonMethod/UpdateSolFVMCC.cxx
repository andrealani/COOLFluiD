// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "NewtonMethod/NewtonMethod.hh"
#include "NewtonMethod/UpdateSolFVMCC.hh"

#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"
#include "Common/CFLog.hh"
#include "Framework/State.hh"
#include "Common/BadValueException.hh"

#include "Framework/SubSystemStatus.hh"
#include "Framework/SpaceMethod.hh"
#include "Framework/SpaceMethodData.hh"
#include "Framework/ConvectiveVarSet.hh"

#include "Framework/ConvergenceMethod.hh"
#include "Framework/ConvergenceMethodData.hh"
#include "Framework/ConvergenceStatus.hh"
#include "NewtonMethod/NewtonIteratorData.hh"
#include "NewtonMethod/NewtonIterator.hh"

#include "Framework/PathAppender.hh"
#include "Framework/GeometricEntity.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace NewtonMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<UpdateSolFVMCC, NewtonIteratorData, NewtonMethodModule> 
updateSolFVMCCProvider("UpdateSolFVMCC");

//////////////////////////////////////////////////////////////////////////////

void UpdateSolFVMCC::defineConfigOptions(Config::OptionList& options)
{
  // options.addConfigOption<vector<CFreal>,Config::DynamicOption<> >("Relaxation","Relaxation factor");
}

//////////////////////////////////////////////////////////////////////////////

UpdateSolFVMCC::UpdateSolFVMCC(const std::string& name) : 
  StdUpdateSol(name),
  socket_stencil("stencil"),
  m_correctedState()
{
  addConfigOptionsTo(this);
  
  // by default enable the validation
  m_validate = true;
}
      
//////////////////////////////////////////////////////////////////////////////

void UpdateSolFVMCC::setup()
{
  StdUpdateSol::setup();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  m_correctedState.resize(nbEqs);
}

//////////////////////////////////////////////////////////////////////////////

void UpdateSolFVMCC::correctUnphysicalStates(const std::vector<CFuint>& badStatesIDs)
{
  CFAUTOTRACE;
  
  const CFuint nbBadStates = badStatesIDs.size();
  CFLog(INFO, "UpdateSolFVMCC::correctUnphysicalStates() => [" << nbBadStates << "] states\n");
  
  SafePtr<SpaceMethod> theSpaceMethod = getMethodData().getCollaborator<SpaceMethod>();
  SafePtr<SpaceMethodData> theSpaceMethodData = theSpaceMethod->getSpaceMethodData();
  SafePtr<ConvectiveVarSet> updateVS = theSpaceMethodData->getUpdateVar();
  DataHandle < State*, GLOBAL > states  = socket_states.getDataHandle(); 
  DataHandle < vector<State*> > stencil = socket_stencil.getDataHandle();
    
  for (CFuint iState = 0; iState < nbBadStates; ++iState) {
    const CFuint stateID = badStatesIDs[iState];
    State* const state = states[stateID]; 
    const vector<State*>& s = stencil[stateID];
    const CFuint nbNeighbors = s.size();
    
    CFuint countStates = 0;
    m_correctedState = 0.;
    for (CFuint is = 0; is < nbNeighbors; ++is) {
      cf_assert(nbNeighbors > 0);
      State *const currState = s[is];
      if (state != currState) {
	const bool isGhost = currState->isGhost();	
	const CFuint localID = currState->getLocalID();
	// as a first attempt, only extrapolate from ghost states or neighboring inner states which are not invalid 
	// do not extrapolate from already corrected states
	const bool isValid = (!isGhost) ? 
	  !binary_search(badStatesIDs.begin(), badStatesIDs.end(), localID) : updateVS->isValid(*currState);
	if (isValid) {
	  m_correctedState += (RealVector&)(*currState);
	  countStates++;
	}
      }
    }
    
    if (countStates > 0) {
      // substitute the current state with an average of the sum of all valid neighboring states
      m_correctedState /= countStates;
      state->copyData(m_correctedState);
    }
    else {
      // as a second attempt, if not even one valid neighbor is present, also corrected states 
      // can be used for extrapolation, therefore do not check if neighbors belong to badStatesIDs
      for (CFuint is = 0; is < nbNeighbors; ++is) {
	cf_assert(nbNeighbors > 0);
	State *const currState = s[is];
	if (state != currState) {
	  if (updateVS->isValid(*currState)) {
	    m_correctedState += (RealVector&)(*currState);
	    countStates++;
	  }
	}
      }
      
      if (countStates > 0) {
	// substitute the current state with an average of the sum of all valid neighboring states
	m_correctedState /= countStates;
	state->copyData(m_correctedState);
      }
      else {
	CFLog(ERROR, "UpdateSolFVMCC::correctUnphysicalStates() => no valid states in stencil!\n");
	abort();
      }
    }
  }
}
      
//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > UpdateSolFVMCC::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result = StdUpdateSol::needsSockets();
  result.push_back(&socket_stencil);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD
