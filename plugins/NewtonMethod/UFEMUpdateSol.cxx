// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "NewtonMethod/NewtonMethod.hh"
#include "UFEMUpdateSol.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"
#include "Common/CFLog.hh"
#include "Framework/State.hh"
#include "Common/BadValueException.hh"

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

MethodCommandProvider<UFEMUpdateSol, NewtonIteratorData, NewtonMethodModule> 
stdUpdateSolProvider("UFEMUpdateSol");

//////////////////////////////////////////////////////////////////////////////

void UFEMUpdateSol::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption<vector<CFreal> >("Relaxation","Relaxation factor");
}

//////////////////////////////////////////////////////////////////////////////

UFEMUpdateSol::UFEMUpdateSol(const std::string& name) : 
  NewtonIteratorCom(name),
  socket_states("states"),
  socket_rhs("rhs"),
  socket_updateCoeff("updateCoeff")
{
   addConfigOptionsTo(this);
   
  _alpha = vector<CFreal>();
  setParameter("Relaxation",&_alpha);
}

//////////////////////////////////////////////////////////////////////////////

void UFEMUpdateSol::setup()
{
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
 
  if (_alpha.size() == 0) {
    _alpha.resize(nbEqs);
    for (CFuint i = 0; i < nbEqs; ++i) {
      _alpha[i] = 1.;
    }
  }

  if (_alpha.size() == 1) {
    const CFreal value = _alpha[0];
    _alpha.resize(nbEqs);
    for (CFuint i = 0; i < nbEqs; ++i) {
      _alpha[i] = value;
    }
  }

  if (_alpha.size() != nbEqs) {
    throw BadValueException (FromHere(),"UFEMUpdateSol::setup() : _alpha.size() != nbEqs");
  }
}

//////////////////////////////////////////////////////////////////////////////

void UFEMUpdateSol::execute()
{
  CFAUTOTRACE;
  
  DataHandle < Framework::State*, Framework::GLOBAL > states  = socket_states.getDataHandle();
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  
  // rhs is the temporary placeholder for the dU
  DataHandle<CFreal>& dU = rhs;

  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint statesSize = states.size();

  SafePtr<FilterState> filterState = getMethodData().getFilterState();
  
  for (CFuint iState = 0; iState < statesSize; ++iState)
  {
    State& cur_state = *states[iState];
    // do the update only if the state is parallel updatable
    if (cur_state.isParUpdatable())
    {
      for (CFuint iEq = 0; iEq < nbEqs; ++iEq)
      {
       cur_state[iEq] += _alpha[iEq] * dU(iState, iEq, nbEqs);
      }
  	  // apply a polymorphic filter to the state
      filterState->filter(cur_state);
    }
  }
  
  // reset to 0 the update coefficient
  updateCoeff = 0.0;
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > UFEMUpdateSol::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result;
  
  result.push_back(&socket_states);
  result.push_back(&socket_rhs);
  result.push_back(&socket_updateCoeff);
  
  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD
