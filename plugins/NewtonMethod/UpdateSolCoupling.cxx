// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "NewtonMethod/NewtonMethod.hh"


#include "UpdateSolCoupling.hh"
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

MethodCommandProvider<UpdateSolCoupling, NewtonIteratorData, NewtonMethodModule> 
updateSolCouplingProvider("UpdateSolCoupling");

//////////////////////////////////////////////////////////////////////////////

UpdateSolCoupling::UpdateSolCoupling(const std::string& name) : 
  StdUpdateSol(name)
{
}

//////////////////////////////////////////////////////////////////////////////

void UpdateSolCoupling::execute()
{
  CFAUTOTRACE;
  
  DataHandle < Framework::State*, Framework::GLOBAL > states  = socket_states.getDataHandle();
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  
  // rhs is the temporary placeholder for the dU
  DataHandle<CFreal>& dU = rhs;
  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
  
  // update only the state components corresponding to the current 
  // equation subsystem
  const EquationSubSysDescriptor& eqSS = 
    PhysicalModelStack::getActive()->getEquationSubSysDescriptor();
  MultiMethodHandle<LinearSystemSolver> lss = getMethodData().getLinearSystemSolver();
  const CFuint currLSS = eqSS.getEqSS();

  // set the data describing the current equation subsystem 
  const vector<CFuint>& currEqs = *lss[currLSS]->getEquationIDs();
  const CFuint nbEqsSS = eqSS.getNbEqsSS();
  const CFuint statesSize = states.size();
  const CFuint totalNbEqs = PhysicalModelStack::getActive()->getNbEq();
  
  for (CFuint iState = 0; iState < statesSize; ++iState) {
    // do the update only if the state is parallel updatable
    if (states[iState]->isParUpdatable()) {
      for (CFuint iEq = 0; iEq < nbEqsSS; ++iEq) {
	const CFuint iVar = currEqs[iEq];
	(*states[iState])[iVar] += m_alpha[iVar]*dU(iState, iVar, totalNbEqs);
      }
    }
    else {
      // reset to 0 the RHS for ghost states in order to avoid 
      // inconsistencies in the parallel L2 norm computation
      cf_assert(!states[iState]->isParUpdatable());
      for (CFuint iEq = 0; iEq < nbEqsSS; ++iEq) {
	const CFuint iVar = currEqs[iEq];
	dU(iState, iVar, totalNbEqs) = 0.;
      }
    }
  }
    
  // reset to 0 the update coefficient
  updateCoeff = 0.0;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD
