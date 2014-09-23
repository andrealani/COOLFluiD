// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "NewtonMethod/NewtonMethod.hh"


#include "NewmarkExplicitUpdateSol.hh"
#include "Framework/MeshData.hh"
#include "Framework/State.hh"
#include "Framework/PhysicalModel.hh"
#include "Common/CFLog.hh"

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

MethodCommandProvider<NewmarkExplicitUpdateSol, NewtonIteratorData, NewtonMethodModule> newmarkExplicitUpdateSolProvider("NewmarkExplicitUpdateSol");

//////////////////////////////////////////////////////////////////////////////

void NewmarkExplicitUpdateSol::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("Alpha","Alpha coef");
   options.addConfigOption< CFreal >("Gamma","Gamma coef");
}

//////////////////////////////////////////////////////////////////////////////

NewmarkExplicitUpdateSol::NewmarkExplicitUpdateSol(const std::string& name)
  : NewtonIteratorCom(name),
    socket_rhs("rhs"),
    socket_updateCoeff("updateCoeff"),
    socket_pastStates("pastStates"),
    socket_pastStatesD("pastStatesD"),
    socket_pastStatesD2("pastStatesD2"),
    socket_states("states")
{
   addConfigOptionsTo(this);
  _alpha = 0.5;
   setParameter("Alpha",&_alpha);

 _gamma = 0.5;
   setParameter("Gamma",&_gamma);
}

//////////////////////////////////////////////////////////////////////////////

void NewmarkExplicitUpdateSol::execute()
{
  CFAUTOTRACE;
  
  DataHandle < Framework::State*, Framework::GLOBAL > states  = socket_states.getDataHandle();
  DataHandle<State*> pastStates  = socket_pastStates.getDataHandle();
  DataHandle<State*> pastStatesD  = socket_pastStatesD.getDataHandle();
  DataHandle<State*> pastStatesD2 = socket_pastStatesD2.getDataHandle();
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint statesSize = states.size();

  for (CFuint iState = 0; iState < statesSize; ++iState) {
    // do the update only if the state is parallel updatable
    if (states[iState]->isParUpdatable()) {
      for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
         (*states[iState])[iEq] = rhs(iState, iEq, nbEqs);
      }
    }
  }

//   if(SubSystemStatusStack::getActive()->isSubIterationLastStep())
//   {
//     // Compute the coefficients
//     const CFreal dt = SubSystemStatusStack::getActive()->getDT();
//     const CFreal a1 = _alpha * dt;
//     const CFreal a2 = (1.- _alpha) * dt;
//     const CFreal a3 = 2./(_gamma * dt * dt);
//     const CFreal a4 = dt * a3;
//     const CFreal a5 = (1./_gamma) - 1.;
//
//     // Compute the first and second derivatives
//     State statesD;
//     State statesD2;
//
//     for (CFuint iState = 0; iState < statesSize; ++iState) {
//       // do the update only if the state is parallel updatable
//       if (states[iState]->isParUpdatable()) {
//         for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
//           statesD2[iEq] =    a3 * ((*states[iState])[iEq] - (*pastStates[iState])[iEq])
//                           - (a4 * (*pastStatesD[iState])[iEq])
//                           - (a5 * (*pastStatesD2[iState])[iEq]);
//           statesD[iEq]  =         (*pastStatesD[iState])[iEq]
//                           + (a2 * (*pastStatesD2[iState])[iEq])
//                           + (a1 * statesD2[iEq]);
//           }
//
//         *pastStatesD[iState] = statesD;
//         *pastStatesD2[iState] = statesD2;
//       }
//     }
//   }

  // reset to 0 the update coefficient
  updateCoeff = 0.0;
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > NewmarkExplicitUpdateSol::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_rhs);
  result.push_back(&socket_updateCoeff);
  result.push_back(&socket_pastStates);
  result.push_back(&socket_pastStatesD);
  result.push_back(&socket_pastStatesD2);
  result.push_back(&socket_states);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD
