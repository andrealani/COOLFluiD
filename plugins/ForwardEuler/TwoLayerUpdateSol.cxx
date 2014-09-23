// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "ForwardEuler/ForwardEuler.hh"
#include "ForwardEuler/ForwardEuler.hh"
#include "ForwardEuler/ForwardEuler.hh"


#include "TwoLayerUpdateSol.hh"
#include "Framework/State.hh"
#include "Framework/PhysicalModel.hh"
#include "Common/CFLog.hh"
#include "MathTools/MathChecks.hh"
#include "Framework/CFL.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  

    namespace ForwardEuler {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<TwoLayerUpdateSol, FwdEulerData, ForwardEulerLib> TwoLayerUpdateSolProvider("TwoLayerUpdateSol");

//////////////////////////////////////////////////////////////////////////////

void TwoLayerUpdateSol::execute()
{
  CFAUTOTRACE;
  
  const CFreal cfl = getMethodData().getCFL()->getCFLValue();

  DataHandle < Framework::State*, Framework::GLOBAL > states  = socket_states.getDataHandle();
  DataHandle<State*> interStates  = socket_interStates.getDataHandle();
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  DataHandle<CFreal> interRhs = socket_interRhs.getDataHandle();
  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
  DataHandle<CFreal> interUpdateCoeff = socket_interUpdateCoeff.getDataHandle();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbStates = states.size();
  const bool isTimeAccurate = getMethodData().isTimeAccurate();
  if(isTimeAccurate)  CFLog(ERROR, "The option Time Accurate should not be put together with TwoLayer SpaceTime.\n");

  CFreal dt = 0.;
  CFreal interDt = 0.;
  RealVector valueVec(2);

  for (CFuint i = 0; i < nbStates; ++i) {

    // do the update only if the state is parallel updatable
    if (states[i]->isParUpdatable()) {

      CFLogDebugMax( "updateCoeff[i] = " << updateCoeff[i] << "\n");
      CFLogDebugMax( "interUpdateCoeff[i] = " << interUpdateCoeff[i] << "\n");

      if (MathChecks::isZero(updateCoeff[i])) {
        CFLogDebugMax( "updateCoeff[i] = ZERO" << "\n");

        dt = 0.;
        bool nullSpeed = true;
        for (CFuint j = 0; j < nbEqs; ++j) {
          if (MathChecks::isNotZero(rhs(i, j, nbEqs))) {
            nullSpeed = false;
          }
        }
        if (!nullSpeed){
          CFLog(ERROR,"Residual contribution to state "
                << i << " with null characteristic speed.\n");
        }
      }
      else {
        if (MathChecks::isZero(interUpdateCoeff[i])) {
          CFLogDebugMax( "interUpdateCoeff[i] = ZERO" << "\n");

          interDt = 0.;
          bool nullSpeed = true;
          for (CFuint j = 0; j < nbEqs; ++j) {
            if (MathChecks::isNotZero(interRhs(i, j, nbEqs))) {
              nullSpeed = false;
            }
          }
          if (!nullSpeed){
            CFLog(ERROR,"Residual contribution to state "
                  << i << " with null characteristic speed.\n");
          }
        }
        else {
          dt = cfl / updateCoeff[i];
          interDt = cfl / interUpdateCoeff[i];
        }
      }

      for (CFuint j = 0; j < nbEqs; ++j) {
        CFLogDebugMax( "dt = " << dt << "\n");
        CFLogDebugMax( "interDt = " << interDt << "\n");
        // update of the state
        (*states[i])[j] += rhs(i,j,nbEqs) *= dt;
        (*interStates[i])[j] += interRhs(i,j,nbEqs) *= dt;
      }

      cf_assert(states[i]->isValid());
    }

    // computation of the norm of the dU (a.k.a rhs)
    CFreal value = 0.0;
    CFreal interValue = 0.0;
    for (CFuint ii = 0; ii < nbStates; ++ii) {
      const CFreal tmp = rhs(ii, getMethodData().getVarID(), nbEqs);
      const CFreal interTmp = interRhs(ii, getMethodData().getVarID(), nbEqs);
      value += tmp*tmp;
      interValue += interTmp*interTmp;
    }

    valueVec[0] = log10(sqrt(interValue));
    valueVec[1] = log10(sqrt(value));

    getMethodData().setNorm(valueVec);

    updateCoeff[i] = 0.0;
    interUpdateCoeff[i] = 0.0;
  }
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
TwoLayerUpdateSol::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_pastStates);
  result.push_back(&socket_interStates);
  result.push_back(&socket_rhs);
  result.push_back(&socket_interRhs);
  result.push_back(&socket_updateCoeff);
  result.push_back(&socket_interUpdateCoeff);
  result.push_back(&socket_states);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace ForwardEuler



} // namespace COOLFluiD
