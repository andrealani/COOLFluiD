// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "RungeKutta2/RungeKutta2.hh"


#include "CorrectorStep.hh"
#include "Framework/State.hh"
#include "Framework/PhysicalModel.hh"
#include "Common/CFLog.hh"
#include "MathTools/MathChecks.hh"
#include "Framework/CFL.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/SubSystemStatus.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace RungeKutta2 {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<CorrectorStep, RK2Data, RungeKutta2Module> correctorStepProvider("CorrectorStep");

//////////////////////////////////////////////////////////////////////////////

void CorrectorStep::execute()
{
  CFAUTOTRACE;
  
  const CFreal CFL = getMethodData().getCFL()->getCFLValue();

  DataHandle < Framework::State*, Framework::GLOBAL > states  = socket_states.getDataHandle();
  DataHandle<RealVector> u0  = socket_u0.getDataHandle();
  DataHandle<CFreal> upCoeff0 = socket_upCoeff0.getDataHandle();
  DataHandle<CFreal> rhs  = socket_rhs.getDataHandle();
  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
  const bool isTimeAccurate = getMethodData().isTimeAccurate();

  DataHandle<CFreal> volumes(CFNULL);
  if(isTimeAccurate){
    volumes = socket_volumes.getDataHandle();
  }

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbStates = states.size();

  bool isTimeStepTooLarge = false;
  CFreal dt = 0.;
  CFreal maxCFL = 1.;

  for (CFuint i = 0; i < nbStates; ++i) {
    // do the update only if the state is parallel updatable
    if (states[i]->isParUpdatable()) {

      if(!isTimeAccurate)
      {
        if (MathChecks::isZero(updateCoeff[i])) {
          dt = 0.;
          bool nullSpeed = true;
          for (CFuint j = 0; j < nbEqs; ++j) {
            if (MathChecks::isNotZero(rhs(i, j, nbEqs))) {
              nullSpeed = false;
            }
          }
          if(!nullSpeed){
            CFLog(ERROR,"Residual contribution to state "
                  << i << " with null characteristic speed." << "\n");
          }
        }
        else {
          dt = CFL / upCoeff0[i];
        }
      }
      else //if time accurate
      {
        // Compute maximum DT
        const CFreal dtmax = 1./updateCoeff[i];
        dt = SubSystemStatusStack::getActive()->getDT() / volumes[i];

        // Compute equivalent CFL
        const CFreal ratio = dt/dtmax;

        if(ratio > 1.){
          isTimeStepTooLarge = true;
          maxCFL = max(maxCFL,ratio);
        }
      }

      (*states[i]) = u0[i];

      for (CFuint j = 0; j < nbEqs; ++j) {
        CFLogDebugMax( "dt = " << dt << "\n");
        // update of the state
        (*states[i])[j] += rhs(i,j,nbEqs) *= dt;
      }
      cf_assert(states[i]->isValid());

    }

    // reset to 0 the update coefficient
    updateCoeff[i] = 0.0;
  }

  if(isTimeAccurate && isTimeStepTooLarge)
      CFLog(WARN, "The chosen time step is too large as it gives a maximum CFL of " << maxCFL <<".\n");

  if(isTimeAccurate)
      SubSystemStatusStack::getActive()->setMaxDT(SubSystemStatusStack::getActive()->getDT()/maxCFL);

 }

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > CorrectorStep::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_rhs);
  result.push_back(&socket_u0);
  result.push_back(&socket_updateCoeff);
  result.push_back(&socket_upCoeff0);
  result.push_back(&socket_states);
  result.push_back(&socket_volumes);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace RungeKutta2

  } // namespace Numerics

} // namespace COOLFluiD
