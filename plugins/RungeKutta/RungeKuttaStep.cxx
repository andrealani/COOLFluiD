// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.



#include "Common/CFLog.hh"
#include "MathTools/MathChecks.hh"
#include "Framework/State.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/CFL.hh"
#include "Framework/MeshData.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/SocketException.hh"
#include "RungeKutta/RungeKutta.hh"
#include "RungeKutta/RungeKuttaStep.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace RungeKutta {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<RungeKuttaStep, RKData, RungeKuttaModule> RungeKuttaStepProvider("RungeKuttaStep");

//////////////////////////////////////////////////////////////////////////////

RungeKuttaStep::RungeKuttaStep(const std::string& name) :
    RKCom(name),
    socket_rhs("rhs"),
    socket_u0("u0"),
    socket_updateCoeff("updateCoeff"),
    socket_tempu("tempu"),
    socket_states("states"),
    socket_volumes("volumes",false)
{
}

//////////////////////////////////////////////////////////////////////////////

RungeKuttaStep::~RungeKuttaStep()
{
}

//////////////////////////////////////////////////////////////////////////////

void RungeKuttaStep::setup()
{
  if (!socket_volumes.isConnected() && getMethodData().isTimeAccurate())
      throw SocketException (FromHere(),"Non essential 'volumes' socket must be plugged for time accurate RungeKutta computation");
}

//////////////////////////////////////////////////////////////////////////////

void RungeKuttaStep::unsetup()
{
}

//////////////////////////////////////////////////////////////////////////////

void RungeKuttaStep::execute()
{
  CFAUTOTRACE;
  
  const CFreal cfl = getMethodData().getCFL()->getCFLValue();

  DataHandle < Framework::State*, Framework::GLOBAL > states  = socket_states.getDataHandle();
  DataHandle<RealVector> u0  = socket_u0.getDataHandle();
  DataHandle<RealVector> tempu  = socket_tempu.getDataHandle();
  DataHandle<CFreal> rhs  = socket_rhs.getDataHandle();
  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();

  const bool isTimeAccurate = getMethodData().isTimeAccurate();

  DataHandle<CFreal> volumes(CFNULL);
  if(isTimeAccurate)
  {
    cf_assert(socket_volumes.isConnected());
    volumes = socket_volumes.getDataHandle();
  }

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbStates = states.size();

  CFreal dt  = 0.;
  CFreal dt1 = 0.;
  bool isTimeStepTooLarge = false;
  CFreal maxCFL = 1.;

  const CFuint step  = getMethodData().getCurrentStep();
  const CFreal alpha = getMethodData().getAlpha(step);
  const CFreal beta  = getMethodData().getBeta(step);

  for (CFuint i = 0; i < nbStates; ++i) {

    // do the update only if the state is parallel updatable
    if (states[i]->isParUpdatable()) {

      if(!isTimeAccurate)
      {

        // Compute pseudo Time Step dt
        if (MathChecks::isZero(updateCoeff[i])) {
          dt = 0.;

          CFLogDebugMed("UpdateCoeff[" << i << "] = " << updateCoeff[i] << "\n");

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
          dt = cfl / updateCoeff[i];
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

      dt1 = dt;

      //
      //Compute U(k+1) = U0 + alpha(k) * dt * H(k)
      //
      if (alpha!=0.){
        (*states[i]) = u0[i];
        dt *= alpha ;
        for (CFuint j = 0; j < nbEqs; ++j) {
          CFLogDebugMax( "dt = " << dt << "\n");
          // update of the state
          (*states[i])[j] += rhs(i,j,nbEqs) * dt;
        }
      }

      //
      //Compute Un+1(k+1) = Un+1(k) + beta(k) * dt * H(k)
      //
      if (beta != 0.){
        dt1 *= beta ;
        for (CFuint j = 0; j < nbEqs; ++j) {
          CFLogDebugMax( "dt = " << dt1 << "\n");
          // update of the state
          (tempu[i])[j] += (rhs(i,j,nbEqs) *= dt1) ;
        }
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

vector<SafePtr<BaseDataSocketSink> > RungeKuttaStep::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_rhs);
  result.push_back(&socket_u0);
  result.push_back(&socket_updateCoeff);
  result.push_back(&socket_tempu);
  result.push_back(&socket_states);
  result.push_back(&socket_volumes);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace RungeKutta

  } // namespace Numerics

} // namespace COOLFluiD
