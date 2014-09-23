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
#include "RungeKuttaLS/RungeKuttaLS.hh"
#include "RungeKuttaLS/RungeKuttaStep.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace RungeKuttaLS {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<RungeKuttaStep, RKLSData, RungeKuttaLSModule> RungeKuttaStepProvider("RungeKuttaStep");

//////////////////////////////////////////////////////////////////////////////

RungeKuttaStep::RungeKuttaStep(const std::string& name) :
    RKLSCom(name),
    socket_rhs("rhs"),
    socket_u0("u0"),
    socket_updateCoeff("updateCoeff"),
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
      throw SocketException (FromHere(),"Non essential 'volumes' socket must be plugged for time accurate RungeKuttaLS computation");
}

//////////////////////////////////////////////////////////////////////////////

void RungeKuttaStep::unsetup()
{
}

//////////////////////////////////////////////////////////////////////////////

void RungeKuttaStep::execute()
{
  CFAUTOTRACE;

  // get datahandles of the states, rhs and update coefficients
  DataHandle < Framework::State*, Framework::GLOBAL > states      = socket_states     .getDataHandle();
  DataHandle<RealVector> u0          = socket_u0         .getDataHandle();
  DataHandle<CFreal>     rhs         = socket_rhs        .getDataHandle();
  DataHandle<CFreal>     updateCoeff = socket_updateCoeff.getDataHandle();


  // get current CFL number
  const CFreal cfl = getMethodData().getCFL()->getCFLValue();

  // get boolean telling whether computation is time accurate
  const bool isTimeAccurate = getMethodData().isTimeAccurate();

  // get datahandle of volumes if necessary
  DataHandle<CFreal> volumes(CFNULL);
  if(isTimeAccurate)
  {
    cf_assert(socket_volumes.isConnected());
    volumes = socket_volumes.getDataHandle();
  }

  // get number of equations and number of states
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbStates = states.size();

  // variables for time step and cfl
  CFreal dt  = 0.;
  bool isTimeStepTooLarge = false;
  CFreal maxCFL = 1.;

  // get current stage and scheme coefficients
  const CFuint step  = getMethodData().getCurrentStep();
  const CFreal alpha = getMethodData().getAlpha(step);
  const CFreal beta  = getMethodData().getBeta(step);
  const CFreal oEminusAlpha = 1.0 - alpha;

  // loop over states
  for (CFuint i = 0; i < nbStates; ++i)
  {

    // do the update only if the state is parallel updatable
    if (states[i]->isParUpdatable())
    {

      // compute time step
      if(!isTimeAccurate)
      {
        // Compute pseudo Time Step dt
        if (MathChecks::isZero(updateCoeff[i]))
        {
          dt = 0.;

          CFLogDebugMed("UpdateCoeff[" << i << "] = " << updateCoeff[i] << "\n");

          bool nullSpeed = true;
          for (CFuint j = 0; j < nbEqs; ++j)
          {
            if (MathChecks::isNotZero(rhs(i, j, nbEqs)))
            {
              nullSpeed = false;
            }
          }
          if(!nullSpeed)
          {
            CFLog(ERROR,"Residual contribution to state "<< i << " with null characteristic speed." << "\n");
          }
        }
        else
        {
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

        if(ratio > 1.)
        {
          isTimeStepTooLarge = true;
          maxCFL = max(maxCFL,ratio);
        }
      }

      // update solution
      // Uk+1 = (1.0-alpha[k])*U0 + alpha[k]*Uk + beta[k]*dt*rhs[k]
      dt *= beta;
      for (CFuint j = 0; j < nbEqs; ++j)
      {
        // update of the state
        (*states[i])[j] = oEminusAlpha*u0[i][j] + alpha*(*states[i])[j] + rhs(i,j,nbEqs) * dt;
      }

      cf_assert(states[i]->isValid());
    }
    // reset to 0 the update coefficient
    updateCoeff[i] = 0.0;
  }

  if(isTimeAccurate && isTimeStepTooLarge)
  {
    CFLog(WARN, "The chosen time step is too large as it gives a maximum CFL of " << maxCFL <<".\n");
  }

  if(isTimeAccurate)
  {
    SubSystemStatusStack::getActive()->setMaxDT(SubSystemStatusStack::getActive()->getDT()/maxCFL);
  }

}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > RungeKuttaStep::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_rhs);
  result.push_back(&socket_u0);
  result.push_back(&socket_updateCoeff);
  result.push_back(&socket_states);
  result.push_back(&socket_volumes);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace RungeKuttaLS

  } // namespace Numerics

} // namespace COOLFluiD
