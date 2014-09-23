// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "NewtonMethod/NewtonMethod.hh"


#include "GReKOUpdateSol.hh"
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

MethodCommandProvider<GReKOUpdateSol, NewtonIteratorData, NewtonMethodModule>
grekoUpdateSolProvider("GReKOUpdateSol");

//////////////////////////////////////////////////////////////////////////////

void GReKOUpdateSol::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption<vector<CFreal> >("Relaxation","Relaxation factor");
  options.addConfigOption<CFreal >("KInlet","Inlet K value");
  options.addConfigOption<CFreal >("OmegaInlet","Inlet Omega value");
  options.addConfigOption<CFreal >("GammaInlet","Inlet Omega value");
  options.addConfigOption<CFreal >("ReThetaInlet","Inlet Omega value");
}

//////////////////////////////////////////////////////////////////////////////

GReKOUpdateSol::GReKOUpdateSol(const std::string& name) :
  NewtonIteratorCom(name),
  socket_states("states"),
  socket_rhs("rhs"),
  socket_updateCoeff("updateCoeff")
{
   addConfigOptionsTo(this);
  _alpha = vector<CFreal>();
   setParameter("Relaxation",&_alpha);

  _KInlet = 1.e-05;
   setParameter("KInlet",&_KInlet);

  _OmegaInlet = 200.;
   setParameter("OmegaInlet",&_OmegaInlet);

  _GammaInlet = 1.;
   setParameter("GammaInlet",&_OmegaInlet);
  
   _ReThetaInlet = 1000.;
   setParameter("ReThetaInlet",&_OmegaInlet);
}

//////////////////////////////////////////////////////////////////////////////

void GReKOUpdateSol::setup()
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
    throw BadValueException (FromHere(),"GReKOUpdateSol::setup() : _alpha.size() != nbEqs");
  }
}

//////////////////////////////////////////////////////////////////////////////

void GReKOUpdateSol::execute()
{
  CFAUTOTRACE;
  
  DataHandle < Framework::State*, Framework::GLOBAL > states  = socket_states.getDataHandle();
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();

  // rhs is the temporary placeholder for the dU
  DataHandle<CFreal>& dU = rhs;

  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbDim = PhysicalModelStack::getActive()->getDim();
  const CFuint statesSize = states.size();

  CFuint negativeK;
  CFuint negativeOmega;
  CFuint negativeGamma;
  CFuint negativeRetheta;
  CFuint negativeP;
  CFuint negativeT;

  negativeP = 0;
  negativeT = 0;
  negativeK = 0;
  negativeOmega = 0;
  negativeGamma = 0;
  negativeRetheta = 0;
  const CFreal Kmin = _KInlet;
  const CFreal OmegaMin = _OmegaInlet;
  const CFreal GammaMin =_GammaInlet;
  const CFreal RethetaMin =  _ReThetaInlet;

  for (CFuint iState = 0; iState < statesSize; ++iState) {
    // do the update only if the state is parallel updatable
    if (states[iState]->isParUpdatable()) {
      for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
        (*states[iState])[iEq] += _alpha[iEq] * dU(iState, iEq, nbEqs);

        if(iEq == 0) {
          if( (*states[iState])[iEq] < 0. ){
            negativeP++;
            (*states[iState])[iEq] = std::max((*states[iState])[iEq], (CFreal)1000.);
          }
        }

        if(iEq == nbDim+1) {
          if( (*states[iState])[iEq] < 0. ){
            negativeT++;
            (*states[iState])[iEq] = std::max((*states[iState])[iEq], (CFreal)298.15);
          }
        }

        if(iEq == nbDim+2) {
          if( (*states[iState])[iEq] < 0. ){
            (*states[iState])[iEq] = std::max((*states[iState])[iEq], Kmin);
            negativeK++;
          }

        }
        if(iEq == nbDim+3) {
          if( (*states[iState])[iEq] < 0. ){
            (*states[iState])[iEq] = std::max((*states[iState])[iEq], OmegaMin);
            negativeOmega++;
          }
        }

        if(iEq == nbDim+4) {
          if( ((*states[iState])[iEq] < 0.)){
            (*states[iState])[iEq] = std::max((*states[iState])[iEq], GammaMin);
            negativeGamma++;
          }
          //if( ((*states[iState])[iEq] > 1.)){
          //  (*states[iState])[iEq] = std::min((*states[iState])[iEq], 0.);
          //  }
        }
        if(iEq == nbDim+5) {
          if( (*states[iState])[iEq] < 0. ){
            (*states[iState])[iEq] = std::max((*states[iState])[iEq], RethetaMin);
            negativeRetheta++;
          }
        }

      }
    }
  }

  CFout << " Negative p appeared :  " << negativeP << " times\n";
  CFout << " Negative T appeared :  " << negativeT << " times\n";
  CFout << " Negative k appeared :  " << negativeK << " times\n";
  CFout << " Negative Omega appeared :  " << negativeOmega << " times\n";
  CFout << " Negative Gamma appeared :  " << negativeGamma << " times\n";
  CFout << " Negative Retheta appeared :  " << negativeRetheta << " times\n";
  // reset to 0 the update coefficient
  updateCoeff = 0.0;
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > GReKOUpdateSol::needsSockets()
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
