// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "NewtonMethod/NewtonMethod.hh"


#include "CombustionTurbUpdateSol.hh"
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

	/// This class represents a solution updater based on TurbUpdateSol, but
	/// valid for rho_i, u, v, T, k, omega variables instead of P, u, v, T, k, omega
	/// Author: Alessandro Mazzetti
	
MethodCommandProvider<CombustionTurbUpdateSol, NewtonIteratorData, NewtonMethodModule>
combustionTurbUpdateSolProvider("CombustionTurbUpdateSol");

//////////////////////////////////////////////////////////////////////////////

void CombustionTurbUpdateSol::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption<vector<CFreal> >("Relaxation","Relaxation factor");
  options.addConfigOption<CFreal >("KInlet","Inlet K value");
  options.addConfigOption<CFreal >("OmegaInlet","Inlet Omega value");
  options.addConfigOption<CFreal >("TClip","Clipping value for temperature T");
}

//////////////////////////////////////////////////////////////////////////////

CombustionTurbUpdateSol::CombustionTurbUpdateSol(const std::string& name) :
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

  _TClip = 6000.;
   setParameter("TClip",&_TClip);
   //CFout << "Clip Temperature  " << _TClip << "\n";

}

//////////////////////////////////////////////////////////////////////////////

void CombustionTurbUpdateSol::setup()
{
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  
  m_library = PhysicalModelStack::getActive()->getImplementor()->getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
  
  //not needed here
  //const CFuint nbSpecies = m_library->getNbSpecies();
  
  //_ysIn.resize(m_library->getNbSpecies());
  //cf_assert(_ysIn.size() > 0);
  
  negativeRho_i.resize(m_library->getNbSpecies());

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
    throw BadValueException (FromHere(),"CombustionTurbUpdateSol::setup() : _alpha.size() != nbEqs");
  }
}

//////////////////////////////////////////////////////////////////////////////

void CombustionTurbUpdateSol::execute()
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
  
  ////Filter from StdUpdateSol (also lines 145 & 184)
  //SafePtr<FilterState> filterState = getMethodData().getFilterState();
  //SafePtr<FilterRHS> filterRHS     = getMethodData().getFilterRHS();
  
  
  //CFout << "CombustionTurbUpdateSol nbEq " << nbEqs << "\n";
  //CFout << "CombustionTurbUpdateSol nbDim " << nbDim << "\n";
  
  const CFuint nbSpecies = m_library->getNbSpecies();
  CFuint negativeK;
  CFuint negativeOmega;
  //CFuint negativeP;
  CFuint negativeT;

  negativeRho_i = 0.;
  //negativeP = 0;
  negativeT = 0.;
  negativeK = 0.;
  negativeOmega = 0.;
  const CFreal Kmin = _KInlet;
  const CFreal OmegaMin = _OmegaInlet;
  const CFreal ClipT = _TClip;

  for (CFuint iState = 0; iState < statesSize; ++iState) {
    // do the update only if the state is parallel updatable
    if (states[iState]->isParUpdatable()) {
      for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
		
		////Filter from StdUpdateSol  
		//filterRHS->filter(iEq, dU(iState, iEq, nbEqs));
		
        (*states[iState])[iEq] += _alpha[iEq] * dU(iState, iEq, nbEqs);
		//CFout << " states vector " << (*states[iState])[9] << "\n";

        if(iEq <= nbSpecies-1) {
          if( (*states[iState])[iEq] < 0. ){
            negativeRho_i[iEq]++;
			//if rho_i is < 0. clip it to 0.
            (*states[iState])[iEq] = std::max((*states[iState])[iEq], 0.);
			//CFout << " (*states[iState])[iEq] " << (*states[iState])[iEq] << "\n";
          }
        }

        // This works only if there is only one Temperature
        if(iEq == (nbSpecies + nbDim)) {
			if( (*states[iState])[iEq] < 0. ){
            negativeT++;
            (*states[iState])[iEq] = std::max((*states[iState])[iEq], 298.15);
          }
			else if( (*states[iState])[iEq] > ClipT ){
            
            (*states[iState])[iEq] = std::min((*states[iState])[iEq], ClipT);
 	  }
        }
        // This works only if there is only one Temperature
        if(iEq == (nbSpecies + nbDim + 1)) {
          if( (*states[iState])[iEq] < 0. ){
            (*states[iState])[iEq] = std::max((*states[iState])[iEq], Kmin);
            negativeK++;
          }

        }
		// This works only if there is only one Temperature
        if(iEq == (nbSpecies + nbDim + 2)) {
          if( (*states[iState])[iEq] < 0. ){
			//CFout << "*states[iState])[iEq] " << (*states[iState])[iEq] << "\n";
            (*states[iState])[iEq] = std::max((*states[iState])[iEq], OmegaMin);
            negativeOmega++;
          }
        }

      }
	  ////Filter from StdUpdateSol
	  //filterState->filter((*states[iState]));
    }
  }
  
  
  CFout << " Negative rho_i appeared :  " << negativeRho_i << " times\n";
  //if (!negativeT==0) {
  //CFout << " Negative T appeared :  " << negativeT << " times\n";
  //}
  if(!negativeK==0) {
  CFout << " Negative k appeared :  " << negativeK << " times\n";
  }
  if(!negativeOmega==0) {
  CFout << " Negative Omega appeared :  " << negativeOmega << " times\n";
  }
  // reset to 0 the update coefficient
  updateCoeff = 0.;
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > CombustionTurbUpdateSol::needsSockets()
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
