// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.



#include "Framework/MeshData.hh"
#include "Framework/State.hh"
#include "Framework/PhysicalModel.hh"
#include "Common/CFLog.hh"

#include "Framework/MethodCommandProvider.hh"

#include "BackwardEuler/BackwardEulerMHD.hh"
#include "BackwardEuler/UpdateSolMHD.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace BackwardEuler {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<UpdateSolMHD, BwdEulerData, BackwardEulerMHDModule> 
updateSolMHDProvider("UpdateSolMHD");

//////////////////////////////////////////////////////////////////////

void UpdateSolMHD::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("pressureCorrectionValue","the correction value for the pressure in negative pressure occurring states, should be a very small positive number.");
}
   

//////////////////////////////////////////////////////////////////////////////

UpdateSolMHD::UpdateSolMHD(const std::string& name) : BwdEulerCom(name),
  socket_rhs("rhs"),
  _model(),
  _pressureCorrectionVal(),
  socket_updateCoeff("updateCoeff"),
  socket_states("states")
{
  addConfigOptionsTo(this);

  _pressureCorrectionVal = MathTools::MathConsts::CFrealEps();
  setParameter("pressureCorrectionValue",&_pressureCorrectionVal);
}

///////////////////////////////////////////////////////////////////////////////

void UpdateSolMHD::setup()
{
  BwdEulerCom::setup();

  _model = PhysicalModelStack::getActive()->getImplementor()->getConvectiveTerm().d_castTo<MHDTerm>();
}

//////////////////////////////////////////////////////////////////////////////

void UpdateSolMHD::execute()
{
  CFAUTOTRACE;
  
  DataHandle < Framework::State*, Framework::GLOBAL > states  = socket_states.getDataHandle();
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();

  // rhs is the temporary placeholder for the dU
  DataHandle<CFreal>& dU = rhs;

  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbStates = states.size();

  CFuint nbStatesWithNegPressure = 0;

  for (CFuint iState = 0; iState < nbStates; ++iState) {
    // do the update only if the state is parallel updatable
    if (states[iState]->isParUpdatable()) {
      for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
        (*states[iState])[iEq] += dU(iState, iEq, nbEqs);
      }
      
      // checking the pressure value in the state and correction in case of negative pressure 
      // (especially necessary for solar wind/planet magnetosphere interaction)
   
      const CFreal gamma = _model->getGamma();
      const CFreal gammaMinus1 = gamma - 1.0;
      const CFreal rho = (*states[iState])[0];
      const CFreal u = (*states[iState])[1]/(*states[iState])[0];
      const CFreal v = (*states[iState])[2]/(*states[iState])[0];
      const CFreal w = (*states[iState])[3]/(*states[iState])[0];
      const CFreal Vsq = u*u + v*v + w*w;
      const CFreal Bsq = (*states[iState])[4]*(*states[iState])[4] +
	      (*states[iState])[5]*(*states[iState])[5] +
	      (*states[iState])[6]*(*states[iState])[6];
      const CFreal rhoE = (*states[iState])[7];
      CFreal p = gammaMinus1*(rhoE-0.5*(rho*Vsq+Bsq));
      
      if (p < 0.0) {
        const CFreal pNegative = p;
	p = _pressureCorrectionVal;
	
	// rhoE is recalculated according to the corrected pressure value
	(*states[iState])[7] = p/gammaMinus1 + 0.5*(rho*Vsq+Bsq);

	nbStatesWithNegPressure += 1;
	
	// position of the negative pressure occurring state is also important
	
	RealVector stateCoord(PhysicalModelStack::getActive()->getDim());
	stateCoord = states[iState]->getCoordinates();
	if (PhysicalModelStack::getActive()->getDim() == 2)
          cout << "Pressure was " << pNegative << " in " << iState << ". state with coordinates (" 
	     << stateCoord[0] << "," << stateCoord[1] << ") and is corrected to " << p << "." << endl;
	else
	  cout << "Pressure was " << pNegative << " in " << iState << ". state with coordinates (" 
	     << stateCoord[0] << "," << stateCoord[1] << "," << stateCoord[2] << ") and is corrected to " << p << "." << endl;
      }
      
    }
  }

  if (nbStatesWithNegPressure > 0)
    cout << "There were " << nbStatesWithNegPressure << " states with negative pressure." << endl;	  

  // reset to 0 the update coefficient
  updateCoeff = 0.0;
}

//////////////////////////////////////////////////////////////////////

void UpdateSolMHD::configure ( Config::ConfigArgs& args )
{
  BwdEulerCom::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
UpdateSolMHD::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_rhs);
  result.push_back(&socket_updateCoeff);
  result.push_back(&socket_states);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace BackwardEuler

  } // namespace Numerics
} // namespace COOLFluiD
