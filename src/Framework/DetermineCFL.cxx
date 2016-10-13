// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Environment/ObjectProvider.hh"
#include "Framework/DetermineCFL.hh"
#include "Framework/CFL.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/Framework.hh"
#include "Framework/MeshData.hh"
#include "Framework/EquationSubSysDescriptor.hh"
#include "Framework/PhysicalModel.hh"
#include "Common/BadValueException.hh"
#include "MathTools/MathConsts.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<DetermineCFL,
         ComputeCFL,FrameworkLib,
         1>
determineCFLProvider("Determine");

//////////////////////////////////////////////////////////////////////////////

void DetermineCFL::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::string >("Def","Definition of the Function.");
  options.addConfigOption< CFuint >("SubSystemID","ID of the subsystem to consider for the CFL calculation.");
}
    
//////////////////////////////////////////////////////////////////////////////

DetermineCFL::DetermineCFL(const std::string& name) :
  ComputeCFL(name)
{
  addConfigOptionsTo(this);
  
  _function = "1.0";
  setParameter("Def",&_function);
  
  _subSystemID = 0;
  setParameter("SubSystemID",&_subSystemID);
}

//////////////////////////////////////////////////////////////////////////////

DetermineCFL::~DetermineCFL()
{
}

//////////////////////////////////////////////////////////////////////////////

void DetermineCFL::operator() (const ConvergenceStatus& m_cstatus)
{

  const std::string nsp = MeshDataStack::getActive()->getPrimaryNamespace();
  std::string datahandleName = nsp + "_updateCoeff";
  DataHandle<CFreal> updateCoeff =
    MeshDataStack::getActive()->getDataStorage()->getData<CFreal>(datahandleName);

  datahandleName = nsp + "_volumes";
  DataHandle<CFreal> volumes =
    MeshDataStack::getActive()->getDataStorage()->getData<CFreal>(datahandleName);

  //cf_assert that we are unsteady
  const CFreal timestep = SubSystemStatusStack::getActive()->getDT();
  if (!(timestep > 0.))
    throw Common::BadValueException (FromHere(),"Using DetermineCFL when simulation is not unsteady");

  //compute the value u dt/dx for each cell
  CFreal maxCFL = 0.;
  CFreal minCFL = MathTools::MathConsts::CFrealMax();
  CFreal avgCFL = 0.;
  
  const EquationSubSysDescriptor& eqSSD = PhysicalModelStack::getActive()->
    getEquationSubSysDescriptor();
  
  const CFuint nbLSS = eqSSD.getTotalNbEqSS();
  const CFuint iLSS = _subSystemID;
  CFLog(VERBOSE, "DetermineCFL::operator() => nbLSS = "<< nbLSS << ", iLSS = "<< iLSS << "\n");
  
  const CFuint nbStates = volumes.size();
  for(CFuint iState = 0; iState < nbStates; iState++)
  {
    const CFuint idx = iState*nbLSS + iLSS;
    cf_assert(idx < updateCoeff.size());
    const CFreal stateCFL = updateCoeff[idx]*timestep/volumes[iState];
    
    maxCFL = max(stateCFL, maxCFL);
    minCFL = min(stateCFL, minCFL);
    avgCFL += stateCFL;
  }

  avgCFL /= nbStates;

  _cfl->setCFLValue(maxCFL);
  
  CFLog(VERBOSE, "CFL Statistics - Max: " << maxCFL<< " - Min: " << minCFL<< " - Avg: " << avgCFL<<"\n");
}

//////////////////////////////////////////////////////////////////////////////

void DetermineCFL::configure ( Config::ConfigArgs& args )
{
  ComputeCFL::configure(args);

}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
