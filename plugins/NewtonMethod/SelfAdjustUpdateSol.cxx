// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "NewtonMethod/NewtonMethod.hh"
#include "SelfAdjustUpdateSol.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/SubSystemStatus.hh" 
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

MethodCommandProvider<SelfAdjustUpdateSol, NewtonIteratorData, NewtonMethodModule> 
selfAdjustUpdateSolProvider("SelfAdjustUpdateSol");

//////////////////////////////////////////////////////////////////////////////

void SelfAdjustUpdateSol::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

SelfAdjustUpdateSol::SelfAdjustUpdateSol(const std::string& name) : 
  StdUpdateSol(name)
{
}

//////////////////////////////////////////////////////////////////////////////

void SelfAdjustUpdateSol::execute()
{
  CFAUTOTRACE;

  const CFreal oldRes = SubSystemStatusStack::getActive()->getResidual();
  getMethodData().updateResidual();
  const CFreal newRes = SubSystemStatusStack::getActive()->getResidual();
  const CFuint nbIter =  SubSystemStatusStack::getActive()->getNbIter();
  const CFreal factor = 1.5;
  
  bool resDecreased = false;
  if  (newRes > 0. && oldRes > 0.) {
    resDecreased = (newRes < oldRes*factor);
  }
  if  (newRes <= 0. && oldRes > 0.) {
    resDecreased = true;
  }
  if  (newRes >= 0. && oldRes < 0.) {
    resDecreased = (newRes < 1.5);
  }
  if  (newRes < 0. && oldRes < 0.) {
    resDecreased = true;
  }
  
  CFout << "NEW residual = " << newRes << ", OLD residual = " << oldRes  << "\n";
  if (newRes > -30. && (resDecreased || (nbIter == 1)) && newRes < 7.) {
    StdUpdateSol::execute();
  }
  else {
    // the current rhs is set to a value that will surely make the outer iteration loop continue 
    socket_rhs.getDataHandle() = 1000000.;
    
    // reset to 0 the update coefficient 
    socket_updateCoeff.getDataHandle() = 0.0; 
    
    // states are not updated, since the iteration will be repeated
    
    // decrease current number of iterations before repeating the same iteration again
    //  const CFuint nbIter =  SubSystemStatusStack::getActive()->getNbIter();
    //  SubSystemStatusStack::getActive()->setNbIter(nbIter-1);
    
    getMethodData().getConvergenceStatus().iter -= 1;
    
    // decrease the CFL
    const CFreal newCFL = 0.5*getMethodData().getCFL()->getCFLValue();
    getMethodData().getCFL()->setCFLValue(newCFL);
    
    CFout << "Current iteration will be recalculated with CFL = " << newCFL << "\n";	
  }	
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD
