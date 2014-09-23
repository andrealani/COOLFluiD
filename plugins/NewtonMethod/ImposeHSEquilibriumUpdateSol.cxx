// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "NewtonMethod/NewtonMethod.hh"

#include "ImposeHSEquilibriumUpdateSol.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"
#include "Common/CFLog.hh"
#include "Framework/State.hh"
#include "Common/BadValueException.hh"

#include "Framework/SubSystemStatus.hh"
#include "Framework/SpaceMethod.hh"
#include "Framework/SpaceMethodData.hh"
#include "Framework/ConvectiveVarSet.hh"

#include "Framework/ConvergenceMethod.hh"
#include "Framework/ConvergenceMethodData.hh"
#include "Framework/ConvergenceStatus.hh"
#include "NewtonMethod/NewtonIteratorData.hh"
#include "NewtonMethod/NewtonIterator.hh"
#include "NavierStokes/EulerTerm.hh"

#include "Framework/PathAppender.hh"
#include "Framework/GeometricEntity.hh"

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

MethodCommandProvider<ImposeHSEquilibriumUpdateSol, NewtonIteratorData, NewtonMethodModule> 
imposeHSEquilUpdateSolProvider("ImposeHSEquilibriumUpdateSol");

//////////////////////////////////////////////////////////////////////////////

void ImposeHSEquilibriumUpdateSol::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

ImposeHSEquilibriumUpdateSol::ImposeHSEquilibriumUpdateSol(const std::string& name) : 
  StdUpdateSol(name),
  socket_volumes("volumes"),
  socket_gravity("gravity")
{
  addConfigOptionsTo(this);
}

//////////////////////////////////////////////////////////////////////////////

void ImposeHSEquilibriumUpdateSol::setup()
{
  CFAUTOTRACE;
  
  StdUpdateSol::setup();
  
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFuint nbCells = socket_states.getDataHandle().size();
  socket_gravity.getDataHandle().resize(nbCells*dim);
}

//////////////////////////////////////////////////////////////////////////////

void ImposeHSEquilibriumUpdateSol::execute()
{
  CFAUTOTRACE;
  
  const CFuint iter = SubSystemStatusStack::getActive()->getNbIter();
  if (iter == 1 /*&& SubSystemStatusStack::getActive()->isSubIterationLastStep()*/) {
    // for the moment the gravity will be recomputed every subiteration
    
    const CFuint totNbEqs = PhysicalModelStack::getActive()->getNbEq();
    DataHandle<CFreal> rhs = socket_rhs.getDataHandle(); 
    DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
    DataHandle<CFreal> volumes = socket_volumes.getDataHandle(); 
    DataHandle<CFreal> gravity = socket_gravity.getDataHandle();
    DataHandle < Framework::State*, Framework::GLOBAL > states  = socket_states.getDataHandle();
    
    const CFuint nbCells = socket_states.getDataHandle().size();
    const CFuint dim   = PhysicalModelStack::getActive()->getDim();
    const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
    
    SafePtr<SpaceMethod> theSpaceMethod = getMethodData().getCollaborator<SpaceMethod>();
    SafePtr<SpaceMethodData> theSpaceMethodData = theSpaceMethod->getSpaceMethodData();
    SafePtr<ConvectiveVarSet> convVarSet = theSpaceMethodData->getUpdateVar();
    
    RealVector physicalData;
    PhysicalModelStack::getActive()->getImplementor()->getConvectiveTerm()->resizePhysicalData(physicalData);
    cf_assert(physicalData.size() > 0);
    
    for (CFuint i = 0; i < nbCells; ++i) {
      cf_assert(i < volumes.size());
      const CFreal dt = updateCoeff[i]*volumes[i]; // check if this is the actual DT // ask Sarp
      convVarSet->computePhysicalData(*states[i], physicalData);
      const CFreal rho = physicalData[Physics::NavierStokes::EulerTerm::RHO];
      for (CFuint d = 0; d < dim; ++d) {
	cf_assert(i*dim + d < gravity.size());
	cf_assert(i*nbEqs+1+d < rhs.size());
	gravity[i*dim + d] = -rhs[i*nbEqs+1+d]/(dt*rho);
      }
    } 
  }
  
  StdUpdateSol::execute();
}
      
//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > ImposeHSEquilibriumUpdateSol::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result = StdUpdateSol::needsSockets();
  result.push_back(&socket_volumes);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSource> > ImposeHSEquilibriumUpdateSol::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result = 
    StdUpdateSol::providesSockets();
  result.push_back(&socket_gravity);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD
