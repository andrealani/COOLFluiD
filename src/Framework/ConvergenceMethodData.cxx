// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/ConvergenceMethodData.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/ConvergenceMethod.hh"

#ifdef CF_HAVE_SINGLE_EXEC
#include "Framework/IdentityFilterState.hh"
#include "Framework/IdentityFilterRHS.hh"
#endif 

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

void ConvergenceMethodData::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::string >("NormRes","Norm type to monitor residual.");
   options.addConfigOption< std::string >("FilterState","Filter to apply on state variables.");
   options.addConfigOption< std::string >("FilterRHS","Filter to apply on the RHS.");
   options.addConfigOption< CFuint >("SolvingRate","Number of steps after which the convergence method has to be applied.");
   options.addConfigOption< bool >("DoComputeJacobian","Flag to tell to compute the jacobian.");
   options.addConfigOption< bool >("DoUpdateSolution","Flag to tell to update the solution.");
   options.addConfigOption< bool >("FreezeJacobian","Flag to tell to freeze the jacobian during the iterative process.");
   options.addConfigOption< bool >("OnlyPreprocessSolution","Flag to tell to only preprocess the solution once.");
   options.addConfigOption< CFint >("NbLSSToSolveAtOnce","Number of LSS to solve at once in this convergence method.");
}
    
//////////////////////////////////////////////////////////////////////////////

ConvergenceMethodData::ConvergenceMethodData(Common::SafePtr<Method> owner)
  : MethodData(owner), m_cstatus(), m_sm()
{
  addConfigOptionsTo(this);

  m_normStr = "L2";
  setParameter("NormRes",&m_normStr);

  m_filterStateStr = "Identity";
  setParameter("FilterState",&m_filterStateStr);

  m_filterRHSStr = "Identity";
  setParameter("FilterRHS",&m_filterRHSStr);

  m_solvingRate = 1;
  setParameter("SolvingRate",&m_solvingRate);  
  
  m_doComputeJacob = true;
  setParameter("DoComputeJacobian",&m_doComputeJacob);  
  
  m_doUpdateSolution = true;
  setParameter("DoUpdateSolution",&m_doUpdateSolution);
  
  m_freezeJacobian = false;
  setParameter("FreezeJacobian",&m_freezeJacobian);

  m_onlyPreprocessSolution = false;
  setParameter("OnlyPreprocessSolution",&m_onlyPreprocessSolution);
  
  m_nbLSSToSolveAtOnce = -1;
  setParameter("NbLSSToSolveAtOnce",&m_nbLSSToSolveAtOnce);
}

//////////////////////////////////////////////////////////////////////////////

ConvergenceMethodData::~ConvergenceMethodData()
{
  CFLog(VERBOSE, "ConvergenceMethodData::~ConvergenceMethodData()\n");
}

//////////////////////////////////////////////////////////////////////////////

void ConvergenceMethodData::configure ( Config::ConfigArgs& args )
{
  MethodData::configure(args);

  configureNested(&m_CFL, args);

  configureStrategy(args, m_normStr, m_normStr, m_computeNorm);

#ifndef CF_HAVE_SINGLE_EXEC
  configureStrategy(args, m_filterStateStr, m_filterStateStr, m_filterState);
  configureStrategy(args, m_filterRHSStr, m_filterRHSStr, m_filterRHS);
#else
  m_filterState = SelfRegistPtr<FilterState>(new IdentityFilterState("Identity"), CFNULL);
  m_filterRHS   = SelfRegistPtr<FilterRHS>(new IdentityFilterRHS("Identity"), CFNULL);
#endif
}

//////////////////////////////////////////////////////////////////////////////

void ConvergenceMethodData::setup()
{
  MethodData::setup();

  m_computeNorm->setup();

  // set monitored var in SubSystemStatus
  Common::SafePtr<SubSystemStatus> subSysStatus = SubSystemStatusStack::getActive();
  subSysStatus->setMonitoredVar(m_computeNorm->getMonitoredVarIndex());
  subSysStatus->setGlobalRes(m_computeNorm->getGlobalRes());
  
  // set the var registry
  ssys_var_regist = SubSystemStatusStack::getActive()->getVarRegistry();
  
  // add space residual to the var registry
  ssys_var_regist->registVar<RealVector>("spaceResidual", new RealVector(getNormComputer()->getComputedNormVarIDs().size()));  
  
}

//////////////////////////////////////////////////////////////////////////////

void ConvergenceMethodData::unsetup()
{
  MethodData::unsetup();
  
  // unregist en destroy space residual
  RealVector * ptr_spaceResidual = ssys_var_regist->unregistVar<RealVector>("spaceResidual");
  deletePtr(ptr_spaceResidual);
}

//////////////////////////////////////////////////////////////////////////////

RealVector ConvergenceMethodData::getSpaceResidual() const
{
  return ssys_var_regist->getVar<RealVector>("spaceResidual");
}

//////////////////////////////////////////////////////////////////////////////

void ConvergenceMethodData::setSpaceResidual(const RealVector& spaceResidual)
{
  ssys_var_regist->setVar<RealVector>("spaceResidual",spaceResidual);
}


//////////////////////////////////////////////////////////////////////////////

void ConvergenceMethodData::computeSpaceResidualNorm()
{
  /// Warning: This must be called before time contributions are added to the rhs
  if (getOwnMethod().d_castTo<ConvergenceMethod>()->outputSpaceResidual()) {
    setSpaceResidual(m_computeNorm.getPtr()->compute());
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvergenceMethodData::updateResidual()
{
  // set monitored var in SubSystemStatus
  Common::SafePtr<SubSystemStatus> subSysStatus = SubSystemStatusStack::getActive();
  subSysStatus->setResidual( m_computeNorm.getPtr()->compute() );
}

//////////////////////////////////////////////////////////////////////////////

void ConvergenceMethodData::setFactoryRegistry(SafePtr<FactoryRegistry> fr) 
{
  MethodData::setFactoryRegistry(fr);
  m_CFL.setFactoryRegistry(fr);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
