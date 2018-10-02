// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Common/NotImplementedException.hh"
#include "Common/BadValueException.hh"
#include "Common/EventHandler.hh"

#include "Environment/CFEnv.hh"

#include "Framework/SpaceMethod.hh"
#include "Framework/SpaceMethodData.hh"
#include "Framework/MeshDataBuilder.hh"
#include "Framework/GlobalJacobianSparsity.hh"
#include "Framework/SubSystemStatus.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

void SpaceMethod::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::string >("Builder","Which MeshDataBuilder should be used with the method.");
   options.addConfigOption< std::string >("JacobianSparsity","Define the Jacobian sparsity that this method produces.");
   options.addConfigOption< bool >("Restart","Option to restart the SpaceMethod from the solution provided.");
}

//////////////////////////////////////////////////////////////////////////////

void SpaceMethod::build_dynamic_functions()
{
  add_dynamic_function("initializeSolution", &SpaceMethod::initializeSolution);
  add_dynamic_function("prepareComputation",&SpaceMethod::prepareComputation);
  add_dynamic_function("applyBC",&SpaceMethod::applyBC);
  add_dynamic_function("postProcessSolution",&SpaceMethod::postProcessSolution);
  add_dynamic_function("preProcessSolution",&SpaceMethod::preProcessSolution);
  add_dynamic_function("extrapolateStatesToNodes",&SpaceMethod::extrapolateStatesToNodes);
}

//////////////////////////////////////////////////////////////////////////////

SpaceMethod::SpaceMethod(const std::string& name)
  : Method(name)
{
  // define which functions might be called dynamic
  build_dynamic_functions();
  // regist which functions might be called by raising Events
  registActionListeners();
  // regist the configuration options
  addConfigOptionsTo(this);
  
  m_restart = false;
  setParameter("Restart",&m_restart);

  m_builder = "";
  setParameter("Builder",&m_builder);

  m_sparsity = "";
  setParameter("JacobianSparsity",&m_sparsity);
}

//////////////////////////////////////////////////////////////////////////////

SpaceMethod::~SpaceMethod()
{
}

//////////////////////////////////////////////////////////////////////////////

void SpaceMethod::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;
    
  Method::configure(args);
  
  m_stored_args = args;
  
  // check that a builder has been chosen
  if (m_builder.empty())
    throw BadValueException (FromHere(),"No builder of mesh data of selected. Set the option 'Builder' in the SpaceMethod");
}

//////////////////////////////////////////////////////////////////////////////

void SpaceMethod::registActionListeners()
{
  CFAUTOTRACE;

  Method::registActionListeners();
  
  Common::SafePtr<EventHandler> event_handler = Environment::CFEnv::getInstance().getEventHandler();
  
  const std::string ssname = SubSystemStatusStack::getCurrentName();   
  event_handler->addListener(event_handler->key(ssname, "CF_ON_MESHADAPTER_BEFOREMESHUPDATE"), 
			     this, &SpaceMethod::beforeMeshUpdateAction);
  event_handler->addListener(event_handler->key(ssname, "CF_ON_MESHADAPTER_AFTERMESHUPDATE"),  
			     this, &SpaceMethod::afterMeshUpdateAction);
  event_handler->addListener(event_handler->key(ssname, "CF_ON_MAESTRO_MODIFYRESTART"),        
			     this, &SpaceMethod::modifyRestartAction);
}
    
//////////////////////////////////////////////////////////////////////////////

Common::Signal::return_t SpaceMethod::modifyRestartAction(Common::Signal::arg_t eModifyRestart)
{
  CFAUTOTRACE;

  m_restart = true;
  
  return Common::Signal::return_t ();
}

//////////////////////////////////////////////////////////////////////////////

void SpaceMethod::initializeSolution()
{
  CFAUTOTRACE;

  CFLog(NOTICE,"-------------------------------------------------------------\n");
  CFLog(NOTICE,"Initializing solution of method [" << getName() << "]\n");
  
  cf_assert(isConfigured());
  cf_assert(isSetup());
  
  pushNamespace();
  
  getSpaceMethodData()->setIsRestart(m_restart); 
  initializeSolutionImpl(m_restart);
  
  popNamespace();
    
  CFLog(NOTICE,"-------------------------------------------------------------\n");
}

//////////////////////////////////////////////////////////////////////////////

void SpaceMethod::prepareComputation()
{
  CFAUTOTRACE;

  cf_assert(isConfigured());
  cf_assert(isSetup());

  pushNamespace();

  prepareComputationImpl();

  popNamespace();
}

//////////////////////////////////////////////////////////////////////////////

void SpaceMethod::computeSpaceResidual(CFreal factor)
{
  CFAUTOTRACE;

  cf_assert(isConfigured());
  cf_assert(isSetup());

  pushNamespace();

  computeSpaceResidualImpl(factor);

  popNamespace();
}

//////////////////////////////////////////////////////////////////////////////

void SpaceMethod::computeTimeResidual(CFreal factor)
{
  CFAUTOTRACE;

  cf_assert(isConfigured());
  cf_assert(isSetup());

  pushNamespace();

  computeTimeResidualImpl(factor);

  popNamespace();
}

//////////////////////////////////////////////////////////////////////////////

void SpaceMethod::applyBC()
{
  CFAUTOTRACE;

  cf_assert(isConfigured());
  cf_assert(isSetup());

  pushNamespace();

  applyBCImpl();

  popNamespace();
}

//////////////////////////////////////////////////////////////////////////////

void SpaceMethod::postProcessSolution()
{
  CFAUTOTRACE;

  cf_assert(isConfigured());
  cf_assert(isSetup());

  pushNamespace();

  postProcessSolutionImpl();

  popNamespace();
}

//////////////////////////////////////////////////////////////////////////////

void SpaceMethod::preProcessSolution()
{
  CFAUTOTRACE;

  cf_assert(isConfigured());
  cf_assert(isSetup());

  pushNamespace();
  
  preProcessSolutionImpl();
  
  popNamespace();
}

//////////////////////////////////////////////////////////////////////////////

// void SpaceMethod::computeSpaceDiagBlockJacobContrib(CFreal factor)
// {
//   CFAUTOTRACE;
//
//   cf_assert(isConfigured());
//   cf_assert(isSetup());
//
//   pushNamespace();
//
//   computeSpaceDiagBlockJacobContribImpl(factor);
//
//   popNamespace();
// }

//////////////////////////////////////////////////////////////////////////////

// void SpaceMethod::computeTimeDiagBlockJacobContrib(CFreal factor)
// {
//   CFAUTOTRACE;
//
//   cf_assert(isConfigured());
//   cf_assert(isSetup());
//
//   pushNamespace();
//
//   computeTimeDiagBlockJacobContribImpl(factor);
//
//   popNamespace();
// }

//////////////////////////////////////////////////////////////////////////////

void SpaceMethod::computeSpaceRhsForStatesSet(CFreal factor)
{
  CFAUTOTRACE;

  cf_assert(isConfigured());
  cf_assert(isSetup());

  pushNamespace();

  computeSpaceRhsForStatesSetImpl(factor);

  popNamespace();
}

//////////////////////////////////////////////////////////////////////////////

void SpaceMethod::computeTimeRhsForStatesSet(CFreal factor)
{
  CFAUTOTRACE;

  cf_assert(isConfigured());
  cf_assert(isSetup());

  pushNamespace();

  computeTimeRhsForStatesSetImpl(factor);

  popNamespace();
}

//////////////////////////////////////////////////////////////////////////////

void SpaceMethod::extrapolateStatesToNodes()
{
  CFAUTOTRACE;

  cf_assert(isConfigured());
  cf_assert(isSetup());

  pushNamespace();

  extrapolateStatesToNodesImpl();

  popNamespace();
}

//////////////////////////////////////////////////////////////////////////////

// void SpaceMethod::computeSpaceDiagBlockJacobContribImpl(CFreal factor)
// {
//   throw Common::NotImplementedException (FromHere(),"computeSpaceDiagBlockJacobContribImpl() not implemented for this SpaceMethod");
// }

//////////////////////////////////////////////////////////////////////////////

// void SpaceMethod::computeTimeDiagBlockJacobContribImpl(CFreal factor)
// {
//   throw Common::NotImplementedException (FromHere(),"computeTimeDiagBlockJacobContribImpl() not implemented for this SpaceMethod");
// }

//////////////////////////////////////////////////////////////////////////////

void SpaceMethod::computeSpaceRhsForStatesSetImpl(CFreal factor)
{
  throw Common::NotImplementedException (FromHere(),"computeSpaceRhsForStatesSetImpl() not implemented for this SpaceMethod");
}

//////////////////////////////////////////////////////////////////////////////

void SpaceMethod::computeTimeRhsForStatesSetImpl(CFreal factor)
{
  throw Common::NotImplementedException (FromHere(),"computeTimeRhsForStatesSetImpl() not implemented for this SpaceMethod");
}

//////////////////////////////////////////////////////////////////////////////

Common::SelfRegistPtr<MeshDataBuilder> SpaceMethod::createMeshBuilder()
{
  CFAUTOTRACE;

  cf_assert(isConfigured());

  // allocate and configure the method mesh data builder

  Common::SelfRegistPtr<MeshDataBuilder> ptr =
    FACTORY_GET_PROVIDER(getFactoryRegistry(), MeshDataBuilder, m_builder)->create(m_builder);
  ptr->setFactoryRegistry(getFactoryRegistry());
  
  configureNested ( ptr.getPtr(), m_stored_args );
  return ptr;
}
    
//////////////////////////////////////////////////////////////////////////////

Common::SelfRegistPtr<GlobalJacobianSparsity> SpaceMethod::createJacobianSparsity()
{
  CFAUTOTRACE;

  cf_assert(isConfigured());
  
  Common::SelfRegistPtr<GlobalJacobianSparsity> ptr =
    FACTORY_GET_PROVIDER(getFactoryRegistry(), GlobalJacobianSparsity, m_sparsity)->create();
  ptr->setFactoryRegistry(getFactoryRegistry());
  
  return ptr;
}

//////////////////////////////////////////////////////////////////////////////

Common::Signal::return_t SpaceMethod::beforeMeshUpdateAction(Common::Signal::arg_t eBefore)
{
  CFAUTOTRACE;

  cf_assert(isConfigured());
  cf_assert(isSetup());

  pushNamespace();

  Common::Signal::return_t ret = beforeMeshUpdateActionImpl(eBefore);

  popNamespace();

  return ret;
}

//////////////////////////////////////////////////////////////////////////////

Common::Signal::return_t SpaceMethod::afterMeshUpdateAction(Common::Signal::arg_t eAfter)
{
  CFAUTOTRACE;

  cf_assert(isConfigured());
  cf_assert(isSetup());

  pushNamespace();

  Common::Signal::return_t ret = afterMeshUpdateActionImpl(eAfter);

  popNamespace();

  return ret;
}

//////////////////////////////////////////////////////////////////////////////

void SpaceMethod::setComputeJacobianFlag(bool flag)
{
  getSpaceMethodData()->setComputeJacobianFlag(flag);
}

//////////////////////////////////////////////////////////////////////////////

void SpaceMethod::setOnlyPreprocessSolution(bool flag)
{
  getSpaceMethodData()->setOnlyPreprocessSolution(flag);
}

//////////////////////////////////////////////////////////////////////////////

void SpaceMethod::setMethodImpl()
{
}

//////////////////////////////////////////////////////////////////////////////

void SpaceMethod::unsetMethodImpl()
{
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
