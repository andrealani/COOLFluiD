// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Environment/ObjectProvider.hh"
#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/FluxReconstructionSolver.hh"
#include "FluxReconstructionMethod/ReconstructStatesFluxReconstruction.hh"
#include "FluxReconstructionMethod/BasePointDistribution.hh"
#include "FluxReconstructionMethod/BaseInterfaceFlux.hh"
#include "FluxReconstructionMethod/ConvBndFaceTermRHSFluxReconstruction.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider< FluxReconstructionSolver,SpaceMethod,FluxReconstructionModule,1 >
  fluxReconstructionSolverProvider("FluxReconstruction");

//////////////////////////////////////////////////////////////////////////////

void FluxReconstructionSolver::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::string >("SolveCom","Command to solve the problem with FluxReconstruction solver.");
  options.addConfigOption< std::string >("UnSetupCom","Command to deallocate FluxReconstruction solver data.");
  options.addConfigOption< std::string >("SetupCom","Command to initialize FluxReconstruction solver data.");
  options.addConfigOption< std::vector<std::string> >("SrcTermComds","Types of the source term commands.");
  options.addConfigOption< std::vector<std::string> >("SrcTermNames","Names of the source term commands.");
  options.addConfigOption< std::vector<std::string> >("InitComds","Types of the initializing commands.");
  options.addConfigOption< std::vector<std::string> >("InitNames","Names of the initializing commands.");
  options.addConfigOption< std::vector<std::string> >("BcNames","Names of the boundary condition commands.");
  options.addConfigOption< std::string >("SpaceRHSJacobCom","Command for the computation of the space discretization contibution to RHS and Jacobian.");
  options.addConfigOption< std::string >("LimiterCom","Command to limit the solution.");
}

//////////////////////////////////////////////////////////////////////////////

FluxReconstructionSolver::FluxReconstructionSolver(const std::string& name) :
  SpaceMethod(name),
  m_setup(),
  m_unsetup(),
  m_solve(),
  m_limiter(),
  m_inits(),
  m_srcTerms(),
  m_convVolTerm(),
  m_convFaceTerm(),
  m_bcs(),
  m_bcsComs()
{
  addConfigOptionsTo(this);
  m_data.reset(new FluxReconstructionSolverData(this));

  cf_assert(m_data.isNotNull());

  // set default value of builder for FluxReconstructionSolver
  // to be FluxReconstructionBuilder
  m_builder = "StdBuilder";
  // set default global jacobian sparsity
  m_sparsity = "None";
  
  m_spaceRHSJacobStr = "RHS";
  setParameter("SpaceRHSJacobCom", &m_spaceRHSJacobStr);
  
  m_limiterStr = "Null";
  setParameter("LimiterCom", &m_limiterStr);

  m_setupStr   = "StdSetup";
  setParameter( "SetupCom",   &m_setupStr );

  m_unsetupStr = "StdUnSetup";
  setParameter( "UnSetupCom", &m_unsetupStr );

  m_solveStr   = "StdSolve";
  setParameter( "SolveCom",   &m_solveStr );
  
  // options for source term commands
  m_srcTermTypeStr = std::vector<std::string>();
  setParameter("SrcTermComds",&m_srcTermTypeStr);

  m_srcTermNameStr = std::vector<std::string>();
  setParameter("SrcTermNames",&m_srcTermNameStr);
  
  // options for initialize commands
  m_initTypeStr = std::vector<std::string>();
  setParameter("InitComds",&m_initTypeStr);

  m_initNameStr = std::vector<std::string>();
  setParameter("InitNames",&m_initNameStr);
  
  // options for bc commands
  m_bcNameStr = std::vector<std::string>();
  setParameter("BcNames",&m_bcNameStr);
}

//////////////////////////////////////////////////////////////////////////////

FluxReconstructionSolver::~FluxReconstructionSolver()
{
}

//////////////////////////////////////////////////////////////////////////////

void FluxReconstructionSolver::setCollaborator( MultiMethodHandle<LinearSystemSolver> lss )
{
  m_data->setLinearSystemSolver(lss);
}

//////////////////////////////////////////////////////////////////////////////

void FluxReconstructionSolver::setCollaborator(MultiMethodHandle<ConvergenceMethod> convMtd)
{
  // convergence method collaborator is made available to the commands through
  // the method data
  m_data->setConvergenceMethod(convMtd);
}

//////////////////////////////////////////////////////////////////////////////

void FluxReconstructionSolver::configure ( Config::ConfigArgs& args )
{
  SpaceMethod::configure(args);

  configureNested ( m_data.getPtr(), args );

  // add here configures to the FluxReconstructionSolver
  configureCommand< FluxReconstructionSolverData,FluxReconstructionSolverCom::PROVIDER >(
    args, m_setup,m_setupStr,m_data );
  configureCommand< FluxReconstructionSolverData,FluxReconstructionSolverCom::PROVIDER >(
    args, m_unsetup,m_unsetupStr,m_data );
  configureCommand< FluxReconstructionSolverData,FluxReconstructionSolverCom::PROVIDER >(
    args, m_solve,m_solveStr,m_data );
  configureCommand< FluxReconstructionSolverData,FluxReconstructionSolverCom::PROVIDER >( args, m_limiter,m_limiterStr,m_data );

  cf_assert(m_setup.isNotNull());
  cf_assert(m_unsetup.isNotNull());
  cf_assert(m_solve.isNotNull());
  cf_assert(m_limiter.isNotNull());
  
  configureSourceTermCommands(args);
  configureInitCommands(args);
  
  CFLog(INFO,"FR: Creating convective volume term command...\n");
    CFLog(INFO,"ConvVolTerm" << m_spaceRHSJacobStr << "\n");
    try
    {
      configureCommand< FluxReconstructionSolverData,FluxReconstructionSolverCom::PROVIDER >( args,
        m_convVolTerm,"ConvVolTerm"+m_spaceRHSJacobStr,m_data );
    }
    catch (Common::NoSuchValueException& e)
    {
      CFLog(INFO, e.what() << "\n");
      CFLog(INFO, "Choosing ConvVolTermRHS instead ...\n");

      configureCommand< FluxReconstructionSolverData,FluxReconstructionSolverCom::PROVIDER >( args, m_convVolTerm,"ConvVolTermRHS",m_data );
    }

    CFLog(INFO,"FR: Creating convective face term command...\n");
    CFLog(INFO,"ConvFaceTerm" << m_spaceRHSJacobStr << "\n");
    try
    {
      configureCommand< FluxReconstructionSolverData,FluxReconstructionSolverCom::PROVIDER >( args, m_convFaceTerm,"ConvFaceTerm"+m_spaceRHSJacobStr,m_data );
    }
    catch (Common::NoSuchValueException& e)
    {
      CFLog(INFO, e.what() << "\n");
      CFLog(INFO, "Choosing ConvFaceTermRHS instead ...\n");

      configureCommand< FluxReconstructionSolverData,FluxReconstructionSolverCom::PROVIDER >( args, m_convFaceTerm,"ConvFaceTermRHS",m_data );
    }
    cf_assert(m_convVolTerm.isNotNull());
    cf_assert(m_convFaceTerm.isNotNull());
    
    configureBcCommands(args);
}

//////////////////////////////////////////////////////////////////////////////

void FluxReconstructionSolver::configureSourceTermCommands ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;
  cf_assert(m_srcTermTypeStr.size() == m_srcTermNameStr.size());

  m_srcTerms.resize(m_srcTermTypeStr.size());

  for(CFuint i = 0; i < m_srcTerms.size(); ++i)
  {

    CFLog(INFO, "SOURCE TERM type = " << m_srcTermTypeStr[i] << "\n");
    CFLog(INFO, "SOURCE TERM name = " << m_srcTermNameStr[i] << "\n");

    configureCommand<FluxReconstructionSolverCom,
                     FluxReconstructionSolverData,
                     FluxReconstructionSolverComProvider>
        (args, m_srcTerms[i], m_srcTermTypeStr[i],m_srcTermNameStr[i], m_data);

    cf_assert(m_srcTerms[i].isNotNull());
  }
}

//////////////////////////////////////////////////////////////////////////////

void FluxReconstructionSolver::configureInitCommands ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;
  cf_assert(m_initTypeStr.size() == m_initNameStr.size());

  m_inits.resize(m_initTypeStr.size());

  for(CFuint i = 0; i < m_inits.size(); ++i)
  {

    CFLog(INFO, "INIT type = " << m_initTypeStr[i] << "\n");
    CFLog(INFO, "INIT name = " << m_initNameStr[i] << "\n");

    configureCommand<FluxReconstructionSolverCom,
      FluxReconstructionSolverData,
      FluxReconstructionSolverComProvider>
      (args, m_inits[i], m_initTypeStr[i],m_initNameStr[i], m_data);

    cf_assert(m_inits[i].isNotNull());
  }
}

//////////////////////////////////////////////////////////////////////////////

void FluxReconstructionSolver::configureBcCommands ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  // get bcStateComputers
  SafePtr< std::vector< SafePtr< BCStateComputer > > > bcStateComputers = m_data->getBCStateComputers();
  cf_assert(m_bcNameStr.size() == bcStateComputers->size());

  // resize vector with bc commands
  m_bcsComs.resize(m_bcNameStr.size());

  // get bc TRS names variable from the method data
  SafePtr< std::vector< std::vector< std::string > > > bcTRSNames = m_data->getBCTRSNameStr();
  bcTRSNames->resize(m_bcNameStr.size());

  // configure commands
    m_bcs.resize(m_bcNameStr.size());

    for(CFuint iBc = 0; iBc < m_bcsComs.size(); ++iBc)
    {
      CFLog(INFO,"FluxReconstruction: Creating convective boundary face term command for boundary condition: "
                  << m_bcNameStr[iBc] << "\n");
      CFLog(INFO,"ConvBndFaceTerm" << m_spaceRHSJacobStr << "\n");
      try
      {
        configureCommand<FluxReconstructionSolverCom,
          FluxReconstructionSolverData,
          FluxReconstructionSolverComProvider>
          (args, m_bcsComs[iBc], "ConvBndFaceTerm"+m_spaceRHSJacobStr,m_bcNameStr[iBc], m_data);
      }
      catch (Common::NoSuchValueException& e)
      {
        CFLog(INFO, e.what() << "\n");
        CFLog(INFO, "Choosing ConvBndFaceTermRHS instead ...\n");

        configureCommand<FluxReconstructionSolverCom,
          FluxReconstructionSolverData,
          FluxReconstructionSolverComProvider>
          (args, m_bcsComs[iBc], "ConvBndFaceTermRHS",m_bcNameStr[iBc], m_data);
      }

      cf_assert(m_bcsComs[iBc].isNotNull());

      // dynamic_cast to ConvBndFaceTermRHSFluxReconstruction
      SafePtr< FluxReconstructionSolverCom > bcComm = m_bcsComs[iBc].getPtr();
      m_bcs[iBc] = bcComm.d_castTo< ConvBndFaceTermRHSFluxReconstruction >();
      cf_assert(m_bcs[iBc].isNotNull());

      // set bcStateComputer corresponding to this bc command
      m_bcs[iBc]->setBcStateComputer((*bcStateComputers)[iBc]);

      // set TRS names corresponding to this command in bcTRSNames and in bcStateComputers
      const std::vector<std::string> currBCTRSNames = m_bcs[iBc]->getTrsNames();
      const CFuint nbrBCTRSs = currBCTRSNames.size();
      (*bcTRSNames)[iBc].resize(nbrBCTRSs);
      for (CFuint iTRS = 0; iTRS < nbrBCTRSs; ++iTRS)
      {
        (*bcTRSNames)[iBc][iTRS] = currBCTRSNames[iTRS];
        (*bcStateComputers)[iBc]->addTRSName(currBCTRSNames[iTRS]);
      }
    }
}

//////////////////////////////////////////////////////////////////////////////

void FluxReconstructionSolver::extrapolateStatesToNodesImpl()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void FluxReconstructionSolver::initializeSolutionImpl(bool isRestart)
{
  CFAUTOTRACE;
  
  cf_assert(isConfigured());
  cf_assert(isSetup());

  if (!isRestart)
  {
    for(CFuint i = 0; i < m_inits.size(); ++i)
    {
      cf_assert(m_inits[i].isNotNull());
      m_inits[i]->execute();
    }
  }

  // apply a limiter to the solution
  cf_assert(m_limiter.isNotNull());
  m_limiter->execute();
}

//////////////////////////////////////////////////////////////////////////////

void FluxReconstructionSolver::computeSpaceResidualImpl(CFreal factor)
{
  CFAUTOTRACE;
  cf_assert(m_solve.isNotNull());
  m_solve->execute();
}

//////////////////////////////////////////////////////////////////////////////

void FluxReconstructionSolver::computeTimeResidualImpl(CFreal factor)
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void FluxReconstructionSolver::applyBCImpl()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void FluxReconstructionSolver::prepareComputationImpl()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

Common::Signal::return_t FluxReconstructionSolver::beforeMeshUpdateActionImpl(Common::Signal::arg_t eBefore)
{
  CFAUTOTRACE;
  return Common::Signal::return_t ();
}

//////////////////////////////////////////////////////////////////////////////

Common::Signal::return_t FluxReconstructionSolver::afterMeshUpdateActionImpl(Common::Signal::arg_t eAfter)
{
  CFAUTOTRACE;
  return Common::Signal::return_t ();
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<NumericalStrategy> > FluxReconstructionSolver::getStrategyList() const
{
  std::vector<Common::SafePtr<NumericalStrategy> > result;

  // add strategies here
  result.push_back(m_data->getStatesReconstructor()  .d_castTo<NumericalStrategy>());
  result.push_back(m_data->getInterfaceFlux()  .d_castTo<NumericalStrategy>());
  
  return result;
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

