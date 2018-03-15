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
#include "FluxReconstructionMethod/BaseCorrectionFunction.hh"
#include "FluxReconstructionMethod/BCStateComputer.hh"
#include "FluxReconstructionMethod/ConvBndCorrectionsRHSFluxReconstruction.hh"
#include "FluxReconstructionMethod/DiffBndCorrectionsRHSFluxReconstruction.hh"
#include "FluxReconstructionMethod/RiemannFlux.hh"

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
  options.addConfigOption< std::string >("ExtrapolateCom","Command to extrapolate the states values to the position at the nodes.");
  options.addConfigOption< std::string >("PrepareCom","Command to prepare before the computation of the residuals.");
  options.addConfigOption< std::string >("ConvSolveCom","Command to solve the convection problem with FluxReconstruction solver.");
  options.addConfigOption< std::string >("DiffSolveCom","Command to solve the diffusion problem with FluxReconstruction solver.");
  options.addConfigOption< std::string >("UnSetupCom","Command to deallocate FluxReconstruction solver data.");
  options.addConfigOption< std::string >("SetupCom","Command to initialize FluxReconstruction solver data.");
  options.addConfigOption< std::vector<std::string> >("SrcTermComds","Types of the source term commands.");
  options.addConfigOption< std::vector<std::string> >("SrcTermNames","Names of the source term commands.");
  options.addConfigOption< std::vector<std::string> >("InitComds","Types of the initializing commands.");
  options.addConfigOption< std::vector<std::string> >("InitNames","Names of the initializing commands.");
  options.addConfigOption< std::vector<std::string> >("BcNames","Names of the boundary condition commands.");
  options.addConfigOption< std::vector<std::string> >("BcNamesDiff","Names of the diffusive boundary condition commands.");
  options.addConfigOption< std::string >("SpaceRHSJacobCom","Command for the computation of the space discretization contribution to RHS and Jacobian.");
  options.addConfigOption< std::string >("TimeRHSJacobCom","Command for the computation of the time discretization contibution to RHS and Jacobian.");
  options.addConfigOption< std::string >("LimiterCom","Command to limit the solution.");
  options.addConfigOption< std::string >("PhysicalityCom","Command to enforce physical soundness of the solution.");
  options.addConfigOption< std::string >("ComputeErrorCom","Command to compute the error of the solution.");
  options.addConfigOption< std::string >("FinalizeRHSCom","Finilaze computation of the RHS.");
  options.addConfigOption< std::string >("ArtificialViscosityCom","Command to add artificial viscosity.");
}

//////////////////////////////////////////////////////////////////////////////

FluxReconstructionSolver::FluxReconstructionSolver(const std::string& name) :
  SpaceMethod(name),
  m_setup(),
  m_unsetup(),
  m_extrapolate(),
  m_prepare(),
  m_convSolve(),
  m_diffSolve(),
  m_limiter(),
  m_srcTerms(),
  m_inits(),
  m_bcs(),
  m_timeRHSJacob(),
  m_bcsComs(),
  m_bcsDiff(),
  m_bcsDiffComs(),
  m_computeError(),
  m_finalizeRHS(),
  m_artificialVisc()
{
  addConfigOptionsTo(this);
  m_data.reset(new FluxReconstructionSolverData(this));

  cf_assert(m_data.isNotNull());

  // set default value of builder for FluxReconstructionSolver to be FluxReconstructionBuilder
  m_builder = "StdBuilder";
  
  // set default global jacobian sparsity
  m_sparsity = "None";
  
  m_spaceRHSJacobStr = "RHS";
  setParameter("SpaceRHSJacobCom", &m_spaceRHSJacobStr);
  
  m_timeRHSJacobStr = "Null";
  setParameter("TimeRHSJacobCom", &m_timeRHSJacobStr);
  
  m_limiterStr = "Null";
  setParameter("LimiterCom", &m_limiterStr);
  
  m_artificialViscStr = "Null";
  setParameter("ArtificialViscosityCom", &m_artificialViscStr);
  
  m_physicalityStr = "Null";
  setParameter("PhysicalityCom", &m_physicalityStr);

  m_setupStr   = "StdSetup";
  setParameter( "SetupCom",   &m_setupStr );

  m_unsetupStr = "StdUnSetup";
  setParameter( "UnSetupCom", &m_unsetupStr );
  
  m_extrapolateStr = "StdExtrapolate";
  setParameter("ExtrapolateCom", &m_extrapolateStr);
  
  m_prepareStr = "StdPrepare";
  setParameter("PrepareCom", &m_prepareStr);

  m_convSolveStr   = "ConvRHS";
  setParameter( "ConvSolveCom",   &m_convSolveStr );
  
  m_diffSolveStr   = "DiffRHS";
  setParameter( "DiffSolveCom",   &m_diffSolveStr );
  
  m_computeErrorStr = "Null";
  setParameter("ComputeErrorCom", &m_computeErrorStr);
  
  m_finalizeRHSStr = "StdFinalize";
  setParameter("FinalizeRHSCom", &m_finalizeRHSStr);
  
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
  
  // options for bc commands
  m_bcNameDiffStr = std::vector<std::string>();
  setParameter("BcNamesDiff",&m_bcNameDiffStr);
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
    args, m_extrapolate,m_extrapolateStr,m_data );
  configureCommand< FluxReconstructionSolverData,FluxReconstructionSolverCom::PROVIDER >( 
    args, m_prepare,m_prepareStr,m_data );
  configureCommand< FluxReconstructionSolverData,FluxReconstructionSolverCom::PROVIDER >(
    args, m_convSolve,m_convSolveStr,m_data );
  configureCommand< FluxReconstructionSolverData,FluxReconstructionSolverCom::PROVIDER >(
    args, m_diffSolve,m_diffSolveStr,m_data );
  configureCommand< FluxReconstructionSolverData,FluxReconstructionSolverCom::PROVIDER >( 
    args, m_limiter,m_limiterStr,m_data );
  configureCommand< FluxReconstructionSolverData,FluxReconstructionSolverCom::PROVIDER >( 
    args, m_artificialVisc,m_artificialViscStr,m_data );
  configureCommand< FluxReconstructionSolverData,FluxReconstructionSolverCom::PROVIDER >( 
    args, m_physicality,m_physicalityStr,m_data );
  configureCommand< FluxReconstructionSolverData,FluxReconstructionSolverCom::PROVIDER >( 
    args, m_timeRHSJacob,m_timeRHSJacobStr,m_data );
  configureCommand< FluxReconstructionSolverData,FluxReconstructionSolverCom::PROVIDER >( 
    args, m_computeError,m_computeErrorStr,m_data );
  configureCommand< FluxReconstructionSolverData,FluxReconstructionSolverCom::PROVIDER >( 
    args, m_finalizeRHS,m_finalizeRHSStr,m_data );

  cf_assert(m_setup.isNotNull());
  cf_assert(m_unsetup.isNotNull());
  cf_assert(m_extrapolate.isNotNull());
  cf_assert(m_prepare.isNotNull());
  cf_assert(m_convSolve.isNotNull());
  cf_assert(m_diffSolve.isNotNull());
  cf_assert(m_limiter.isNotNull());
  cf_assert(m_artificialVisc.isNotNull());
  cf_assert(m_physicality.isNotNull());
  cf_assert(m_timeRHSJacob.isNotNull());
  cf_assert(m_computeError.isNotNull());
  cf_assert(m_finalizeRHS.isNotNull());
  
  configureSourceTermCommands(args);
  configureInitCommands(args);
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
  if (m_bcNameDiffStr.size() != m_bcNameStr.size())
  {
    CFLog(NOTICE,"Number of Diffusive BCs doesn't match convective BCs: should only happen when there is no diffusive term!\n");
    m_bcNameDiffStr.resize(m_bcNameStr.size());
    for (CFuint iBC = 0; iBC < m_bcNameStr.size(); ++iBC)
    {
      m_bcNameDiffStr[iBC] = m_bcNameStr[iBC];
    }
  }
  CFLog(VERBOSE, "Nbr of BCs: " << bcStateComputers->size() << "\n");

  // resize vector with bc commands
  m_bcsComs.resize(m_bcNameStr.size());

  // get bc TRS names variable from the method data
  SafePtr< std::vector< std::vector< std::string > > > bcTRSNames = m_data->getBCTRSNameStr();
  bcTRSNames->resize(m_bcNameStr.size());

  // configure commands
    m_bcs.resize(m_bcNameStr.size());
    m_bcsDiffComs.resize(m_bcNameStr.size());
    m_bcsDiff.resize(m_bcNameStr.size());
    std::vector<std::string> TRSsConv;

    for(CFuint iBc = 0; iBc < m_bcsComs.size(); ++iBc)
    {
      CFLog(INFO,"FluxReconstruction: Creating convective boundary correction command for boundary condition: "
                  << m_bcNameStr[iBc] << "\n");
      CFLog(INFO,"ConvBndCorrections" << m_spaceRHSJacobStr << "\n");

      try
      {
        configureCommand<FluxReconstructionSolverCom,
          FluxReconstructionSolverData,
          FluxReconstructionSolverComProvider>
          (args, m_bcsComs[iBc], "ConvBndCorrections"+m_spaceRHSJacobStr,m_bcNameStr[iBc], m_data);
      }
      catch (Common::NoSuchValueException& e)
      {
        CFLog(INFO, e.what() << "\n");
        CFLog(INFO, "Choosing ConvBndCorrectionsRHS instead ...\n");

        configureCommand<FluxReconstructionSolverCom,
          FluxReconstructionSolverData,
          FluxReconstructionSolverComProvider>
          (args, m_bcsComs[iBc], "ConvBndCorrectionsRHS",m_bcNameStr[iBc], m_data);
      }

      cf_assert(m_bcsComs[iBc].isNotNull());

      // dynamic_cast to ConvBndCorrectionsRHSFluxReconstruction
      SafePtr< FluxReconstructionSolverCom > bcComm = m_bcsComs[iBc].getPtr();
      m_bcs[iBc] = bcComm.d_castTo< ConvBndCorrectionsRHSFluxReconstruction >();
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
      
      CFLog(INFO,"FluxReconstruction: Creating diffusive boundary correction command for boundary condition: "
                  << m_bcNameDiffStr[iBc] << "\n");
      CFLog(INFO,"DiffBndCorrections" << m_spaceRHSJacobStr << "\n");
      try
      {
        configureCommand<FluxReconstructionSolverCom,
          FluxReconstructionSolverData,
          FluxReconstructionSolverComProvider>
          (args, m_bcsDiffComs[iBc], "DiffBndCorrections"+m_spaceRHSJacobStr,m_bcNameDiffStr[iBc], m_data);
      }
      catch (Common::NoSuchValueException& e)
      {
        CFLog(INFO, e.what() << "\n");
        CFLog(INFO, "Choosing DiffBndCorrectionsRHS instead ...\n");

        configureCommand<FluxReconstructionSolverCom,
          FluxReconstructionSolverData,
          FluxReconstructionSolverComProvider>
          (args, m_bcsDiffComs[iBc], "DiffBndCorrectionsRHS",m_bcNameDiffStr[iBc], m_data);
      }

      cf_assert(m_bcsDiffComs[iBc].isNotNull());

      // dynamic_cast to DiffBndCorrectionsRHSFluxReconstruction
      SafePtr< FluxReconstructionSolverCom > bcDiffComm = m_bcsDiffComs[iBc].getPtr();
      m_bcsDiff[iBc] = bcDiffComm.d_castTo< DiffBndCorrectionsRHSFluxReconstruction >();
      cf_assert(m_bcsDiff[iBc].isNotNull());

      // set bcStateComputer corresponding to this bc command
      m_bcsDiff[iBc]->setBcStateComputer((*bcStateComputers)[iBc]);
      
      TRSsConv = m_bcsComs[iBc]->getTrsNames();
      for (CFuint i = 0; i < TRSsConv.size(); ++i)
      {
	CFLog(VERBOSE, "Conv TRS " << i << ": " << TRSsConv[i] << "\n");
      }
      const std::vector<std::string> TRSsDiff = m_bcsDiffComs[iBc]->getTrsNames();
      for (CFuint i = 0; i < TRSsDiff.size(); ++i)
      {
	CFLog(VERBOSE, "Diff TRS " << i << ": " << TRSsDiff[i] << "\n");
      }
    }
}

//////////////////////////////////////////////////////////////////////////////

void FluxReconstructionSolver::extrapolateStatesToNodesImpl()
{
  CFAUTOTRACE;
  
  cf_assert(m_extrapolate.isNotNull());
  m_extrapolate->execute();
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
  
  cf_assert(isConfigured());
  cf_assert(isSetup());
  
  // set the residual factor in the MethodData
  m_data->setResFactor(factor);
  
  // prepare the computation (reset rhs, updatecoeff, grads)
  cf_assert(m_prepare.isNotNull());
  m_prepare->execute();
  
  // enforce physicality of the solution
  cf_assert(m_physicality.isNotNull());
  m_physicality->execute();
  
  // apply the boundary conditions (this function is in SpaceMethod and is not called anywhere else)
  applyBC();
  
  cf_assert(m_convSolve.isNotNull());
  m_convSolve->execute();

  // if there is a diffusive term, compute the diffusive contributions to the residual
  if (m_data->hasDiffTerm())
  {
    // add the diffusive boundary fluxes
    applyBCDiffImpl();

    cf_assert(m_diffSolve.isNotNull());
    m_diffSolve->execute();
  }
  
  if (m_artificialViscStr != "Null")
  {
    cf_assert(m_artificialVisc.isNotNull());
    m_artificialVisc->execute();
  }
  
  // add source terms
  addSourceTermsImpl();
  
  cf_assert(m_computeError.isNotNull());
  m_computeError->execute();
  
  if (m_data->isResidualTransformationNeeded())
  {
    m_finalizeRHS->execute();
  }
  
  cf_assert(m_limiter.isNotNull());
  //m_limiter->execute();

//   // divide by volume/Jacobian determinant
//   m_divideRHSByCellVol->execute();
}

//////////////////////////////////////////////////////////////////////////////

void FluxReconstructionSolver::computeTimeResidualImpl(CFreal factor)
{
  CFAUTOTRACE;
  
  cf_assert(isConfigured());
  cf_assert(isSetup());

  // set the residual factor in the MethodData
  m_data->setResFactor(factor);

  // compute the time contribution to the jacobian
  m_timeRHSJacob->execute();
}

//////////////////////////////////////////////////////////////////////////////

void FluxReconstructionSolver::applyBCImpl()
{
  CFAUTOTRACE;
  
  const CFuint nbrBcs = m_bcsComs.size();
  for(CFuint iBc = 0; iBc < nbrBcs; ++iBc)
  {
    cf_assert(m_bcsComs[iBc].isNotNull());
    m_bcsComs[iBc]->execute();
  }
}

//////////////////////////////////////////////////////////////////////////////

void FluxReconstructionSolver::applyBCDiffImpl()
{
  CFAUTOTRACE;

  const CFuint nbrBcs = m_bcsDiffComs.size();
  for(CFuint iBc = 0; iBc < nbrBcs; ++iBc)
  {
    CFLog(VERBOSE,"Executing BC " << iBc << "\n");
    cf_assert(m_bcsDiffComs[iBc].isNotNull());
    m_bcsDiffComs[iBc]->execute();
  }
}

//////////////////////////////////////////////////////////////////////////////

void FluxReconstructionSolver::addSourceTermsImpl()
{
  CFAUTOTRACE;

  const CFuint nbrSrcTerms = m_srcTerms.size();
  for(CFuint iSrc = 0; iSrc < nbrSrcTerms; ++iSrc)
  {
    cf_assert(m_srcTerms[iSrc].isNotNull());
    m_srcTerms[iSrc]->execute();
  }
}

//////////////////////////////////////////////////////////////////////////////

void FluxReconstructionSolver::prepareComputationImpl()
{
  CFAUTOTRACE;

  cf_assert(m_prepare.isNotNull());
  m_prepare->execute();
}

//////////////////////////////////////////////////////////////////////////////

void FluxReconstructionSolver::postProcessSolutionImpl()
{
  CFAUTOTRACE;

  cf_assert(m_limiter.isNotNull());
  m_limiter->execute();
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
  result.push_back(m_data->getSolPntDistribution()   .d_castTo<NumericalStrategy>());
  result.push_back(m_data->getFluxPntDistribution()  .d_castTo<NumericalStrategy>());
  result.push_back(m_data->getCorrectionFunction()   .d_castTo<NumericalStrategy>());
  
  // add BCStateComputers
  SafePtr< std::vector< SafePtr< BCStateComputer > > > bcStateComputers = m_data->getBCStateComputers();
  for (CFuint iBC = 0; iBC < bcStateComputers->size(); ++iBC)
  {
    result.push_back((*bcStateComputers)[iBC]        .d_castTo<NumericalStrategy>());
  }
  
  // these have to be added last!!! (need data from previous strategies)
  result.push_back(m_data->getRiemannFlux()        .d_castTo<NumericalStrategy>());
  
  return result;
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

