// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Environment/ObjectProvider.hh"
#include "EmptySpaceMethod/Empty.hh"
#include "EmptySpaceMethod/EmptySolver.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace EmptySpaceMethod {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider< EmptySolver,SpaceMethod,EmptyModule,1 >
  emptySolverProvider("EmptySolver");

//////////////////////////////////////////////////////////////////////////////

void EmptySolver::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::string >("SolveCom","Command to solve the problem with Empty solver.");
  options.addConfigOption< std::string >("UnSetupCom","Command to deallocate Empty solver data.");
  options.addConfigOption< std::string >("SetupCom","Command to initialize Empty solver data.");
}

//////////////////////////////////////////////////////////////////////////////

EmptySolver::EmptySolver(const std::string& name) :
  SpaceMethod(name),
  m_setup(),
  m_unsetup(),
  m_solve()
{
  addConfigOptionsTo(this);
  m_data.reset(new EmptySolverData(this));

  cf_assert(m_data.isNotNull());

  // set default value of builder for EmptySolver
  // to be EmptyMeshDataBuilder
  m_builder = "Empty";
  // set default global jacobian sparsity
  m_sparsity = "None";

  m_setupStr   = "StdSetup";
  setParameter( "SetupCom",   &m_setupStr );

  m_unsetupStr = "StdUnSetup";
  setParameter( "UnSetupCom", &m_unsetupStr );

  m_solveStr   = "StdSolve";
  setParameter( "SolveCom",   &m_solveStr );
}

//////////////////////////////////////////////////////////////////////////////

EmptySolver::~EmptySolver()
{
}

//////////////////////////////////////////////////////////////////////////////

void EmptySolver::setCollaborator( MultiMethodHandle<LinearSystemSolver> lss )
{
  m_data->setLinearSystemSolver(lss);
}

//////////////////////////////////////////////////////////////////////////////

void EmptySolver::setCollaborator(MultiMethodHandle<ConvergenceMethod> convMtd)
{
  // convergence method collaborator is made available to the commands through
  // the method data
  m_data->setConvergenceMethod(convMtd);
}

//////////////////////////////////////////////////////////////////////////////

void EmptySolver::configure ( Config::ConfigArgs& args )
{
  SpaceMethod::configure(args);

  configureNested ( m_data.getPtr(), args );

  // add here configures to the EmptySolver
  configureCommand< EmptySolverData,EmptySolverCom::PROVIDER >(
    args, m_setup,m_setupStr,m_data );
  configureCommand< EmptySolverData,EmptySolverCom::PROVIDER >(
    args, m_unsetup,m_unsetupStr,m_data );
  configureCommand< EmptySolverData,EmptySolverCom::PROVIDER >(
    args, m_solve,m_solveStr,m_data );

  cf_assert(m_setup.isNotNull());
  cf_assert(m_unsetup.isNotNull());
  cf_assert(m_solve.isNotNull());
}

//////////////////////////////////////////////////////////////////////////////

void EmptySolver::setMethodImpl()
{
  CFAUTOTRACE;

  SpaceMethod::setMethodImpl();

  setupCommandsAndStrategies();
  cf_assert(m_setup.isNotNull());
  m_setup->execute();
}

//////////////////////////////////////////////////////////////////////////////

void EmptySolver::unsetMethodImpl()
{
  CFAUTOTRACE;

  cf_assert(m_unsetup.isNotNull());
  m_unsetup->execute();
  unsetupCommandsAndStrategies();

  SpaceMethod::setMethodImpl();
}

//////////////////////////////////////////////////////////////////////////////

void EmptySolver::extrapolateStatesToNodesImpl()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void EmptySolver::initializeSolutionImpl(bool isRestart)
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void EmptySolver::computeSpaceResidualImpl(CFreal factor)
{
  CFAUTOTRACE;
  cf_assert(m_solve.isNotNull());
  m_solve->execute();
}

//////////////////////////////////////////////////////////////////////////////

void EmptySolver::computeTimeResidualImpl(CFreal factor)
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void EmptySolver::applyBCImpl()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void EmptySolver::prepareComputationImpl()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

Common::Signal::return_t EmptySolver::beforeMeshUpdateActionImpl(Common::Signal::arg_t eBefore)
{
  CFAUTOTRACE;
  return Common::Signal::return_t ();
}

//////////////////////////////////////////////////////////////////////////////

Common::Signal::return_t EmptySolver::afterMeshUpdateActionImpl(Common::Signal::arg_t eAfter)
{
  CFAUTOTRACE;
  return Common::Signal::return_t ();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace EmptySpaceMethod
}  // namespace COOLFluiD

