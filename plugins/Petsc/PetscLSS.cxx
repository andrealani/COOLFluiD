// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Petsc/PetscHeaders.hh" // must come before any header

#include "Common/PE.hh"
#include "Common/Stopwatch.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/BlockAccumulator.hh"

#include "Petsc/Petsc.hh"
#include "Petsc/PetscLSS.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace Petsc {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<PetscLSS,
               LinearSystemSolver,
               PetscModule,
               1>
petscLSSMethodProvider("PETSC");

//////////////////////////////////////////////////////////////////////////////

void PetscLSS::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::string >( "SetupCom",
     "Setup Command to run. This command seldomly needs overriding." );
   options.addConfigOption< std::string >( "UnSetupCom",
     "UnSetup Command to run. This command seldomly needs overriding." );
   options.addConfigOption< std::string >( "SysSolver",
     "Command that solves the linear system." );
}

//////////////////////////////////////////////////////////////////////////////

PetscLSS::PetscLSS(const std::string& name)
  : LinearSystemSolver(name)
{
  addConfigOptionsTo(this);

  m_data.reset(new PetscLSSData(getMaskArray(),
               getNbSysEquations(),
               this));

  cf_assert(m_data.getPtr() != CFNULL);

  m_setupStr = "Null";
  setParameter("SetupCom",&m_setupStr);

  m_unSetupStr = "Null";
  setParameter("UnSetupCom",&m_unSetupStr);

  m_solveSysStr = "Null";
  setParameter("SysSolver",&m_solveSysStr);
}

//////////////////////////////////////////////////////////////////////////////

PetscLSS::~PetscLSS()
{
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<MethodData> PetscLSS::getMethodData () const
{
  return m_data.getPtr();
}

//////////////////////////////////////////////////////////////////////////////

void PetscLSS::configure ( Config::ConfigArgs& args )
{
  LinearSystemSolver::configure(args);
  
  m_data->setFactoryRegistry(getFactoryRegistry());
  configureNested ( m_data.getPtr(), args );

  const bool isParallel = Common::PE::GetPE().IsParallel ();

  if (m_setupStr == "Null") {
    m_setupStr = (isParallel) ? "NewParSetup" : "NewSeqSetup";
  }

  if (m_unSetupStr == "Null") {
    m_unSetupStr = (isParallel) ? "StdParUnSetup" : "StdSeqUnSetup";
  }

  if (m_solveSysStr == "Null") {
    m_solveSysStr = (isParallel) ? "StdParSolveSys" : "StdSeqSolveSys";
  }

  configureCommand<PetscLSSData,PetscLSSComProvider>( args, m_setup,m_setupStr,m_data);

  configureCommand<PetscLSSData,PetscLSSComProvider>( args, m_unSetup,m_unSetupStr,m_data);

  configureCommand<PetscLSSData,PetscLSSComProvider>( args, m_solveSys,m_solveSysStr,m_data);
}

//////////////////////////////////////////////////////////////////////////////

void PetscLSS::solveSysImpl()
{
  /// @TODO each processor should print in separate files
  if(m_data->isSaveSystemToFile()) {
    CFLog(NOTICE, "Printing system to files... " << CFendl);
    std::string prefix = "system-";
    std::string suffix =  "-" + getName() + ".dat";
    printToFile(prefix,suffix);
    CFLog(NOTICE, "Done !!!\n");
  }

  Common::Stopwatch<Common::WallTime> stopTimer;
  stopTimer.start();

  CFLog(DEBUG_MAX, "Solving LSS: " << getName() << CFendl);
  m_solveSys->execute();

  stopTimer.stop ();

  CFLog(VERBOSE, "PetscLSS::solveSys() WallTime: " << stopTimer << "s\n");
}

//////////////////////////////////////////////////////////////////////////////

BlockAccumulator* PetscLSS::
createBlockAccumulator(const CFuint nbRows,
                       const CFuint nbCols,
                       const CFuint subBlockSize,
		       CFreal* ptr) const
{
  if (!m_data->useBlockPreconditionerMatrix()) { 
    return new BlockAccumulator
      (nbRows,nbCols,subBlockSize, m_lssData->getLocalToGlobalMapping(), ptr);
  }
  return new BlockAccumulator
    (nbRows,nbCols,subBlockSize, m_lssData->getLocalToLocallyUpdatableMapping(), ptr);
}

//////////////////////////////////////////////////////////////////////////////

void PetscLSS::printToFile(const std::string prefix, const std::string suffix)
{
  cf_assert(isSetup());
  cf_assert(isConfigured());

  PetscMatrix& mat = m_data->getMatrix();
  PetscVector& rhs = m_data->getRhsVector();
  PetscVector& sol = m_data->getSolVector();

  std::string matStr = prefix + "mat" + suffix;
  std::string rhsStr = prefix + "rhs" + suffix;
  std::string solStr = prefix + "sol" + suffix;

  mat.printToFile(matStr.c_str());
  rhs.printToFile(rhsStr.c_str());
  sol.printToFile(solStr.c_str());
}

//////////////////////////////////////////////////////////////////////////////

void PetscLSS::setMethodImpl()
{

  LinearSystemSolver::setMethodImpl();

//  setupCommandsAndStrategies();

  m_setup->setup();
  m_setup->execute();

  m_solveSys->setup();
  m_unSetup->setup();
  
  m_data->getShellPreconditioner()->setup();
}

//////////////////////////////////////////////////////////////////////////////

void PetscLSS::unsetMethodImpl()
{
  m_unSetup->execute();
  unsetupCommandsAndStrategies();

  LinearSystemSolver::unsetMethodImpl();
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<NumericalStrategy> > PetscLSS::getStrategyList() const
{
  vector<Common::SafePtr<NumericalStrategy> > result;

  // add strategies here
  result.push_back(m_data->getShellPreconditioner().d_castTo<NumericalStrategy>());
  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace Petsc

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

