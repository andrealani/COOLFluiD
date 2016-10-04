// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Common/PE.hh"
#include "Common/Stopwatch.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/BlockAccumulator.hh"

#include "Paralution/Paralution.hh"
#include "Paralution/ParalutionLSS.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace Paralution {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<ParalutionLSS,
               LinearSystemSolver,
               ParalutionModule,
               1>
paralutionLSSMethodProvider("PARALUTION");

//////////////////////////////////////////////////////////////////////////////

void ParalutionLSS::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::string >
    ( "SetupCom", "Setup Command to run. This command seldomly needs overriding." );
  options.addConfigOption< std::string >
    ( "UnSetupCom", "UnSetup Command to run. This command seldomly needs overriding." );
  options.addConfigOption< std::string >
    ( "SysSolver", "Command that solves the linear system." );
}

//////////////////////////////////////////////////////////////////////////////

ParalutionLSS::ParalutionLSS(const std::string& name)
  : LinearSystemSolver(name)
{
  addConfigOptionsTo(this);

  m_data.reset(new ParalutionLSSData(getMaskArray(),
				     getNbSysEquations(),
				     this));
  
  cf_assert(m_data.getPtr() != CFNULL);
  
  m_setupStr = "StdSetup";
  setParameter("SetupCom",&m_setupStr);
  
  m_unSetupStr = "Null";
  setParameter("UnSetupCom",&m_unSetupStr);

  m_solveSysStr = "Null";
  setParameter("SysSolver",&m_solveSysStr);
}

//////////////////////////////////////////////////////////////////////////////

ParalutionLSS::~ParalutionLSS()
{
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<MethodData> ParalutionLSS::getMethodData () const
{
  return m_data.getPtr();
}

//////////////////////////////////////////////////////////////////////////////

void ParalutionLSS::configure ( Config::ConfigArgs& args )
{
  LinearSystemSolver::configure(args);
  configureNested ( m_data.getPtr(), args );
  
  // const bool isParallel = Common::PE::GetPE().IsParallel ();
  
  // if (m_setupStr == "Null") {
  //   m_setupStr = (isParallel) ? "NewParSetup" : "NewSeqSetup";
  // }
  
  // if (m_unSetupStr == "Null") {
  //   m_unSetupStr = (isParallel) ? "StdParUnSetup" : "StdSeqUnSetup";
  // }

  // if (m_solveSysStr == "Null") {
  //   m_solveSysStr = (isParallel) ? "StdParSolveSys" : "StdSeqSolveSys";
  // }

  configureCommand<ParalutionLSSData,ParalutionLSSComProvider>( args, m_setup,m_setupStr,m_data);

  configureCommand<ParalutionLSSData,ParalutionLSSComProvider>( args, m_unSetup,m_unSetupStr,m_data);
  
  configureCommand<ParalutionLSSData,ParalutionLSSComProvider>( args, m_solveSys,m_solveSysStr,m_data);
}

//////////////////////////////////////////////////////////////////////////////

void ParalutionLSS::solveSysImpl()
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
  
  CFLog(VERBOSE, "ParalutionLSS::solveSys() WallTime: " << stopTimer << "s\n");
}

//////////////////////////////////////////////////////////////////////////////

BlockAccumulator* ParalutionLSS::
createBlockAccumulator(const CFuint nbRows,
                       const CFuint nbCols,
                       const CFuint subBlockSize,
		       CFreal* ptr) const
{
  return new BlockAccumulator
    (nbRows,nbCols,subBlockSize, m_lssData->getLocalToLocallyUpdatableMapping(), ptr);
}

//////////////////////////////////////////////////////////////////////////////

void ParalutionLSS::printToFile(const std::string prefix, const std::string suffix)
{
  cf_assert(isSetup());
  cf_assert(isConfigured());

  ParalutionMatrix& mat = m_data->getMatrix();
  ParalutionVector& rhs = m_data->getRhsVector();
  ParalutionVector& sol = m_data->getSolVector();
  
  std::string matStr = prefix + "mat" + suffix;
  std::string rhsStr = prefix + "rhs" + suffix;
  std::string solStr = prefix + "sol" + suffix;
  
  mat.printToFile(matStr.c_str());
  rhs.printToFile(rhsStr.c_str());
  sol.printToFile(solStr.c_str());
}

//////////////////////////////////////////////////////////////////////////////

void ParalutionLSS::setMethodImpl()
{
  LinearSystemSolver::setMethodImpl();

  //  setupCommandsAndStrategies();
  
  m_setup->setup();
  m_setup->execute();
  
  m_solveSys->setup();
  m_unSetup->setup();
}

//////////////////////////////////////////////////////////////////////////////

void ParalutionLSS::unsetMethodImpl()
{
  m_unSetup->execute();
  unsetupCommandsAndStrategies();
  
  LinearSystemSolver::unsetMethodImpl();
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<NumericalStrategy> > ParalutionLSS::getStrategyList() const
{
  vector<Common::SafePtr<NumericalStrategy> > result;
  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace Paralution

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

