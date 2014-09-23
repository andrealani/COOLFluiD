// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "TrilinosLSS.hh"
#include "Environment/ObjectProvider.hh"
#include "Common/PE.hh"
#include "Trilinos/Trilinos.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace Trilinos {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<TrilinosLSS,
         LinearSystemSolver,
               TrilinosModule,
         1>
trilinosLSSMethodProvider("TRILINOS");

//////////////////////////////////////////////////////////////////////////////

void TrilinosLSS::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::string >("SetupCom","Setup Command to run. This command seldomly needs overriding.");
   options.addConfigOption< std::string >("UnSetupCom","UnSetup Command to run. This command seldomly needs overriding.");
   options.addConfigOption< std::string >("SysSolver","Command that solves the linear system.");
}

//////////////////////////////////////////////////////////////////////////////

TrilinosLSS::TrilinosLSS(const std::string& name)
  : LinearSystemSolver(name)
{
   addConfigOptionsTo(this);
/* no matching constructor
  m_data.reset(new TrilinosLSSData(getLocalToGlobalMapping(),
				  getMaskArray(),
				  getNbSysEquations(),
				  this));
*/
  m_data.reset(new TrilinosLSSData(getMaskArray(),
				   getNbSysEquations(),
				   this));

  cf_assert(m_data.getPtr() != CFNULL);

  const bool isParallel = Common::PE::GetPE().IsParallel ();


  m_setupStr = "StdParSetup";
  setParameter("SetupCom",&m_setupStr);

  m_unSetupStr = "StdParUnSetup";
  setParameter("UnSetupCom",&m_unSetupStr);

  m_solveSysStr = "StdParSolveSys";
  setParameter("SysSolver",&m_solveSysStr);

}

//////////////////////////////////////////////////////////////////////////////

TrilinosLSS::~TrilinosLSS()
{
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<MethodData> TrilinosLSS::getMethodData() const
{
  return m_data.getPtr();
}

//////////////////////////////////////////////////////////////////////////////

void TrilinosLSS::configure ( Config::ConfigArgs& args )
{
  LinearSystemSolver::configure(args);
  configureNested ( m_data.getPtr(), args );

  // add here configures to the TrilinosLSS
  configureCommand<TrilinosLSSData,TrilinosLSSComProvider>(args,m_setup,m_setupStr,m_data);
  configureCommand<TrilinosLSSData,TrilinosLSSComProvider>(args,m_unSetup,m_unSetupStr,m_data);
  configureCommand<TrilinosLSSData,TrilinosLSSComProvider>(args,m_solveSys,m_solveSysStr,m_data);
}

//////////////////////////////////////////////////////////////////////////////

void TrilinosLSS::solveSysImpl()
{
  m_solveSys->execute();
}

//////////////////////////////////////////////////////////////////////////////

void TrilinosLSS::printToFile(const std::string prefix, const std::string suffix)
{
  cf_assert(isSetup());
  cf_assert(isConfigured());

  TrilinosMatrix& mat = *m_data->getMatrix();
  TrilinosVector& rhs = *m_data->getRhsVector();
  TrilinosVector& sol = *m_data->getSolVector();

  std::string matStr = prefix + "mat" + suffix;
  std::string rhsStr = prefix + "rhs" + suffix;
  std::string solStr = prefix + "sol" + suffix;

  mat.printToFile(matStr.c_str());
  rhs.printToFile(rhsStr.c_str());
  sol.printToFile(solStr.c_str());
}

//////////////////////////////////////////////////////////////////////////////

void TrilinosLSS::setMethodImpl()
{
  LinearSystemSolver::setMethodImpl();

  m_setup->setup();
  m_setup->execute();

  m_solveSys->setup();
  m_unSetup->setup();
}

//////////////////////////////////////////////////////////////////////////////

void TrilinosLSS::unsetMethodImpl()
{
  m_unSetup->execute();
  unsetupCommandsAndStrategies();

  LinearSystemSolver::unsetMethodImpl();
}

//////////////////////////////////////////////////////////////////////////////

Framework::BlockAccumulator* TrilinosLSS::
createBlockAccumulator(const CFuint nbRows,
                       const CFuint nbCols,
                       const CFuint subBlockSize,  
                       CFreal* ptr) const
{
  Framework::BlockAccumulator *accu =
    new Framework::BlockAccumulator(nbRows,
                                    nbCols,
                                    subBlockSize,
                                    m_lssData->getLocalToGlobalMapping(),
                                    ptr);
  accu->getWorkspaceNC().resize((nbRows+2*nbCols)*subBlockSize);
  return accu;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace Trilinos

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
