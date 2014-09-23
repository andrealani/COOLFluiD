#include <fstream>
#include <sstream>

#include "Common/PE.hh"

#include "Common/ProcessInfo.hh"
#include "Common/OSystem.hh"
#include "Environment/FileHandlerOutput.hh"

#include "Environment/CFEnv.hh"
#include "Environment/DirPaths.hh"
#include "Environment/SingleBehaviorFactory.hh"

#include "Framework/PathAppender.hh"
#include "Framework/ConvergenceMethod.hh"
#include "Framework/ConvergenceMethodData.hh"

#include "Common/Stopwatch.hh"
#include "MathTools/MathConsts.hh"

#include "Environment/ObjectProvider.hh"
#include "Common/CFLog.hh"

#include "Framework/SubSystemStatus.hh"
#include "Framework/SpaceMethod.hh"
#include "Framework/CFL.hh"

#include "LUSGSMethod/LUSGSMethod.hh"
#include "LUSGSMethod/LUSGSBDF2.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace LUSGSMethod {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<LUSGSBDF2,
               ConvergenceMethod,
               LUSGSMethodModule,
               1>
lusgsBDF2Provider("NonlinearLUSGSBDF2");

//////////////////////////////////////////////////////////////////////////////

LUSGSBDF2::LUSGSBDF2(const std::string& name)
  : LUSGSIterator(name)
{
}

//////////////////////////////////////////////////////////////////////////////

LUSGSBDF2::~LUSGSBDF2()
{
}

//////////////////////////////////////////////////////////////////////////////

void LUSGSBDF2::prepare()
{
  CFAUTOTRACE;

  LUSGSIterator::prepare();
}

//////////////////////////////////////////////////////////////////////////////

void LUSGSBDF2::takeStepImpl()
{
  CFAUTOTRACE;

  // do the prepare step, usually backing up the solution to past states
  prepare();

  // get subsystemstatus
  Common::SafePtr<SubSystemStatus> subSysStatus = SubSystemStatusStack::getActive();

  // create a timer
  Common::Stopwatch<WallTime> timer;
  Common::Stopwatch<WallTime> total_timer;
  total_timer.restart();

  // get time steps
  const CFreal newDT = subSysStatus->getDT();
  const CFreal oldDT = subSysStatus->getPreviousDT();

  // compute factors related to BDF2 scheme
  const CFreal alpha = oldDT / newDT;
  CFreal xi = 1. / ( 1.+ alpha);
  CFreal theta = 1.;

  // for the first step, backward Euler is used
  if(subSysStatus->getNbIter() == 1)
  {
    xi = 0.;
    theta = 1.0;
  }

  // SET CURRENT TIME TO TIME AT N+1 (IN THE FOLLOWING ONLY (APPROXIMATIONS OF) THE RHS AT N+1 ARE COMPUTED
  subSysStatus->updateCurrentTime();

  // Compute the block Jacobian matrices at the first time-step and every m_JacobFreezFreq time-steps
  if((subSysStatus->getNbIter() == 1) || ((subSysStatus->getNbIter())%(m_data->getJacobFreezFreq(0))==0))
  {
    timer.restart();

    CFout << "Compute the block Jacobian matrices ...\n";

    // MAKE THE SPACE METHOD COMPUTE THE BLOCK JACOBIAN MATRICES.
    CFLog(VERBOSE,"Start computation of diagonal block Jacobian matrices\n");
    CFLogDebugMed( "LUSGSBDF2::takeStep(): preparing Computation\n");
    m_data->getCollaborator<SpaceMethod>()->prepareComputation();
    CFLogDebugMed( "LUSGSBDF2::takeStep(): computing Space Residual\n");
    m_data->getCollaborator<SpaceMethod>()->computeSpaceResidual(theta);
    CFLogDebugMed( "LUSGSBDF2::takeStep(): computing the Time Residual\n");
    m_data->getCollaborator<SpaceMethod>()->computeTimeResidual(xi);

    // FACTORIZES THE BLOCK JACOBIAN MATRICES.
    CFLog(VERBOSE,"Starting matrix factorizations\n");
    m_luFactorization->execute();

    // output factorization time
    CFLog(VERBOSE,"Computing and factorizing diagonal block matrices took: " << timer << "s\n");

  }

  CFout << "Starting LU-SGS iterator ...\n";

  // DO THE SGS SWEEPS
  CFuint k = 0;
  for(; k < m_data->getNbMaxSweeps(0) && !m_data->isAchieved(); ++k)
  {
    timer.restart();

    // set subsystem status
    subSysStatus->setFirstStep( k == 0 );
    subSysStatus->setMaxDT(MathTools::MathConsts::CFrealMax());

    // Do forward sweep
    CFLog(VERBOSE,"LUSGSBDF2::takeStep(): starting forward sweep\n");
    m_data->setForwardSweep(true);
    m_data->setStopSweep(false);
    // Update states set index (it is equal to -1 at this point)
    m_updateStatesSetIndex->execute();
    for (;!m_data->stopSweep();)
    {
      // Compute space residual for the current states set
      m_data->getCollaborator<SpaceMethod>()->computeSpaceRhsForStatesSet(theta);

      // Compute time residual for the current states set
      m_data->getCollaborator<SpaceMethod>()->computeTimeRhsForStatesSet(xi);

      // Compute the solution update for the current states set
      m_computeStatesSetUpdate->execute();

      // Update the solution for the current states set
      m_updateSol->execute();

      // Update states set index
      m_updateStatesSetIndex->execute();
    }

    // Syncronize the states
    ConvergenceMethod::syncGlobalDataComputeResidual(false);

    // reset the residual norms in the local processor
    m_data->getLUSGSNormComputer()->resetResiduals();

    // Do backward sweep
    CFLog(VERBOSE,"LUSGSBDF2::takeStep(): starting backward sweep\n");
    m_data->setForwardSweep(false);
    m_data->setStopSweep(false);
    // Update states set index (it is equal to the number of states sets at this point)
    m_updateStatesSetIndex->execute();
    for (;!m_data->stopSweep();)
    {
      // Compute space residual for the current states set
      m_data->getCollaborator<SpaceMethod>()->computeSpaceRhsForStatesSet(theta);

      // Compute time residual for the current states set
      m_data->getCollaborator<SpaceMethod>()->computeTimeRhsForStatesSet(xi);

      // Compute the solution update for the current states set
      m_computeStatesSetUpdate->execute();

      // Update the solution for the current states set
      m_updateSol->execute();

      // add contribution of current states set to the residual norms in the local processor
      m_data->getLUSGSNormComputer()->addStatesSetContribution();

      // Update states set index
      m_updateStatesSetIndex->execute();
    }

    // Syncronize the states
    ConvergenceMethod::syncGlobalDataComputeResidual(true);

    // output SGS time
    CFLog (VERBOSE, "SGS sweep took: " << timer << "s\n");

    // Gets the residual norm
    getConvergenceMethodData()->getConvergenceStatus().res = subSysStatus->getResidual();

    // Display info over each step of the Nonlinear LU-SGS iterator
    // The L2 norm of dU shown here is the global residual (spatial residual + time residual)
    if (m_data->isPrintHistory())
    {
      CFout << "Nonlinear LU-SGS sweep: " << k+1 << " L2 dU: " << getConvergenceMethodData()->getConvergenceStatus().res  << "\n";
    }
    
    // Simple Stop Condition (L2Norm reached)
    // HERE, if you stop because of residual, then you will do one step too much...
    m_data->setAchieved(subSysStatus->getResidual() < m_data->getMaxNorm());

    // set subsystem status
    subSysStatus->setFirstStep(false);
    subSysStatus->setLastStep(m_data->isAchieved());
  }

  // BACK UP THE TIME CONTRIBUTION
  m_data->setForwardSweep(true);
  m_data->setStopSweep(false);

  // Update states set index (it is equal to -1 at this point)
  m_updateStatesSetIndex->execute();
  for (;!m_data->stopSweep();)
  {
    // Back up the time rhs for the current states set
    m_data->getCollaborator<SpaceMethod>()->computeTimeRhsForStatesSet(1.0);

    // Update states set index
    m_updateStatesSetIndex->execute();
  }

  /// @todo MP: The following commands compute the spatial residuals for each computed variables. If n is the number of cpu used, then n files are written but their spatial residuals are exactly the same. This is a preliminary implementation.
  // BACK UP THE L2 NORMS OF THE ALL GLOBAL RESIDUALS
  RealVector tmpRes = subSysStatus->getAllResiduals();

  // COMPUTE THE SPACE RESIDUAL AND WRITE IT TO A FILE
  // Reset the residual norms in the local processor
  m_data->getLUSGSNormComputer()->resetResiduals();

  //Reset the global residual norms
  subSysStatus->resetResidual();

  m_data->setForwardSweep(false);
  m_data->setStopSweep(false);

  // Update states set index (it is equal to -1 at this point)
  m_updateStatesSetIndex->execute();
  for (;!m_data->stopSweep();)
  {
    // Compute space residual for the current states set
    m_data->getCollaborator<SpaceMethod>()->computeSpaceRhsForStatesSet(1.0);

    // Add contribution of current states set to the residual norms in the local processor
    m_data->getLUSGSNormComputer()->addStatesSetContribution();

    // Update states set index
    m_updateStatesSetIndex->execute();
  }

  // Compute all the spatial residual norms
  // Syncronize the states
  ConvergenceMethod::syncGlobalDataComputeResidual(true);

  // Open the Spatial-Residual.plt file in the results direcoty
  m_nameConvergenceFile = "Spatial-Residual.plt";
  fpath = m_nameConvergenceFile;
  fpath = Environment::DirPaths::getInstance().getResultsDir() /
          Framework::PathAppender::getInstance().appendParallel( fpath );

  SelfRegistPtr<Environment::FileHandlerOutput> fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();

  // Get the number of equations
  m_nbrEqs = PhysicalModelStack::getActive()->getNbEq();

  // Write Tecplot header at the first iteration
  if(subSysStatus->getNbIter() == 1)
  {
    ofstream& convergenceFile = fhandle->open(fpath);
    convergenceFile << "TITLE  =  Convergence of the spatial residuals" << "\n";
    convergenceFile << "VARIABLES = Iter";

    for (m_var_itr = 0; m_var_itr < tmpRes.size(); ++m_var_itr)
    {
      convergenceFile << " Res[" << m_var_itr << "]";
    }

    convergenceFile << "\n";
    fhandle->close();
  }

  ofstream& convergenceFile = fhandle->open(fpath,ios_base::app);

  // Write iteration index and the spatial L2 norm of the residuals of all the computed variables
  convergenceFile << subSysStatus->getNbIter() << " ";

  RealVector spatialRes = subSysStatus->getAllResiduals();

  for (m_var_itr = 0; m_var_itr < tmpRes.size(); ++m_var_itr)
  {
    //CFreal res = log10(sqrt(m_data->getLUSGSNormComputer()->getSpaceResValue(m_var_itr)));

    convergenceFile << spatialRes[m_var_itr] << " ";
  }

  convergenceFile << "\n";

  // Close the file Spatial-Residual.plt
  fhandle->close();

  // SET BACK ALL THE GLOBAL RESIDUAL NORMS
  subSysStatus->resetResidual();
  subSysStatus->setResidual(tmpRes);

  CFLog (VERBOSE, "Nonlinear LU-SGS convergence took: " << total_timer << "s\n");
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace LUSGSMethod

  } // namespace Numerics

} // namespace COOLFluiD
