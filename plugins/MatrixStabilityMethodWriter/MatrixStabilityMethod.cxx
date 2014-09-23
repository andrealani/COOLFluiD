#include "Environment/ObjectProvider.hh"
#include "Environment/DirPaths.hh"
#include "Framework/SpaceMethod.hh"

#include "MatrixStabilityMethodWriter/MatrixStabilityMethodWriter.hh"
#include "MatrixStabilityMethodWriter/MatrixStabilityMethod.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace MatrixStabilityMethodWriter {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<MatrixStabilityMethod,
               ConvergenceMethod,
               MatrixStabilityMethodWriterModule,
               1>
MatrixStabilityMethodProvider("MatrixStabilityMethod");

//////////////////////////////////////////////////////////////////////////////

void MatrixStabilityMethod::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::string >("SetupCom"             ,"SetupCommand to run. This command seldomly needs overriding.");
  options.addConfigOption< std::string >("UnSetupCom"           ,"UnSetupCommand to run. This command seldomly needs overriding.");
  options.addConfigOption< std::string >("SetStates"            ,"SetStates command to run. This command seldomly needs overriding.");
  options.addConfigOption< std::string >("AddMatrixColumnToFile","AddMatrixColumnToFile SetStates command to run. This command seldomly needs overriding.");
}

//////////////////////////////////////////////////////////////////////////////

MatrixStabilityMethod::MatrixStabilityMethod(const std::string& name) :
  ConvergenceMethod(name)
{
   addConfigOptionsTo(this);

   m_data.reset(new MatrixStabilityMethodData(this));

   m_setupStr = "StdSetup";
   setParameter("SetupCom",&m_setupStr);

   m_unSetupStr = "StdUnSetup";
   setParameter("UnSetupCom",&m_unSetupStr);

   m_setStatesStr = "SetStates";
   setParameter("SetStates",&m_setStatesStr);

   m_addMatrixColumnToFileStr = "AddMatrixColumnToFile";
   setParameter("AddMatrixColumnToFile",&m_addMatrixColumnToFileStr);
}

//////////////////////////////////////////////////////////////////////////////

MatrixStabilityMethod::~MatrixStabilityMethod()
{
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<MethodData> MatrixStabilityMethod::getMethodData() const
{
  return m_data.getPtr();
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<ConvergenceMethodData> MatrixStabilityMethod::getConvergenceMethodData()
{
  return m_data.getPtr();
}

//////////////////////////////////////////////////////////////////////////////

void MatrixStabilityMethod::configure ( Config::ConfigArgs& args )
{
  ConvergenceMethod::configure(args);
  configureNested ( m_data.getPtr(), args );

  // add here configures to the MatrixStabilityMethod's
  configureCommand<MatrixStabilityMethodData,MatrixStabilityMethodComProvider>( args, m_setup,m_setupStr,m_data);

  configureCommand<MatrixStabilityMethodData,MatrixStabilityMethodComProvider>( args, m_unSetup,m_unSetupStr,m_data);

  configureCommand<MatrixStabilityMethodData,MatrixStabilityMethodComProvider>( args, m_setStates,m_setStatesStr,m_data);

  configureCommand<MatrixStabilityMethodData,MatrixStabilityMethodComProvider>( args, m_addMatrixColumnToFile,
                                                                               m_addMatrixColumnToFileStr,m_data);
}

//////////////////////////////////////////////////////////////////////////////

void MatrixStabilityMethod::setMethodImpl()
{
  //call the parent
  ConvergenceMethod::setMethodImpl();

  setupCommandsAndStrategies();
  m_setup->execute();
}

//////////////////////////////////////////////////////////////////////////////

void MatrixStabilityMethod::unsetMethodImpl()
{
  m_unSetup->execute();
  unsetupCommandsAndStrategies();
  //call the parent
  ConvergenceMethod::unsetMethodImpl();
}

//////////////////////////////////////////////////////////////////////////////

void MatrixStabilityMethod::takeStepImpl()
{
  CFAUTOTRACE;

  if (m_data->getSetStateToZero())
  // command has already been executed --> more than one iteration, not necessary
  {
    CF_DEBUG_EXIT;
  }

  // get number of states
  const CFuint nbrStates = m_data->getNbrStates();

  // create file
  const std::string fileName =
      Environment::DirPaths::getInstance().getResultsDir().string() + "/"
      + m_data->getOutputFileName();
  SafePtr< ofstream > outputFile = m_data->getOutputFile();
  if (m_data->writeBinary())
  {
    outputFile->open(fileName.c_str(), ios::binary | ios::out | ios::trunc);
    outputFile->write(reinterpret_cast< const char* >(&nbrStates),sizeof(CFuint));
  }
  else
  {
    outputFile->open(fileName.c_str());
    *outputFile << nbrStates << "\n";
  }

  // set all states to zero
  m_data->setAllStatesToZero(true);
  m_setStates->execute();
  m_data->setAllStatesToZero(false);

  // loop over the states
  for (CFuint i = 0; i < nbrStates; ++i)
  {
    CFLog(INFO,"State idx: " << i << "\n");

    // set the current state to one
    m_data->setStateIdx(i);
    m_data->setStateToZero(false);
    m_setStates->execute();

    // Compute the RHS
    m_data->getCollaborator<SpaceMethod>()->prepareComputation();
    m_data->getCollaborator<SpaceMethod>()->computeSpaceResidual(1.0);
    m_data->getCollaborator<SpaceMethod>()->computeTimeResidual(1.0);

    // add column to the file
    m_addMatrixColumnToFile->execute();

    // reset the current state to zero
    m_data->setStateToZero(true);
    m_setStates->execute();
  }

  // close the file
  outputFile->close();
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace MatrixStabilityMethodWriter

  } // namespace Numerics

} // namespace COOLFluiD
