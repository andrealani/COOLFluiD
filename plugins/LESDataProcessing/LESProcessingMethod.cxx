#include "LESDataProcessing.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/Method.hh"
//#include "Utils/Event.hh"
#include "Framework/SubSystemStatus.hh"
#include "LESProcessingMethod.hh"
#include "Framework/DataProcessingMethod.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

namespace COOLFluiD {
  using namespace Common;
  using namespace Framework;

  namespace Numerics {

    namespace LESDataProcessing {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<LESProcessing,
               Framework::DataProcessingMethod,
               LESDataProcessingModule,
               1>
               LESProcessingMethodProvider("LESDataProcessing");

//////////////////////////////////////////////////////////////////////////////

void LESProcessing::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption<std::vector<std::string> >("ProcessComds","Types of the LESProcessing process commands.");
}

//////////////////////////////////////////////////////////////////////////////

LESProcessing::LESProcessing(const std::string& name)
  : DataProcessingMethod(name),
    m_prepare(),
    m_processes(),
    m_processTypeStr()
{
   addConfigOptionsTo(this);

  m_data.reset(new LESProcessingData(this));

  m_processTypeStr = vector<std::string>();
  setParameter("ProcessComds",&m_processTypeStr);
}

//////////////////////////////////////////////////////////////////////////////

LESProcessing::~LESProcessing()
{
  // clearPrepareComs();
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<Framework::MethodData> LESProcessing::getMethodData() const
{
  return m_data.getPtr();
}

//////////////////////////////////////////////////////////////////////////////

void LESProcessing::configure ( Config::ConfigArgs& args )
{  
  Method::configure(args);
  
  CFLog(INFO," +++ Data     " << m_data->getClassName() << "::configure() \n");
  configureNested ( m_data.getPtr(), args );

  m_processes.resize(m_processTypeStr.size());
  for(CFuint i=0; i<m_processes.size(); ++i) {
    CFLog(INFO, " +++ Command  " << m_processTypeStr[i] << "::configure() \n");
    configureCommand<LESProcessingComBase,
                     LESProcessingData,
                     LESProcessingComProvider>(args,
                                m_processes[i],
                                m_processTypeStr[i],
                                m_processTypeStr[i],
                                m_data);
    cf_assert(m_processes[i].isNotNull());
  }
  
  

}

//////////////////////////////////////////////////////////////////////////////

void LESProcessing::setMethodImpl()
{
  CFAUTOTRACE;

  DataProcessingMethod::setMethodImpl();

  // The commands in this method have the ProcessRate. This Method-ProcessRate
  // must therefore be set to 1, and allow the commands to decide.
  m_processRate = 1;

  // Setup
  setupCommandsAndStrategies();
}

//////////////////////////////////////////////////////////////////////////////

void LESProcessing::unsetMethodImpl()
{
  CFAUTOTRACE;

  unsetupCommandsAndStrategies();
  
  DataProcessingMethod::unsetMethodImpl();
}

//////////////////////////////////////////////////////////////////////////////

void LESProcessing::processDataImpl()
{
  CFAUTOTRACE;

  cf_assert(isConvergenceMethodSet());

  // The prepare command should only be run once.
  // The iteration should be checked inside the prepare command
  // m_prepare->execute();

  // execute all processes.
  // Each command should have a built-in processRate which can be configured
  for(CFuint i=0; i<m_processes.size(); ++i) {
    cf_assert(m_processes[i].isNotNull());
    m_processes[i]->execute();
  }

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace LESProcessing

  } // namespace Numerics

} // namespace COOLFluiD

