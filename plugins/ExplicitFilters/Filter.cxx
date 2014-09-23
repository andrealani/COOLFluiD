#include "Filter.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/Method.hh"
//#include "Utils/Event.hh"
#include "Framework/SubSystemStatus.hh"
#include "ExplicitFilters/ExplicitFilters.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

namespace COOLFluiD {
  using namespace Common;
  using namespace Framework;

  namespace Numerics {

    namespace ExplicitFilters {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<Filter,
               DataProcessingMethod,
               ExplicitFiltersModule,
               1>
               filterMethodProvider("Filter");

//////////////////////////////////////////////////////////////////////////////

void Filter::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption<std::string>("PrepareCmd","Preparation command for the filter process commands.");
   options.addConfigOption<std::vector<std::string> >("ProcessComds","Types of the filter process commands.");
   options.addConfigOption<std::vector<std::string> >("ProcessNames","Names of the filter process commands.");
}

//////////////////////////////////////////////////////////////////////////////

Filter::Filter(const std::string& name)
  : DataProcessingMethod(name),
    m_prepare(),
    m_processes(),
    m_processTypeStr(),
    m_processNameStr()
{
   addConfigOptionsTo(this);

  m_data.reset(new FilterData(this));

  m_prepareTypeStr = "Prepare";
  setParameter("PrepareCmd",&m_prepareTypeStr);

  m_processTypeStr = vector<std::string>();
  setParameter("ProcessComds",&m_processTypeStr);

  m_processNameStr = vector<std::string>();
  setParameter("ProcessNames",&m_processNameStr);

}

//////////////////////////////////////////////////////////////////////////////

Filter::~Filter()
{
  // clearPrepareComs();
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<Framework::MethodData> Filter::getMethodData() const
{
  return m_data.getPtr();
}

//////////////////////////////////////////////////////////////////////////////
//
// void Filter::clearPrepareComs()
// {
//   // _prepares.clear();
//   // vector<SelfRegistPtr<FilterCom> >().swap(_prepares);
// }
//

//////////////////////////////////////////////////////////////////////////////

void Filter::configure ( Config::ConfigArgs& args )
{
  
  CFLog(INFO, "Configuring Filter... \n");
  
  Method::configure(args);
  configureNested ( m_data.getPtr(), args );

  configureCommand<FilterCom,FilterData,FilterComProvider>(args,m_prepare,m_prepareTypeStr,m_prepareTypeStr,m_data);

  CFLog(INFO, "Configuring FilterCommands (" << m_processTypeStr.size() <<")  \n");
  cf_assert(m_processTypeStr.size() == m_processNameStr.size());
  m_processes.resize(m_processTypeStr.size());
  for(CFuint i=0; i<m_processes.size(); ++i) {
    CFLog(INFO, "COMMAND type = " << m_processTypeStr[i] << "\n");
    CFLog(INFO, "COMMAND name = " << m_processNameStr[i] << "\n");
    configureCommand<FilterCom,
                     FilterData,
                     FilterComProvider>(args,
                                m_processes[i],
                                m_processTypeStr[i],
                                m_processNameStr[i],
                                m_data);
    cf_assert(m_processes[i].isNotNull());
  }
  
  

}

//////////////////////////////////////////////////////////////////////////////

void Filter::setMethodImpl()
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

void Filter::unsetMethodImpl()
{
  CFAUTOTRACE;

  CFLog(INFO, "unsetup filter commands and strategies \n");
  unsetupCommandsAndStrategies();
  
  DataProcessingMethod::unsetMethodImpl();
  
  CFLog(INFO, "unsetup finished \n");
}

//////////////////////////////////////////////////////////////////////////////

void Filter::processDataImpl()
{

  CFAUTOTRACE;

  cf_assert(isConvergenceMethodSet());

  // The prepare command should only be run once.
  // The iteration should be checked inside the prepare command
  m_prepare->execute();

  // execute all processes.
  // Each command should have a built-in processRate which can be configured
  for(CFuint i=0; i<m_processes.size(); ++i) {
    cf_assert(m_processes[i].isNotNull());
    m_processes[i]->execute();
  }

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace ExplicitFilters

  } // namespace Numerics

} // namespace COOLFluiD

