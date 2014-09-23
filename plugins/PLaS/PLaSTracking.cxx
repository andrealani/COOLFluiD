
#include "Environment/ObjectProvider.hh"
#include "PLaS/StgImplementation.hh"
#include "PLaS/PLaSModule.hh"
#include "PLaS/PLaSTracking.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;

namespace COOLFluiD {
  namespace PLaS {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider< PLaSTracking,DataProcessingMethod,PLaSModule,1 > plasTrackingProvider("PLaSTracking");

//////////////////////////////////////////////////////////////////////////////

void PLaSTracking::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::string >("SetupCom","Setup command (default \"StdSetup\")");
  options.addConfigOption< std::string >("ProcessCom","Process command (default \"StdProcess\")");
  options.addConfigOption< std::string >("UnSetupCom","UnSetup command (default \"StdUnSetup\")");
}

//////////////////////////////////////////////////////////////////////////////

PLaSTracking::PLaSTracking(const std::string& name) : DataProcessingMethod(name)
{
  addConfigOptionsTo(this);

  m_data.reset(new PLaSTrackingData(this));
  cf_assert(m_data.getPtr() != CFNULL);

  // commands
  m_setupStr   = "StdSetup";
  m_processStr = "StdProcess";
  m_unsetupStr = "StdUnSetup";
  setParameter("SetupCom",&m_setupStr);
  setParameter("ProcessCom",&m_processStr);
  setParameter("UnSetupCom",&m_unsetupStr);
}

//////////////////////////////////////////////////////////////////////////////

PLaSTracking::~PLaSTracking()
{
}

//////////////////////////////////////////////////////////////////////////////

void PLaSTracking::configure(Config::ConfigArgs& args)
{
  CFAUTOTRACE;

  // configure parent and data
  DataProcessingMethod::configure(args);

  configureNested(m_data.getPtr(),args);

  // configure setup, unsetup and process commands
  configureCommand< PLaSTrackingCom,PLaSTrackingData,PLaSTrackingComProvider >(args,m_setup,m_setupStr,m_setupStr,m_data);
  configureCommand< PLaSTrackingCom,PLaSTrackingData,PLaSTrackingComProvider >(args,m_process,m_processStr,m_processStr,m_data);
  configureCommand< PLaSTrackingCom,PLaSTrackingData,PLaSTrackingComProvider >(args,m_unsetup,m_unsetupStr,m_unsetupStr,m_data);
  CFLog(INFO,"PLaS command: " << m_setup->getName()   << "\n");
  CFLog(INFO,"PLaS command: " << m_process->getName() << "\n");
  CFLog(INFO,"PLaS command: " << m_unsetup->getName() << "\n");
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace PLaS
} // namespace COOLFluiD

