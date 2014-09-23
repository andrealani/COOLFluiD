#include "AnalyticEE.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/SubSystemStatus.hh"
#include "AnalyticalEE/AnalyticalEE.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

namespace COOLFluiD {

  using namespace Framework;

    namespace AnalyticalEE {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<AnalyticEE,
               ErrorEstimatorMethod,
               AnalyticalEEModule,
               1>
analyticEEProvider("AnalyticalEE");

//////////////////////////////////////////////////////////////////////////////

void AnalyticEE::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::string >("SetupCom","Command to setup the method. This command seldomly needs overriding.");
   options.addConfigOption< std::string >("UnSetupCom","Command to unsetup the method. This command seldomly needs overriding.");
   options.addConfigOption< std::string >("ComputeCom","Command to compute the analytic error function.");
}

//////////////////////////////////////////////////////////////////////////////

AnalyticEE::AnalyticEE(const std::string& name)
  : ErrorEstimatorMethod(name)
{
   addConfigOptionsTo(this);

  m_data.reset( new AnalyticEEData(this) );

  m_setupStr = "StdSetup";
  setParameter("SetupCom",&m_setupStr);

  m_unSetupStr = "StdUnSetup";
  setParameter("UnSetupCom",&m_unSetupStr);

  m_computeStr = "ComputeDiscreteError";
  setParameter("ComputeCom",&m_computeStr);
}

//////////////////////////////////////////////////////////////////////////////

AnalyticEE::~AnalyticEE()
{
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<MethodData> AnalyticEE::getMethodData () const
{
  return m_data.getPtr();
}

//////////////////////////////////////////////////////////////////////////////

void AnalyticEE::configure ( Config::ConfigArgs& args )
{
  ErrorEstimatorMethod::configure(args);
  configureNested ( m_data.getPtr(), args );

  // add configures to the AnalyticEECom's

  configureCommand<AnalyticEEData,AnalyticEEComProvider>( args, m_setup,m_setupStr,m_data);

  configureCommand<AnalyticEEData,AnalyticEEComProvider>( args, m_unSetup,m_unSetupStr,m_data);

  configureCommand<AnalyticEEData,AnalyticEEComProvider>( args, m_compute,m_computeStr,m_data);

}

//////////////////////////////////////////////////////////////////////////////

void AnalyticEE::setMethodImpl()
{
  ErrorEstimatorMethod::setMethodImpl();

  setupCommandsAndStrategies();
  m_setup->execute();
}

//////////////////////////////////////////////////////////////////////////////

void AnalyticEE::unsetMethodImpl()
{
  m_unSetup->execute();
  unsetupCommandsAndStrategies();

  ErrorEstimatorMethod::unsetMethodImpl();
}

//////////////////////////////////////////////////////////////////////////////

void AnalyticEE::estimateImpl()
{

  CFAUTOTRACE;

  if(!(SubSystemStatusStack::getActive()->getNbIter() % _estimateRate))
  {
    // set the nodal states at first
    getMethodData()->getCollaborator<SpaceMethod>()->extrapolateStatesToNodes();

    m_compute->execute();
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace AnalyticalEE

} // namespace COOLFluiD

