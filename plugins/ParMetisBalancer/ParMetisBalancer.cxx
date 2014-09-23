#include "Environment/ObjectProvider.hh"

#include "Common/BadValueException.hh"

#include "ParMetisBalancer/ParMetisBalancer.hh"
#include "ParMetisBalancer/ParMetisBalancerModule.hh"
#include "ParMetisBalancer/ParMetisBalancerData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace ParMetisBalancer {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<ParMetisBalancer,
               DynamicBalancerMethod,
               ParMetisBalancerModule,
               1>
parMetisBalancerDynamicBalancerMethodProvider("ParMETIS_Based");

//////////////////////////////////////////////////////////////////////////////

void ParMetisBalancer::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::string > ("SetupCom","SetupCommand to run. This command seldomly needs overriding.");
   options.addConfigOption< std::string > ("UnSetupCom","UnSetupCommand to run. This command seldomly needs overriding.");
   options.addConfigOption< std::string > ("StdRepart","UnSetupCommand to run. This command seldomly needs overriding.");
}

//////////////////////////////////////////////////////////////////////////////

ParMetisBalancer::ParMetisBalancer(const std::string& name)
  : DynamicBalancerMethod(name)
{
  addConfigOptionsTo(this);
  m_data.reset(new ParMetisBalancerData(this));
  cf_assert(m_data.isNotNull());

  m_setupStr = "StdSetup";
  setParameter("SetupCom",&m_setupStr);

  m_unSetupStr = "StdUnSetup";
  setParameter("UnSetupCom",&m_unSetupStr);

  m_repartStr = "StdRepart";
  setParameter("StdRepart",&m_repartStr);
}

//////////////////////////////////////////////////////////////////////////////

ParMetisBalancer::~ParMetisBalancer()
{
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<MethodData> ParMetisBalancer::getMethodData() const
{
  return m_data.getPtr();
}

//////////////////////////////////////////////////////////////////////////////

void ParMetisBalancer::doDynamicBalanceImpl()
{
  m_repart->execute();
}

//////////////////////////////////////////////////////////////////////////////

void ParMetisBalancer::configure(Config::ConfigArgs& args)
{
  CFAUTOTRACE;

  DynamicBalancerMethod::configure(args);
  configureNested(m_data.getPtr(), args);

  /// commands are configured here
  configureCommand< ParMetisBalancerData,ParMetisBalancerComProvider >(args, m_setup  , m_setupStr  , m_data );
  configureCommand< ParMetisBalancerData,ParMetisBalancerComProvider >(args, m_unSetup, m_unSetupStr,m_data );
  configureCommand< ParMetisBalancerData,ParMetisBalancerComProvider >(args, m_repart , m_repartStr , m_data );
}

//////////////////////////////////////////////////////////////////////////////

void ParMetisBalancer::setMethodImpl()
{
  DynamicBalancerMethod::setMethodImpl();

  setupCommandsAndStrategies();
  cf_assert(m_setup.isNotNull());
  m_setup->execute();
}

//////////////////////////////////////////////////////////////////////////////

void ParMetisBalancer::unsetMethodImpl()
{
  cf_assert(m_unSetup.isNotNull());
  m_unSetup->execute();
  unsetupCommandsAndStrategies();

  DynamicBalancerMethod::unsetMethodImpl();
}

//////////////////////////////////////////////////////////////////////////////


    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
