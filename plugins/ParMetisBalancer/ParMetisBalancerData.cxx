#include "Common/CFPrintContainer.hh"


#include "Framework/MethodCommandProvider.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Framework/OptionMethodStrategy.hh"

#include "ParMetisBalancer/ParMetisBalancer.hh"




#include "ParMetisBalancer/ParMetisBalancerModule.hh"
#include "ParMetisBalancer/ParMetisBalancerData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace ParMetisBalancer {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<NullMethodCommand<ParMetisBalancerData>,
		      ParMetisBalancerData, ParMetisBalancerModule>
aNullParMetisBalancerComProvider("Null");

//////////////////////////////////////////////////////////////////////////////

void ParMetisBalancerData::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

ParMetisBalancerData::ParMetisBalancerData(Common::SafePtr<Framework::Method> owner)
  : DynamicBalancerMethodData(owner)
{
  CFAUTOTRACE;

  addConfigOptionsTo(this);

}

//////////////////////////////////////////////////////////////////////////////

ParMetisBalancerData::~ParMetisBalancerData()
{
}

//////////////////////////////////////////////////////////////////////////////

void ParMetisBalancerData::configure(Config::ConfigArgs& args)
{
  CFAUTOTRACE;

  CFLog(VERBOSE,"Configuring ParMetisBalancerData\n");

  DynamicBalancerMethodData::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void ParMetisBalancerData::setup()
{
  CFAUTOTRACE;

  DynamicBalancerMethodData::setup();
}

//////////////////////////////////////////////////////////////////////////////

void ParMetisBalancerData::unsetup()
{
  DynamicBalancerMethodData::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

