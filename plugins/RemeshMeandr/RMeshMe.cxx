#include "RMeshMe.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/SubSystemStatus.hh"
#include "RemeshMeandr/RemeshMeandr.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

namespace COOLFluiD {

  using namespace Framework;

  namespace Numerics {

    namespace RemeshMeandros {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider< RMeshMe, MeshAdapterMethod, RemeshMeandrModule, 1>  meandrosRemeshProvider("RemeshMeandros");

//////////////////////////////////////////////////////////////////////////////

void RMeshMe::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::string >("SetupCom","SetupCommand to run. This command seldomly needs overriding.");
   options.addConfigOption< std::string >("MeandrosCallCom","MeandrosCallCom to run. This command seldomly needs overriding.");
   options.addConfigOption< std::string >("WriteControlSpcCom","WriteControlSpc to run. This command seldomly needs overriding.");
   options.addConfigOption< std::string >("UnSetupCom","UnSetupCommand to run. This command seldomly needs overriding.");
   options.addConfigOption< std::string >("InterpolCom","InterpolCom to run. This command seldomly needs overriding.");
}

//////////////////////////////////////////////////////////////////////////////

RMeshMe::RMeshMe(const std::string& name)
  : MeshAdapterMethod(name)
{
   addConfigOptionsTo(this);

  _data.reset( new RMeshMeData(this) );

  _setupStr = "StdSetup";
   setParameter("SetupCom",&_setupStr);



  _unSetupStr = "StdUnSetup";
   setParameter("UnSetupCom",&_unSetupStr);



  _writeControlSpcStr = "StdWriteControlSpc";
   setParameter("WriteControlSpcCom",&_writeControlSpcStr);



  _meandrosCallStr = "StdMeandrosCall";
   setParameter("MeandrosCallCom",&_meandrosCallStr);



  _interpolStr = "StdInterpol";
   setParameter("InterpolCom",&_interpolStr);



}

//////////////////////////////////////////////////////////////////////////////

RMeshMe::~RMeshMe()
{
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<MethodData> RMeshMe::getMethodData() const
{
  return _data.getPtr();
}

//////////////////////////////////////////////////////////////////////////////

void RMeshMe::configure ( Config::ConfigArgs& args )
{
  ConfigObject::configure(args);
  configureNested ( _data.getPtr(), args );

  // add configures to the RMeshMeCom's

  configureCommand<RMeshMeData,RMeshMeComProvider>( args, _setup,_setupStr,_data);

  configureCommand<RMeshMeData,RMeshMeComProvider>( args, _unSetup,_unSetupStr,_data);

  configureCommand<RMeshMeData,RMeshMeComProvider>( args, _writeControlSpcCom,_writeControlSpcStr,_data);

  configureCommand<RMeshMeData,RMeshMeComProvider>( args, _meandrosCallCom,_meandrosCallStr,_data);

  configureCommand<RMeshMeData,RMeshMeComProvider>( args, _interpolCom,_interpolStr,_data);

}

//////////////////////////////////////////////////////////////////////////////

void RMeshMe::setMethodImpl()
{
  MeshAdapterMethod::setMethodImpl();

  setupCommandsAndStrategies();
  _setup->execute();
}

//////////////////////////////////////////////////////////////////////////////

void RMeshMe::unsetMethodImpl()
{
  _unSetup->execute();
  unsetupCommandsAndStrategies();

  MeshAdapterMethod::unsetMethodImpl();
}

//////////////////////////////////////////////////////////////////////////////

void RMeshMe::adaptMeshImpl()
{
  CFAUTOTRACE;

  if(!(SubSystemStatusStack::getActive()->getNbIter() % _adaptRate))
  {
    CFLog( ERROR,"HELLO " << " FUNCT: " <<  __FUNCTION__ << " LINE: " << __LINE__ << "\n" );

    _writeControlSpcCom->execute();

    _meandrosCallCom->execute();

    _interpolCom->execute();

    //if ( _data->isHessianSmooth() )
    //    _smoother->execute();

    //_computeMetric->execute();

    /// @todo Add this later: _smoother->execute();
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace RemeshMeandros

  } // namespace Numerics

} // namespace COOLFluiD

