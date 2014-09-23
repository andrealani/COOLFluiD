////////////////////////////////////////////////////////////////////////////////
//#include "Common/OwnedObject.hh"
//#include "Common/SetupObject.hh"
//#include "Common/NonCopyable.hh"
//#include "Environment/ConcreteProvider.hh"
//#include "RadiativeTransfer/RadiativeTransferModule.hh"
//#include "RadiationDistribution.hh"
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////

//namespace COOLFluiD {

//namespace RadiativeTransfer {

////////////////////////////////////////////////////////////////////////////////

//class RadiationPhysicsProvider : public Common::OwnedObject,
//                                 public Config::ConfigObject,
//                                 public Common::NonCopyable<RadiationPhysicsProvider>

//{
//public:
//    static std::string getClassName() { return "RadiationPhysicsProvider"; }

//    void configure(Config::ConfigArgs& args){}
//    static void defineConfigOptions(Config::OptionList& options){}

//    typedef Environment::ConcreteProvider<RadiationPhysicsProvider,1> PROVIDER;
//    typedef const std::string& ARG1;

//    RadiationPhysicsProvider(const std::string& name);
//    ~RadiationPhysicsProvider(){}
//};

////////////////////////////////////////////////////////////////////////////////


//RadiationPhysicsProvider::RadiationPhysicsProvider(const std::string& name):
//    Common::OwnedObject(),
//    ConfigObject(name)
//{
//  addConfigOptionsTo(this);
//}

