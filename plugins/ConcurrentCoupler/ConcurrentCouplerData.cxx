#include "Framework/MethodCommandProvider.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/ConvectiveVarSet.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/NamespaceSwitcher.hh"

#include "ConcurrentCoupler/ConcurrentCoupler.hh"
#include "ConcurrentCoupler/ConcurrentCouplerData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace ConcurrentCoupler {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<NullMethodCommand<ConcurrentCouplerData>, 
		      ConcurrentCouplerData, ConcurrentCouplerModule> 
nullCouplerComProvider("Null");

//////////////////////////////////////////////////////////////////////////////

void ConcurrentCouplerData::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< string >
    ("SendToRecvVariableTransformer", "Transform send to recv variables.");
}
      
//////////////////////////////////////////////////////////////////////////////

ConcurrentCouplerData::ConcurrentCouplerData(SafePtr<Method> owner) : 
  CouplerData(owner),
  _stdTrsGeoBuilder(),
  _faceTrsGeoBuilder(),
  _spaceMethod(),
  _sendToRecvVecTrans()
{
  addConfigOptionsTo(this);
  
  _sendToRecvVecTransStr = "Identity";
  setParameter("SendToRecvVariableTransformer",&_sendToRecvVecTransStr);
}      

//////////////////////////////////////////////////////////////////////////////

ConcurrentCouplerData::~ConcurrentCouplerData()
{
}

//////////////////////////////////////////////////////////////////////////////

void ConcurrentCouplerData::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;
  
  CouplerData::configure(args);
  
  const string name = getNamespace();
  SafePtr<Namespace> nsp = 
    NamespaceSwitcher::getInstance(SubSystemStatusStack::getCurrentName()).getNamespace(name);
  SafePtr<PhysicalModel> physModel = PhysicalModelStack::getInstance().getEntryByNamespace(nsp);
  SafePtr<VarSetTransformer::PROVIDER> vecTransProv = CFNULL;
  
  CFLog(VERBOSE, "Configuring VarSet Transformer: " << _sendToRecvVecTransStr << "\n");
  
  try {
    vecTransProv = Environment::Factory<VarSetTransformer>::getInstance().getProvider
      (_sendToRecvVecTransStr);
  }
  catch (Common::NoSuchValueException& e) {
    _sendToRecvVecTransStr = "Identity";
    
    CFLog(VERBOSE, e.what() << "\n");
    CFLog(VERBOSE, "Choosing IdentityVarSetTransformer instead ..." << "\n");
    vecTransProv = Environment::Factory<VarSetTransformer>::getInstance().getProvider
      (_sendToRecvVecTransStr);
  }
  
  cf_assert(vecTransProv.isNotNull());
  _sendToRecvVecTrans.reset(vecTransProv->create(physModel->getImplementor()));
  cf_assert(_sendToRecvVecTrans.getPtr() != CFNULL);
}
      
//////////////////////////////////////////////////////////////////////////////

void ConcurrentCouplerData::setup()
{
  // set up the GeometricEntity builders
  _stdTrsGeoBuilder.setup();
  _faceTrsGeoBuilder.setup();
  
  const string name = getNamespace();
  SafePtr<Namespace> nsp = 
    NamespaceSwitcher::getInstance(SubSystemStatusStack::getCurrentName()).getNamespace(name);
  SafePtr<PhysicalModel> physModel = PhysicalModelStack::getInstance().getEntryByNamespace(nsp);
  
  // set up the variable transformer
  _sendToRecvVecTrans->setup(1);
}
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace SubSystem

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
