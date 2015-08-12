#include "PhysicalModelDummy/CouplingModelDummy.hh"
#include "PhysicalModelDummy/CouplingModelDummySendToRecv.hh"
#include "PhysicalModelDummy/Dummy.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace PhysicalModelDummy {
    
//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<CouplingModelDummySendToRecv, 
			    VarSetTransformer, 
			    DummyModule, 1> 
couplingModelDummySendToRecvProvider("CouplingModelDummySendToRecv");

//////////////////////////////////////////////////////////////////////////////

CouplingModelDummySendToRecv::CouplingModelDummySendToRecv
(SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetTransformer(model),
  m_model(model.d_castTo<CouplingModelDummy>())
{
}

//////////////////////////////////////////////////////////////////////////////

CouplingModelDummySendToRecv::~CouplingModelDummySendToRecv()
{
}

//////////////////////////////////////////////////////////////////////////////

void CouplingModelDummySendToRecv::transform(const State& state, State& result)
{
  result = 0.;
  const vector<CFuint>& sendIDs = m_model->getSendIDs();
  const vector<CFuint>& recvIDs = m_model->getRecvIDs();
  cf_assert(sendIDs.size() == recvIDs.size());
  
  for (CFuint i = 0; i < sendIDs.size(); ++i) {
    cf_assert(recvIDs[i] < result.size());
    cf_assert(sendIDs[i] < state.size());
    result[recvIDs[i]] = state[sendIDs[i]];
  }
}

//////////////////////////////////////////////////////////////////////////////

void CouplingModelDummySendToRecv::transformFromRef(const RealVector& data, State& result)
{
  throw NotImplementedException
    (FromHere(),"CouplingModelDummySendToRecv::transformFromRef() not implemented!");
}
    
//////////////////////////////////////////////////////////////////////////////

    } // namespace PhysicalModelDummy

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
