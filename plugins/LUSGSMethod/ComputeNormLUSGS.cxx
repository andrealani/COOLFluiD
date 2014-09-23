#include "Framework/PhysicalModel.hh"

#include "LUSGSMethod/ComputeNormLUSGS.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace LUSGSMethod {

//////////////////////////////////////////////////////////////////////////////

ComputeNormLUSGS::ComputeNormLUSGS(const std::string & name) :
  ComputeNorm(name),
  socket_statesSetIdx("statesSetIdx"),
  socket_statesSetStateIDs("statesSetStateIDs"),
  socket_isStatesSetParUpdatable("isStatesSetParUpdatable"),
  m_localResiduals(),
  m_nbrEqs()
{
}

//////////////////////////////////////////////////////////////////////////////

ComputeNormLUSGS::~ComputeNormLUSGS()
{
}

//////////////////////////////////////////////////////////////////////////////

void ComputeNormLUSGS::setup()
{
  CFAUTOTRACE;

  ComputeNorm::setup();

  m_localResiduals.resize(m_residuals.size());

  m_nbrEqs = PhysicalModelStack::getActive()->getNbEq();
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > ComputeNormLUSGS::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result = ComputeNorm::needsSockets();

  result.push_back(&socket_statesSetIdx);
  result.push_back(&socket_statesSetStateIDs);
  result.push_back(&socket_isStatesSetParUpdatable);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace LUSGSMethod

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
