#include "LUSGSMethod/LUSGSMethod.hh"
#include "StdUnSetup.hh"
#include "Framework/MeshData.hh"
#include "Framework/State.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace LUSGSMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StdUnSetup, LUSGSIteratorData, LUSGSMethodModule> stdUnSetupProvider("StdUnSetup");

//////////////////////////////////////////////////////////////////////////////

StdUnSetup::StdUnSetup(std::string name) : LUSGSIteratorCom(name),
  socket_updateCoeff("updateCoeff"),
  socket_pastStates("pastStates"),
  socket_statesSetIdx("statesSetIdx"),
  socket_pivotLUFactorization("pivotLUFactorization")
{
}

//////////////////////////////////////////////////////////////////////////////

StdUnSetup::~StdUnSetup()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
StdUnSetup::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_updateCoeff);
  result.push_back(&socket_pastStates  );
  result.push_back(&socket_statesSetIdx);
  result.push_back(&socket_pivotLUFactorization);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void StdUnSetup::execute()
{
  CFAUTOTRACE;

  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
  updateCoeff.resize(0);

  DataHandle<State*> pastStates = socket_pastStates.getDataHandle();
  // deallocate the states in the storage
  for(CFuint i = 0; i < pastStates.size(); ++i) {
    delete pastStates[i];
  }
  pastStates.resize(0);

  DataHandle< CFint > statesSetIdx = socket_statesSetIdx.getDataHandle();
  statesSetIdx.resize(0);

  DataHandle< vector< CFuint > > pivotLUFactorization = socket_pivotLUFactorization.getDataHandle();
  pivotLUFactorization.resize(0);

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace LUSGSMethod

  } // namespace Numerics

} // namespace COOLFluiD
