#include "LUSGSMethod/LUSGSMethod.hh"
#include "LUSGSMethod/BDF3Setup.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MeshData.hh"
#include "Framework/State.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace LUSGSMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<BDF3Setup, LUSGSIteratorData, LUSGSMethodModule> bdf3SetupProvider("BDF3Setup");

//////////////////////////////////////////////////////////////////////////////

BDF3Setup::BDF3Setup(const std::string& name) :
  StdSetup(name),
  socket_pastTimeRhs("pastTimeRhs"),
  socket_pastPastTimeRhs("pastPastTimeRhs")
{
}

//////////////////////////////////////////////////////////////////////////////

BDF3Setup::~BDF3Setup()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
BDF3Setup::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result = StdSetup::providesSockets();

  result.push_back(&socket_pastTimeRhs);
  result.push_back(&socket_pastPastTimeRhs);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > BDF3Setup::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result = StdSetup::needsSockets();

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void BDF3Setup::execute()
{
  CFAUTOTRACE;

  // execute of parent class
  StdSetup::execute();

  // get the states datahandle
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

  // get number of states and number of equations
  const CFuint nbStates = states.size();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

  // resize pastRhs
  DataHandle<CFreal> pastTimeRhs = socket_pastTimeRhs.getDataHandle();
  pastTimeRhs.resize(nbStates*nbEqs);

  // resize pastPastRhs
  DataHandle<CFreal> pastPastTimeRhs = socket_pastPastTimeRhs.getDataHandle();
  pastPastTimeRhs.resize(nbStates*nbEqs);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace LUSGSMethod

  } // namespace Numerics

} // namespace COOLFluiD
