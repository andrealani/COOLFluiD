#include "LUSGSMethod/LUSGSMethod.hh"
#include "LUSGSMethod/BDF2Setup.hh"
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

MethodCommandProvider<BDF2Setup, LUSGSIteratorData, LUSGSMethodModule> bdf2SetupProvider("BDF2Setup");

//////////////////////////////////////////////////////////////////////////////

BDF2Setup::BDF2Setup(const std::string& name) :
  StdSetup(name),
  socket_pastTimeRhs("pastTimeRhs")
{
}

//////////////////////////////////////////////////////////////////////////////

BDF2Setup::~BDF2Setup()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
BDF2Setup::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result = StdSetup::providesSockets();

  result.push_back(&socket_pastTimeRhs);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > BDF2Setup::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result = StdSetup::needsSockets();

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void BDF2Setup::execute()
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
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace LUSGSMethod

  } // namespace Numerics

} // namespace COOLFluiD
