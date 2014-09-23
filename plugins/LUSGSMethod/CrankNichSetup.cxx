#include "LUSGSMethod/LUSGSMethod.hh"
#include "LUSGSMethod/CrankNichSetup.hh"
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

MethodCommandProvider<CrankNichSetup, LUSGSIteratorData, LUSGSMethodModule> crankNichSetupProvider("CrankNichSetup");

//////////////////////////////////////////////////////////////////////////////

CrankNichSetup::CrankNichSetup(const std::string& name) :
  StdSetup(name),
  socket_pastRhs("pastRhs")
{
}

//////////////////////////////////////////////////////////////////////////////

CrankNichSetup::~CrankNichSetup()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
CrankNichSetup::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result = StdSetup::providesSockets();

  result.push_back(&socket_pastRhs);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > CrankNichSetup::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result = StdSetup::needsSockets();

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void CrankNichSetup::execute()
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
  DataHandle<CFreal> pastRhs = socket_pastRhs.getDataHandle();
  pastRhs.resize(nbStates*nbEqs);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace LUSGSMethod

  } // namespace Numerics

} // namespace COOLFluiD
