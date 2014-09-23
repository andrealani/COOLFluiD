#include "FiniteVolume/FiniteVolume.hh"
#include "StdLinSetup.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StdLinSetup, CellCenterFVMData, FiniteVolumeModule> StdLinSetupProvider("StdLinSetup");

//////////////////////////////////////////////////////////////////////////////

StdLinSetup::StdLinSetup(const std::string& name) :
  CellCenterFVMCom(name),
  socket_linearizedGhostStates("linearizedGhostStates"),
  socket_pastGhostStates("pastGhostStates"),
  socket_gstates("gstates")
{
}

//////////////////////////////////////////////////////////////////////////////

StdLinSetup::~StdLinSetup()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
StdLinSetup::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result;

  result.push_back(&socket_linearizedGhostStates);
  result.push_back(&socket_pastGhostStates);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> >
StdLinSetup::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_gstates);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void StdLinSetup::execute()
{
  CFAUTOTRACE;

  // Create the storage for the linearized ghost states
  DataHandle<State*> gstates = socket_gstates.getDataHandle();

  const CFuint nbGhostStates = gstates.size();

  DataHandle<State*> pastGhostStates = socket_pastGhostStates.getDataHandle();
  pastGhostStates.resize(nbGhostStates);

  DataHandle<State*> linearizedGhostStates = socket_linearizedGhostStates.getDataHandle();
  linearizedGhostStates.resize(nbGhostStates);

  for (CFuint i = 0; i < nbGhostStates; ++i) {
    pastGhostStates[i] = new State();
    linearizedGhostStates[i] = new State();
  }

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
