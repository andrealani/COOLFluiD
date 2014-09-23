#include "FluctSplit/FluctSplitSpaceTimeNavierStokes.hh"
#include "TwoLayerStrongMirrorEuler2DCons.hh"
#include "CreateBoundaryNodalNormals.hh"
#include "InwardNormalsData.hh"
#include "Framework/MeshData.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/SubSystemStatus.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<TwoLayerStrongMirrorEuler2DCons, FluctuationSplitData, FluctSplitSpaceTimeNavierStokesModule> TwoLayerStrongMirrorEuler2DConsProvider("TwoLayerStrongMirrorEuler2DCons");

//////////////////////////////////////////////////////////////////////////////

TwoLayerStrongMirrorEuler2DCons::TwoLayerStrongMirrorEuler2DCons(const std::string& name) :
  StrongMirrorEuler2DCons(name),
  socket_interRhs("interRhs"),
  socket_interStates("interStates")
{
}

//////////////////////////////////////////////////////////////////////////////

TwoLayerStrongMirrorEuler2DCons::~TwoLayerStrongMirrorEuler2DCons()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
TwoLayerStrongMirrorEuler2DCons::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result = StrongMirrorEuler2DCons::needsSockets();

  result.push_back(&socket_interRhs);
  result.push_back(&socket_interStates);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void TwoLayerStrongMirrorEuler2DCons::executeOnTrs()
{
  vector<RealVector> *const bcNormalsInTrs = &m_bcNormals[getCurrentTrsID()];

  Common::SafePtr<std::vector<CFuint> > statesIdx = getCurrentTRS()->getStatesInTrs();

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<State*> interStates  = socket_interStates.getDataHandle();
  DataHandle<bool> isUpdated = socket_isUpdated.getDataHandle();
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  DataHandle<CFreal> interRhs = socket_interRhs.getDataHandle();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  bool not_init_phase = !getMethodData().isInitializationPhase();

  for (CFuint iState = 0; iState < statesIdx->size(); ++iState) {
    const CFuint stateID = (*statesIdx)[iState];
    State *const state = states[stateID];
    const RealVector* bcNormal = &(*bcNormalsInTrs)[stateID];

    if(not_init_phase)
    {
      if (!isUpdated[stateID])
      {
        const CFreal nx = (*bcNormal)[0];
        const CFreal ny = (*bcNormal)[1];

        CFreal phyNormal = interRhs(stateID, 1, nbEqs)*nx + interRhs(stateID, 2, nbEqs)*ny;

        interRhs(stateID, 1, nbEqs) -= phyNormal*nx;
        interRhs(stateID, 2, nbEqs) -= phyNormal*ny;

        phyNormal = rhs(stateID, 1, nbEqs)*nx + rhs(stateID, 2, nbEqs)*ny;

        rhs(stateID, 1, nbEqs) -= phyNormal*nx;
        rhs(stateID, 2, nbEqs) -= phyNormal*ny;

      }
      isUpdated[stateID] = true; // flagging is important!!!!!
    }
    else
    {
      const CFreal nx = (*bcNormal)[0];
      const CFreal ny = (*bcNormal)[1];
      const CFreal vNormal = (*state)[1]*nx + (*state)[2]*ny;

      (*state)[1] -= vNormal*nx;
      (*state)[2] -= vNormal*ny;

      *interStates[stateID] = *state;
    }

  }

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
