#include "FluctSplit/FluctSplitNavierStokes.hh"
#include "StrongMirrorEuler2DCons.hh"
#include "CreateBoundaryNodalNormals.hh"
#include "InwardNormalsData.hh"
#include "Framework/MeshData.hh"
#include "Framework/MethodCommandProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StrongMirrorEuler2DCons,
                      FluctuationSplitData,
                      FluctSplitNavierStokesModule>
aStrongMirrorEuler2DConsProvider("StrongMirrorEuler2DCons");

//////////////////////////////////////////////////////////////////////////////

StrongMirrorEuler2DCons::StrongMirrorEuler2DCons(const std::string& name) :
  FluctuationSplitCom(name),
  socket_rhs("rhs"),
  socket_states("states"),
  socket_isUpdated("isUpdated"),
  socket_normals("normals"),
  socket_faceNeighCell("faceNeighCell"),
  m_bcNormals()
{
}

//////////////////////////////////////////////////////////////////////////////

StrongMirrorEuler2DCons::~StrongMirrorEuler2DCons()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
StrongMirrorEuler2DCons::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_rhs);
  result.push_back(&socket_states);
  result.push_back(&socket_isUpdated);
  result.push_back(&socket_normals);
  result.push_back(&socket_faceNeighCell);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void StrongMirrorEuler2DCons::setup()
{

  FluctuationSplitCom::setup();

  m_bcNormals.resize(getTrsList().size());

  CreateBoundaryNodalNormals obj(getMethodData().getStdTrsGeoBuilder());
  obj.setDataSockets(socket_normals,socket_faceNeighCell);
  obj.create(getTrsList(), m_bcNormals);

}

//////////////////////////////////////////////////////////////////////////////

void StrongMirrorEuler2DCons::executeOnTrs()
{
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<bool> isUpdated = socket_isUpdated.getDataHandle();
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();

  vector<RealVector> *const bcNormalsInTrs = &m_bcNormals[getCurrentTrsID()];

  Common::SafePtr< vector<CFuint> > const statesIdx = getCurrentTRS()->getStatesInTrs();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

  bool not_init_phase = !getMethodData().isInitializationPhase();

  for (CFuint iState = 0; iState < statesIdx->size(); ++iState)
  {
    State *const state = states[(*statesIdx)[iState]];
    const CFuint stateID = state->getLocalID();

    const RealVector* bcNormal = &(*bcNormalsInTrs)[stateID];
    const CFreal nx = (*bcNormal)[0];
    const CFreal ny = (*bcNormal)[1];

    if(not_init_phase)
    {
      if (!isUpdated[stateID])
      {
        const CFreal phyNormal =
          rhs(stateID, 1, nbEqs)*nx +
          rhs(stateID, 2, nbEqs)*ny;

        rhs(stateID, 1, nbEqs) -= phyNormal*nx;
        rhs(stateID, 2, nbEqs) -= phyNormal*ny;
      }
      isUpdated[stateID] = true; // flagging is important!!!!!
    }
    else
    {
      const CFreal vNormal = (*state)[1]*nx + (*state)[2]*ny;

      (*state)[1] -= vNormal*nx;
      (*state)[2] -= vNormal*ny;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
