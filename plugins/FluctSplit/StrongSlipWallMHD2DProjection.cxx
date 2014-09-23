#include "FluctSplit/FluctSplitMHD.hh"
#include "StrongSlipWallMHD2DProjection.hh"
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

MethodCommandProvider<StrongSlipWallMHD2DProjection, FluctuationSplitData, FluctSplitMHDModule> strongSlipWallMHD2DProjectionProvider("StrongSlipWallMHD2DProjection");

//////////////////////////////////////////////////////////////////////////////

StrongSlipWallMHD2DProjection::StrongSlipWallMHD2DProjection(const std::string& name) :
  FluctuationSplitCom(name),
  socket_rhs("rhs"),
  socket_updateCoeff("updateCoeff"),
  socket_states("states"),
  socket_isUpdated("isUpdated"),
  socket_normals("normals"),
  socket_faceNeighCell("faceNeighCell"),
  _bcNormals(),
  _useForInitialization(false)
{
}

//////////////////////////////////////////////////////////////////////////////

StrongSlipWallMHD2DProjection::~StrongSlipWallMHD2DProjection()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
StrongSlipWallMHD2DProjection::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_rhs);
  result.push_back(&socket_updateCoeff);
  result.push_back(&socket_states);
  result.push_back(&socket_isUpdated);
  result.push_back(&socket_normals);
  result.push_back(&socket_faceNeighCell);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void StrongSlipWallMHD2DProjection::setup()
{
  FluctuationSplitCom::setup();

  _bcNormals.resize(getTrsList().size());

  CreateBoundaryNodalNormals obj(getMethodData().getStdTrsGeoBuilder());
  obj.setDataSockets(socket_normals,socket_faceNeighCell);
  obj.create(getTrsList(), _bcNormals);

}

//////////////////////////////////////////////////////////////////////////////

void StrongSlipWallMHD2DProjection::executeOnTrs()
{
  vector<RealVector>* bcNormalsInTrs = &_bcNormals[getCurrentTrsID()];

  DataHandle< CFreal> rhs = socket_rhs.getDataHandle();
  DataHandle< CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle< bool> isUpdated = socket_isUpdated.getDataHandle();

  Common::SafePtr< vector<CFuint> > const statesIdx = getCurrentTRS()->getStatesInTrs();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  _useForInitialization = getMethodData().isInitializationPhase();

  for (CFuint iState = 0; iState < statesIdx->size(); ++iState) {
    const CFuint stateID = (*statesIdx)[iState];
    const RealVector* bcNormal = &(*bcNormalsInTrs)[stateID];

    if(!_useForInitialization){
      if (!isUpdated[stateID]) {
        const CFreal updateCoeffValue = updateCoeff[stateID];
        if (std::abs(updateCoeffValue) > MathTools::MathConsts::CFrealEps()) {
          const CFreal nx = (*bcNormal)[0];
          const CFreal ny = (*bcNormal)[1];

          CFreal phyNormal = rhs(stateID, 1, nbEqs)*nx +
            rhs(stateID, 2, nbEqs)*ny;

          rhs(stateID, 1, nbEqs) -= phyNormal*nx;
          rhs(stateID, 2, nbEqs) -= phyNormal*ny;

          phyNormal = rhs(stateID, 4, nbEqs)*nx +
            rhs(stateID, 5, nbEqs)*ny;

          rhs(stateID, 4, nbEqs) -= phyNormal*nx;
          rhs(stateID, 5, nbEqs) -= phyNormal*ny;
        }
        isUpdated[stateID] = true; // flagging is important!!!!!
      }
    }
    else
    {
      State *const state = states[(*statesIdx)[iState]];

      CFreal nx = (*bcNormal)[0];
      CFreal ny = (*bcNormal)[1];
      CFreal vNormal = (*state)[1]*nx + (*state)[2]*ny;

      (*state)[1] -= vNormal*nx;
      (*state)[2] -= vNormal*ny;

      vNormal = (*state)[4]*nx + (*state)[5]*ny;

      (*state)[4] -= vNormal*nx;
      (*state)[5] -= vNormal*ny;

      //  (*state)[7] = 1./(5./3.-1.) + 0.5*(((*state)[1]*(*state)[1] +
      //                                  (*state)[2]*(*state)[2] +
      //                                  (*state)[3]*(*state)[3])/(*state)[0] +
      //                                 (*state)[4]*(*state)[4] +
      //                                 (*state)[5]*(*state)[5] +
      //                                 (*state)[6]*(*state)[6]);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void StrongSlipWallMHD2DProjection::unsetup()
{
  DataHandle< CFreal> rhs = socket_rhs.getDataHandle();
  DataHandle< CFreal> updateCoeff = socket_updateCoeff.getDataHandle();

  rhs = DataHandle<CFreal>(CFNULL);
  updateCoeff = DataHandle<CFreal>(CFNULL);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
