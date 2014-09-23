#include "FluctSplit/FluctSplitMHD.hh"
#include "StrongSlipWallMHD3DProjection.hh"
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

MethodCommandProvider<StrongSlipWallMHD3DProjection, FluctuationSplitData, FluctSplitMHDModule> strongSlipWallMHD3DProjectionProvider("StrongSlipWallMHD3DProjection");

//////////////////////////////////////////////////////////////////////////////

StrongSlipWallMHD3DProjection::StrongSlipWallMHD3DProjection(const std::string& name) :
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

StrongSlipWallMHD3DProjection::~StrongSlipWallMHD3DProjection()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
StrongSlipWallMHD3DProjection::needsSockets()
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

void StrongSlipWallMHD3DProjection::setup()
{

  FluctuationSplitCom::setup();

  _bcNormals.resize(getTrsList().size());

  CreateBoundaryNodalNormals obj(getMethodData().getStdTrsGeoBuilder());
  obj.setDataSockets(socket_normals,socket_faceNeighCell);
  obj.create(getTrsList(), _bcNormals);

}

//////////////////////////////////////////////////////////////////////////////

void StrongSlipWallMHD3DProjection::executeOnTrs()
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
          const CFreal nz = (*bcNormal)[2];

          CFreal phyNormal = rhs(stateID, 1, nbEqs)*nx +
            rhs(stateID, 2, nbEqs)*ny + rhs(stateID, 3, nbEqs)*nz;

          rhs(stateID, 1, nbEqs) -= phyNormal*nx;
          rhs(stateID, 2, nbEqs) -= phyNormal*ny;
          rhs(stateID, 3, nbEqs) -= phyNormal*nz;

          phyNormal = rhs(stateID, 4, nbEqs)*nx +
            rhs(stateID, 5, nbEqs)*ny + rhs(stateID, 6, nbEqs)*nz ;

          rhs(stateID, 4, nbEqs) -= phyNormal*nx;
          rhs(stateID, 5, nbEqs) -= phyNormal*ny;
          rhs(stateID, 6, nbEqs) -= phyNormal*nz;
        }
        isUpdated[stateID] = true; // flagging is important!!!!!
      }
    }
    else
    {
      State *const state = states[(*statesIdx)[iState]];

      CFreal nx = (*bcNormal)[0];
      CFreal ny = (*bcNormal)[1];
      CFreal nz = (*bcNormal)[2];
      CFreal vNormal = (*state)[1]*nx + (*state)[2]*ny + (*state)[3]*nz;

      (*state)[1] -= vNormal*nx;
      (*state)[2] -= vNormal*ny;
      (*state)[3] -= vNormal*nz;

      vNormal = (*state)[4]*nx + (*state)[5]*ny + (*state)[6]*nz;

      (*state)[4] -= vNormal*nx;
      (*state)[5] -= vNormal*ny;
      (*state)[6] -= vNormal*nz;

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

void StrongSlipWallMHD3DProjection::unsetup()
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
