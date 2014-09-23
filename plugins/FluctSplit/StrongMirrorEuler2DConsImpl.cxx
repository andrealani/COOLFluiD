#include "FluctSplit/FluctSplitNavierStokes.hh"
#include "StrongMirrorEuler2DConsImpl.hh"
#include "CreateBoundaryNodalNormals.hh"
#include "InwardNormalsData.hh"
#include "Framework/CFL.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/MeshData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StrongMirrorEuler2DConsImpl, FluctuationSplitData, FluctSplitNavierStokesModule> strongMirrorEuler2DConsImplProvider("StrongMirrorEuler2DConsImpl");

//////////////////////////////////////////////////////////////////////////////

StrongMirrorEuler2DConsImpl::StrongMirrorEuler2DConsImpl(const std::string& name) :
  FluctuationSplitCom(name),
  socket_rhs("rhs"),
  socket_states("states"),
  socket_updateCoeff("updateCoeff"),
  socket_isUpdated("isUpdated"),
  socket_bStatesNeighbors("bStatesNeighbors"),
  socket_normals("normals"),
  socket_faceNeighCell("faceNeighCell"),
  _bcNormals(),
  _block(2,4),
  _jacobElem(2,4),
  _im(2),
  _in(4)
{
}

//////////////////////////////////////////////////////////////////////////////

StrongMirrorEuler2DConsImpl::~StrongMirrorEuler2DConsImpl()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
StrongMirrorEuler2DConsImpl::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_rhs);
  result.push_back(&socket_states);
  result.push_back(&socket_updateCoeff);
  result.push_back(&socket_isUpdated);
  result.push_back(&socket_bStatesNeighbors);
  result.push_back(&socket_normals);
  result.push_back(&socket_faceNeighCell);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void StrongMirrorEuler2DConsImpl::setup()
{
  FluctuationSplitCom::setup();

  _bcNormals.resize(getTrsList().size());

  CreateBoundaryNodalNormals obj(getMethodData().getStdTrsGeoBuilder());
  obj.setDataSockets(socket_normals,socket_faceNeighCell);
  obj.create(getTrsList(), _bcNormals);

}

//////////////////////////////////////////////////////////////////////////////

void StrongMirrorEuler2DConsImpl::executeOnTrs()
{
  vector<RealVector> *const bcNormalsInTrs = &_bcNormals[getCurrentTrsID()];

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle< std::valarray<State*> > bStatesNeighbors =
    socket_bStatesNeighbors.getDataHandle();
  DataHandle<bool> isUpdated = socket_isUpdated.getDataHandle();
  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();

  Common::SafePtr< vector<CFuint> > const statesIdx = getCurrentTRS()->getStatesInTrs();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

  SafePtr<LinearSystemSolver> lss =
    getMethodData().getLinearSystemSolver()[0];

  SafePtr<LSSMatrix> jacobMatrix = lss->getMatrix();

  jacobMatrix->finalAssembly();

  const CFreal cfl = getMethodData().getCFL()->getCFLValue();

  const CFuint nbVarsInvolved = 2;
  const CFuint rhouVarID = 1;
  const CFuint rhovVarID = 2;

  for (CFuint iState = 0; iState < statesIdx->size(); ++iState) {
    const CFuint localStateID = (*statesIdx)[iState];
    State *const state = states[localStateID];

    if (state->isParUpdatable()) {

      const RealVector* bcNormal = &(*bcNormalsInTrs)[localStateID];
      if (!isUpdated[localStateID]) {

        const CFreal nx = (*bcNormal)[XX];
        const CFreal ny = (*bcNormal)[YY];
        const CFreal phiNormal = rhs(localStateID, 1, nbEqs)*nx +
          rhs(localStateID, 2, nbEqs)*ny;

        rhs(localStateID, 1, nbEqs) -= phiNormal*nx;
        rhs(localStateID, 2, nbEqs) -= phiNormal*ny;
if (getMethodData().doComputeJacobian()) {
        const CFuint nbEntries = bStatesNeighbors[localStateID].size();
        cf_assert(nbEntries > 0);

        // the values on the block diagonal have to be reduced
        // by the "update" coefficient 1/dt
        const CFreal coeff = updateCoeff[localStateID]/cfl;

        for (CFuint i = 0; i < nbEntries; ++i) {
          const CFuint entryID = bStatesNeighbors[localStateID][i]->getLocalID();
          const CFuint nStart = nbEqs*entryID;

          // set the row ids
          _im[XX] = nStart + rhouVarID;
          _im[YY] = nStart + rhovVarID;
          // set the column ids
          for (CFuint iVar = 0; iVar < nbEqs; ++iVar) {
            _in[0] = nStart + iVar;
          }

          // get the elements of the jacobian matrix
          jacobMatrix->getValues(nbVarsInvolved,
				 &_im[0],
				 nbEqs,
				 &_in[0],
				 &_jacobElem[0]);

	  jacobMatrix->finalAssembly();

          // reset the values in the block to zero
          _block = 0.0;

          if (entryID == state->getLocalID()) {
            // is this needed ?
            _jacobElem(XX,1) = coeff - _jacobElem(XX,1);
            _jacobElem(YY,2) = coeff - _jacobElem(YY,2);
          }

          for (CFuint iVar = 0; iVar < nbEqs; ++iVar) {
            // d(Normal residual)/dU_iVar
            const CFreal dPhiNormalDu =
              _jacobElem(XX, iVar)*nx + _jacobElem(YY, iVar)*ny;

            _block(XX, iVar) = -dPhiNormalDu*nx;
            _block(YY, iVar) = -dPhiNormalDu*ny;
          }

          jacobMatrix->addValues(nbVarsInvolved,
                                 &_im[0],
                                 nbEqs,
                                 &_in[0],
                                 &_block[0]);

	  jacobMatrix->finalAssembly();
        }
 }
        isUpdated[localStateID] = true; // flagging is important!!!!!
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
