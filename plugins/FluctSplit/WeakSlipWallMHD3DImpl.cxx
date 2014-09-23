#include "FluctSplit/FluctSplit.hh"
#include "WeakSlipWallMHD3DImpl.hh"
#include "InwardNormalsData.hh"
#include "Framework/CFL.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/BlockAccumulator.hh"
#include "Framework/LSSMatrix.hh"
#include "MHD/MHD3DVarSet.hh"
#include "Framework/MeshData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::MHD;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<WeakSlipWallMHD3DImpl, FluctuationSplitData, FluctSplitModule> weakSlipWallMHD3DImplProvider("WeakSlipWallMHD3DImpl");

//////////////////////////////////////////////////////////////////////////////

void WeakSlipWallMHD3DImpl::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("alpha","distribution coefficient");
}

//////////////////////////////////////////////////////////////////////////////

WeakSlipWallMHD3DImpl::WeakSlipWallMHD3DImpl(const std::string& name) :
  FluctuationSplitCom(name),
  socket_rhs("rhs"),
  socket_states("states"),
  socket_isUpdated("isUpdated"),
  socket_normals("normals"),
  socket_faceNeighCell("faceNeighCell"),
  _varSet(CFNULL),
  _physicalData(),
  _im0(),
  _in0(),
  _im1(),
  _in1(),
  _flagState(),
  _tJacob()
{
   addConfigOptionsTo(this);
  _alpha = 0.75;
   setParameter("alpha",&_alpha);
}

//////////////////////////////////////////////////////////////////////////////

WeakSlipWallMHD3DImpl::~WeakSlipWallMHD3DImpl()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
WeakSlipWallMHD3DImpl::needsSockets()
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

void WeakSlipWallMHD3DImpl::setup()
{

  _im0.resize(PhysicalModelStack::getActive()->getNbEq());
  _in0.resize(PhysicalModelStack::getActive()->getNbEq());
  _im1.resize(PhysicalModelStack::getActive()->getNbEq());
  _in1.resize(PhysicalModelStack::getActive()->getNbEq());
  _tJacob.resize(PhysicalModelStack::getActive()->getNbEq(), PhysicalModelStack::getActive()->getNbEq());
  _tJacob = 0.;

  _varSet = getMethodData().getUpdateVar().d_castTo<MHD3DVarSet>();

  // set up the physical data
  _varSet->getModel()->resizePhysicalData(_physicalData);

  // get the data handle for the states
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  _flagState.resize(states.size());
  _flagState = false;
 }

//////////////////////////////////////////////////////////////////////////////

void WeakSlipWallMHD3DImpl::executeOnTrs()
{

  DataHandle<bool> isUpdated = socket_isUpdated.getDataHandle();
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFreal third = 1./3.;
  const CFreal halfOEminAlphaThird = 0.5*(1. - _alpha)*third;
  const CFreal alphaThird = _alpha*third;

  RealVector normal(dim);
  RealVector flux0(0.0, nbEqs);
  RealVector flux1(0.0, nbEqs);
  RealVector flux2(0.0, nbEqs);
  RealMatrix fluxJacob0(nbEqs, nbEqs, 0.0);
  RealMatrix fluxJacob1(nbEqs, nbEqs, 0.0);
  RealMatrix fluxJacob2(nbEqs, nbEqs, 0.0);

  SafePtr<LinearSystemSolver> lss =
    getMethodData().getLinearSystemSolver()[0];

  SafePtr<LSSMatrix> jacobMatrix = lss->getMatrix();

  // block accumulator 3*3
  auto_ptr<BlockAccumulator> acc(lss->createBlockAccumulator(3, 3, nbEqs));

  Common::SafePtr<GeometricEntityPool<StdTrsGeoBuilder> >
    geoBuilder = getMethodData().getStdTrsGeoBuilder();

  StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.trs = getCurrentTRS();

  const CFuint nbFaces = getCurrentTRS()->getLocalNbGeoEnts();
  for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
    geoData.idx = iFace;
    GeometricEntity *const currFace = geoBuilder->buildGE();

    // set the face normal
    setFaceNormal(currFace->getID(), normal);

    State *const state0 = currFace->getState(0);
    State *const state1 = currFace->getState(1);
    State *const state2 = currFace->getState(2);
    const CFuint stateID0 = state0->getLocalID();
    const CFuint stateID1 = state1->getLocalID();
    const CFuint stateID2 = state2->getLocalID();

    // set the row - column
    acc->setRowColIndex(0, stateID0);
    acc->setRowColIndex(1, stateID1);
    acc->setRowColIndex(2, stateID2);

    // compute the normal fluxes corrections for both the states
    // of this cell
    computeNormalFluxAndJacob(*state0, normal, flux0, fluxJacob0);
    computeNormalFluxAndJacob(*state1, normal, flux1, fluxJacob1);
    computeNormalFluxAndJacob(*state2, normal, flux2, fluxJacob2);

    // distribute contributions to the nodes
      for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
        rhs(stateID0, iEq, nbEqs) -=
          alphaThird*flux0[iEq] + halfOEminAlphaThird*
          (flux1[iEq] + flux2[iEq]);

        rhs(stateID1, iEq, nbEqs) -=
            alphaThird*flux1[iEq] + halfOEminAlphaThird*
          (flux0[iEq] + flux2[iEq]);

        rhs(stateID2, iEq, nbEqs) -=
          alphaThird*flux2[iEq] + halfOEminAlphaThird*
          (flux0[iEq] + flux1[iEq]);
      }

    // compute the jacobian contributions for the involved states
    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      for (CFuint jEq = 0; jEq < nbEqs; ++jEq) {

        if (state0->isParUpdatable()) {
          acc->setValue(0, 0, iEq, jEq,
                        alphaThird*fluxJacob0(iEq, jEq));

          acc->setValue(0, 1, iEq, jEq,
                        halfOEminAlphaThird*fluxJacob1(iEq, jEq));

          acc->setValue(0, 2, iEq, jEq,
                        halfOEminAlphaThird*fluxJacob2(iEq, jEq));
        }

        if (state1->isParUpdatable()) {
          acc->setValue(1, 0, iEq, jEq,
                        halfOEminAlphaThird*fluxJacob0(iEq, jEq));

          acc->setValue(1, 1, iEq, jEq,
                        alphaThird*fluxJacob1(iEq, jEq));

          acc->setValue(1, 2, iEq, jEq,
                        halfOEminAlphaThird*fluxJacob2(iEq, jEq));
        }

        if (state2->isParUpdatable()) {
          acc->setValue(2, 0, iEq, jEq,
                        halfOEminAlphaThird*fluxJacob0(iEq, jEq));

          acc->setValue(2, 1, iEq, jEq,
                        halfOEminAlphaThird*fluxJacob1(iEq, jEq));

          acc->setValue(2, 2, iEq, jEq,
                        alphaThird*fluxJacob2(iEq, jEq));
        }
      }
    }

    // add the contributions to the jacobian
    jacobMatrix->addValues(*acc);

    acc->reset();

    // release the face
    geoBuilder->releaseGE();
  }
}

//////////////////////////////////////////////////////////////////////////////

void WeakSlipWallMHD3DImpl::computeNormalFluxAndJacob
(const State& state,
 const RealVector& normal,
 RealVector& flux,
 RealMatrix& fluxJacob)
{
  State& ss = *(const_cast<State*>(&state));
  _varSet->computePhysicalData(ss, _physicalData);

   const CFreal nx = normal[XX];
   const CFreal ny = normal[YY];
   const CFreal nz = normal[ZZ];
   const CFreal u = _physicalData[MHDTerm::VX];
   const CFreal v = _physicalData[MHDTerm::VY];
   const CFreal w = _physicalData[MHDTerm::VZ];
   const CFreal Bx = _physicalData[MHDTerm::BX];
   const CFreal By = _physicalData[MHDTerm::BY];
   const CFreal Bz = _physicalData[MHDTerm::BZ];
   const CFreal a = _physicalData[MHDTerm::A];

   const CFreal un = u*nx + v*ny + w*nz;
   // unused // const CFreal Bn = Bx*nx + By*ny + Bz*nz;
   const CFreal gamma = _varSet->getModel()->getGamma();
   const CFreal V2 = u*u + v*v + w*w;
   const CFreal B2 = Bx*Bx + By*By + Bz*Bz;
   // unused // const CFreal VB = u*Bx + v*By + w*Bz;
   const CFreal rho = _physicalData[MHDTerm::RHO];
   const CFreal a2 = a*a;

   flux[0] = rho*un;
   flux[1] = un*rho*u;
   flux[2] = un*rho*v;
   flux[3] = un*rho*w;
   flux[4] = un*Bx;
   flux[5] = un*By;
   flux[6] = un*Bz;
   flux[7] = un*(gamma/(gamma - 1.)*_physicalData[MHDTerm::P] + 0.5*rho*V2 + B2);

   if (!getMethodData().isResidualTransformationNeeded()) {
     // the analytical jacobian of the normal fluxes
     fluxJacob(0,1) = nx;
     fluxJacob(0,2) = ny;
     fluxJacob(0,3) = nz;

     fluxJacob(1,0) = -u*un;
     fluxJacob(1,1) = un + u*nx;
     fluxJacob(1,2) = u*ny;
     fluxJacob(1,3) = u*nz;

     fluxJacob(2,0) = -v*un;
     fluxJacob(2,1) = v*nx;
     fluxJacob(2,2) = un + v*ny;
     fluxJacob(2,3) = v*nz;

     fluxJacob(3,0) = -w*un;
     fluxJacob(3,1) = w*nx;
     fluxJacob(3,2) = w*ny;
     fluxJacob(3,3) = un + w*nz;

     fluxJacob(4,0) = -un*Bx/rho;
     fluxJacob(4,1) = Bx/rho*nx;
     fluxJacob(4,2) = Bx/rho*ny;
     fluxJacob(4,3) = Bx/rho*nz;
     fluxJacob(4,4) = un;

     fluxJacob(5,0) = -un*By/rho;
     fluxJacob(5,1) = By/rho*nx;
     fluxJacob(5,2) = By/rho*ny;
     fluxJacob(5,3) = By/rho*nz;
     fluxJacob(5,5) = un;

     fluxJacob(6,0) = -un*Bz/rho;
     fluxJacob(6,1) = Bz/rho*nx;
     fluxJacob(6,2) = Bz/rho*ny;
     fluxJacob(6,3) = Bz/rho*nz;
     fluxJacob(6,6) = un;

     fluxJacob(7,0) = -a2*un/(gamma - 1.) + 0.5*(gamma - 2.)*V2*un - un*B2/rho;
     fluxJacob(7,1) = (a2/(gamma - 1.) + V2/2. + B2/rho)*nx - (gamma - 1.)*u*un;
     fluxJacob(7,2) = (a2/(gamma - 1.) + V2/2. + B2/rho)*ny - (gamma - 1.)*v*un;
     fluxJacob(7,3) = (a2/(gamma - 1.) + V2/2. + B2/rho)*nz - (gamma - 1.)*w*un;
     fluxJacob(7,4) = (2. - gamma)*Bx*un;
     fluxJacob(7,5) = (2. - gamma)*By*un;
     fluxJacob(7,6) = (2. - gamma)*Bz*un;
     fluxJacob(7,7) = un*gamma;
   }

   else {
     // the analytical jacobian of the normal fluxes
     _tJacob(0,1) = nx;
     _tJacob(0,2) = ny;
     _tJacob(0,3) = nz;

     _tJacob(1,0) = -u*un;
     _tJacob(1,1) = un + u*nx;
     _tJacob(1,2) = u*ny;
     _tJacob(1,3) = u*nz;

     _tJacob(2,0) = -v*un;
     _tJacob(2,1) = v*nx;
     _tJacob(2,2) = un + v*ny;
     _tJacob(2,3) = v*nz;

     _tJacob(3,0) = -w*un;
     _tJacob(3,1) = w*nx;
     _tJacob(3,2) = w*ny;
     _tJacob(3,3) = un + w*nz;

     _tJacob(4,0) = -un*Bx/rho;
     _tJacob(4,1) = Bx/rho*nx;
     _tJacob(4,2) = Bx/rho*ny;
     _tJacob(4,3) = Bx/rho*nz;
     _tJacob(4,4) = un;

     _tJacob(5,0) = -un*By/rho;
     _tJacob(5,1) = By/rho*nx;
     _tJacob(5,2) = By/rho*ny;
     _tJacob(5,3) = By/rho*nz;
     _tJacob(5,5) = un;

     _tJacob(6,0) = -un*Bz/rho;
     _tJacob(6,1) = Bz/rho*nx;
     _tJacob(6,2) = Bz/rho*ny;
     _tJacob(6,3) = Bz/rho*nz;
     _tJacob(6,6) = un;

     _tJacob(7,0) = -a2*un/(gamma - 1.) + 0.5*(gamma - 2.)*V2*un - un*B2/rho;
     _tJacob(7,1) = (a2/(gamma - 1.) + V2/2. + B2/rho)*nx - (gamma - 1.)*u*un;
     _tJacob(7,2) = (a2/(gamma - 1.) + V2/2. + B2/rho)*ny - (gamma - 1.)*v*un;
     _tJacob(7,3) = (a2/(gamma - 1.) + V2/2. + B2/rho)*nz - (gamma - 1.)*w*un;
     _tJacob(7,4) = (2. - gamma)*Bx*un;
     _tJacob(7,5) = (2. - gamma)*By*un;
     _tJacob(7,6) = (2. - gamma)*Bz*un;
     _tJacob(7,7) = un*gamma;

     // set the transformation from update to solution in update
     SafePtr<VarSetMatrixTransformer> updateToSolInUpdate =
       getMethodData().getUpdateToSolutionInUpdateMatTrans();

     updateToSolInUpdate->setMatrix(state);
     const RealMatrix& tMatrix = *updateToSolInUpdate->getMatrix();

     fluxJacob = _tJacob*tMatrix;
   }
}

//////////////////////////////////////////////////////////////////////////////

void WeakSlipWallMHD3DImpl::setFaceNormal
(const CFuint faceID, RealVector& normal)
{
  // get the data handle for the states
  DataHandle< InwardNormalsData*> normals = socket_normals.getDataHandle();
  DataHandle<Common::Trio<CFuint, CFuint, Common::SafePtr<Framework::TopologicalRegionSet> > >
    faceNeighCell = socket_faceNeighCell.getDataHandle();

  const CFuint cellTrsID = faceNeighCell[faceID].first;
  const CFuint iFaceLocal = faceNeighCell[faceID].second;
  Common::SafePtr<TopologicalRegionSet> cellTrs = faceNeighCell[faceID].third;
  const CFuint cellLocalID = cellTrs->getLocalGeoID(cellTrsID);

  normal[XX] = -normals[cellLocalID]->getFaceNormComp(iFaceLocal, XX);
  normal[YY] = -normals[cellLocalID]->getFaceNormComp(iFaceLocal, YY);
  normal[ZZ] = -normals[cellLocalID]->getFaceNormComp(iFaceLocal, ZZ);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
