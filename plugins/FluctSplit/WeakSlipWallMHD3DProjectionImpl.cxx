#include "FluctSplit/FluctSplit.hh"
#include "WeakSlipWallMHD3DProjectionImpl.hh"
#include "InwardNormalsData.hh"
#include "Framework/CFL.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/BlockAccumulator.hh"
#include "Framework/LSSMatrix.hh"
#include "MHD/MHD3DProjectionVarSet.hh"
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

MethodCommandProvider<WeakSlipWallMHD3DProjectionImpl, FluctuationSplitData, FluctSplitModule> weakSlipWallMHD3DProjectionImplProvider("WeakSlipWallMHD3DProjectionImpl");

//////////////////////////////////////////////////////////////////////////////

void WeakSlipWallMHD3DProjectionImpl::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("alpha","distribution coefficient");
   options.addConfigOption< CFreal >("refPhi","Reference phi value imposed at the wall.");
}

//////////////////////////////////////////////////////////////////////////////

WeakSlipWallMHD3DProjectionImpl::WeakSlipWallMHD3DProjectionImpl(const std::string& name) :
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

  _refPhi = 0.;
   setParameter("refPhi",&_refPhi);
}

//////////////////////////////////////////////////////////////////////////////

WeakSlipWallMHD3DProjectionImpl::~WeakSlipWallMHD3DProjectionImpl()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
WeakSlipWallMHD3DProjectionImpl::needsSockets()
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

void WeakSlipWallMHD3DProjectionImpl::setup()
{

  _im0.resize(PhysicalModelStack::getActive()->getNbEq());
  _in0.resize(PhysicalModelStack::getActive()->getNbEq());
  _im1.resize(PhysicalModelStack::getActive()->getNbEq());
  _in1.resize(PhysicalModelStack::getActive()->getNbEq());
  _tJacob.resize(PhysicalModelStack::getActive()->getNbEq(), PhysicalModelStack::getActive()->getNbEq());
  _tJacob = 0.;

  _varSet = getMethodData().getUpdateVar().d_castTo<MHD3DProjectionVarSet>();

  // set up the physical data
  _varSet->getModel()->resizePhysicalData(_physicalData);

  // get the data handle for the states
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  _flagState.resize(states.size());
  _flagState = false;
 }

//////////////////////////////////////////////////////////////////////////////

void WeakSlipWallMHD3DProjectionImpl::executeOnTrs()
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

void WeakSlipWallMHD3DProjectionImpl::computeNormalFluxAndJacob
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
   const CFreal u = _physicalData[MHDProjectionTerm::VX];
   const CFreal v = _physicalData[MHDProjectionTerm::VY];
   const CFreal w = _physicalData[MHDProjectionTerm::VZ];
   const CFreal Bx = _physicalData[MHDProjectionTerm::BX];
   const CFreal By = _physicalData[MHDProjectionTerm::BY];
   const CFreal Bz = _physicalData[MHDProjectionTerm::BZ];

   const CFreal un = u*nx + v*ny + w*nz;
   const CFreal Bn = Bx*nx + By*ny + Bz*nz;
   const CFreal gamma = _varSet->getModel()->getGamma();
   const CFreal beta = _varSet->getModel()->getRefSpeed();
   const CFreal V2 = u*u + v*v + w*w;
   const CFreal B2 = Bx*Bx + By*By + Bz*Bz;
   const CFreal VB = u*Bx + v*By + w*Bz;
   const CFreal rho = _physicalData[MHDProjectionTerm::RHO];
   const CFreal a = _physicalData[MHDProjectionTerm::A];
   const CFreal a2 = a*a;
//   const CFreal _refPhi = _physicalData[MHDProjectionTerm::PHI];

   flux[0] = rho*un;
   flux[1] = un*rho*u - Bn*Bx;
   flux[2] = un*rho*v - Bn*By;
   flux[3] = un*rho*w - Bn*Bz;
   flux[4] = un*Bx - Bn*u + nx*_refPhi;
   flux[5] = un*By - Bn*v + ny*_refPhi;
   flux[6] = un*Bz - Bn*w + nz*_refPhi;
   flux[7] = un*(gamma/(gamma - 1.)*_physicalData[MHDProjectionTerm::P] + 0.5*rho*V2 + B2) - Bn*VB;
   flux[8] = Bn*beta*beta;

     // the analytical jacobian of the normal fluxes
     fluxJacob(0,1) = nx;
     fluxJacob(0,2) = ny;
     fluxJacob(0,3) = nz;

     fluxJacob(1,0) = -u*un;
     fluxJacob(1,1) = un + u*nx;
     fluxJacob(1,2) = u*ny;
     fluxJacob(1,3) = u*nz;
     fluxJacob(1,4) = -Bn - Bx*nx;
     fluxJacob(1,5) = -Bx*ny;
     fluxJacob(1,6) = -Bx*nz;

     fluxJacob(2,0) = -v*un;
     fluxJacob(2,1) = v*nx;
     fluxJacob(2,2) = un + v*ny;
     fluxJacob(2,3) = v*nz;
     fluxJacob(2,4) = -By*nx;
     fluxJacob(2,5) = -Bn - By*ny;
     fluxJacob(2,6) = -By*nz;

     fluxJacob(3,0) = -w*un;
     fluxJacob(3,1) = w*nx;
     fluxJacob(3,2) = w*ny;
     fluxJacob(3,3) = un + w*nz;
     fluxJacob(3,4) = -Bz*nx;
     fluxJacob(3,5) = -Bz*ny;
     fluxJacob(3,6) = -Bn - Bz*nz;

     fluxJacob(4,0) = -un*Bx/rho + u*Bn/rho;
     fluxJacob(4,1) = -By/rho*ny - Bz/rho*nz;
     fluxJacob(4,2) = Bx/rho*ny;
     fluxJacob(4,3) = Bx/rho*nz;
     fluxJacob(4,4) = v*ny + w*nz;
     fluxJacob(4,5) = -u*ny;
     fluxJacob(4,6) = -u*nz;
//     fluxJacob(4,8) = nx;

     fluxJacob(5,0) = -un*By/rho + v*Bn/rho;
     fluxJacob(5,1) = By/rho*nx;
     fluxJacob(5,2) = -Bx/rho*nx - Bz/rho*nz;
     fluxJacob(5,3) = By/rho*nz;
     fluxJacob(5,4) = -v*nx;
     fluxJacob(5,5) = u*nx + w*nz;
     fluxJacob(5,6) = -v*nz;
//     fluxJacob(5,8) = ny;

     fluxJacob(6,0) = -un*Bz/rho + w*Bn/rho;
     fluxJacob(6,1) = Bz/rho*nx;
     fluxJacob(6,2) = Bz/rho*ny;
     fluxJacob(6,3) = -Bx/rho*nx - By/rho*ny;
     fluxJacob(6,4) = -w*nx;
     fluxJacob(6,5) = -w*ny;
     fluxJacob(6,6) = u*nx + v*ny;
//     fluxJacob(6,8) = nz;

     fluxJacob(7,0) = -a2*un/(gamma - 1.) - (2. - gamma)/2.*V2*un - un*B2/rho + Bn/rho*VB;
     fluxJacob(7,1) = (a2/(gamma - 1.) + V2/2. + B2/rho)*nx - (gamma - 1.)*u*un - Bx/rho*Bn;
     fluxJacob(7,2) = (a2/(gamma - 1.) + V2/2. + B2/rho)*ny - (gamma - 1.)*v*un - By/rho*Bn;
     fluxJacob(7,3) = (a2/(gamma - 1.) + V2/2. + B2/rho)*nz - (gamma - 1.)*w*un - Bz/rho*Bn;
     fluxJacob(7,4) = (2. - gamma)*Bx*un - VB*nx - u*Bn;
     fluxJacob(7,5) = (2. - gamma)*By*un - VB*ny - v*Bn;
     fluxJacob(7,6) = (2. - gamma)*Bz*un - VB*nz - w*Bn;
     fluxJacob(7,7) = un*gamma;

     fluxJacob(8,4) = beta*beta*nx;
     fluxJacob(8,5) = beta*beta*ny;
     fluxJacob(8,6) = beta*beta*nz;
}

//////////////////////////////////////////////////////////////////////////////

void WeakSlipWallMHD3DProjectionImpl::setFaceNormal
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
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
