#include "FluctSplit/FluctSplit.hh"
#include "WeakSlipWallMHD3DProjection.hh"
#include "InwardNormalsData.hh"
#include "Framework/CFL.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MethodCommandProvider.hh"
#include "MHD/MHD3DProjectionVarSet.hh"
#include "Framework/MeshData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::MHD;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<WeakSlipWallMHD3DProjection, FluctuationSplitData, FluctSplitModule> weakSlipWallMHD3DProjectionProvider("WeakSlipWallMHD3DProjection");

//////////////////////////////////////////////////////////////////////////////

void WeakSlipWallMHD3DProjection::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("alpha","distribution coefficient");
   options.addConfigOption< CFreal >("refPhi","Reference phi value imposed at the wall.");
}

//////////////////////////////////////////////////////////////////////////////

WeakSlipWallMHD3DProjection::WeakSlipWallMHD3DProjection(const std::string& name) :
  FluctuationSplitCom(name),
  socket_rhs("rhs"),
  socket_states("states"),
  socket_isUpdated("isUpdated"),
  socket_normals("normals"),
  socket_faceNeighCell("faceNeighCell"),
  _flagState(),
  _varSet(CFNULL),
  _physicalData()
{
   addConfigOptionsTo(this);
  _alpha = 0.75;
   setParameter("alpha",&_alpha);

  _refPhi = 0.;
   setParameter("refPhi",&_refPhi);
}

//////////////////////////////////////////////////////////////////////////////

WeakSlipWallMHD3DProjection::~WeakSlipWallMHD3DProjection()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
WeakSlipWallMHD3DProjection::needsSockets()
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

void WeakSlipWallMHD3DProjection::setup()
{
  _varSet = getMethodData().getUpdateVar().d_castTo<MHD3DProjectionVarSet>();

  // set up the physical data
  _varSet->getModel()->resizePhysicalData(_physicalData);

  // get the data handle for the states
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  _flagState.resize(states.size());
  _flagState = false;
 }

//////////////////////////////////////////////////////////////////////////////

void WeakSlipWallMHD3DProjection::executeOnTrs()
{

  DataHandle<bool> isUpdated = socket_isUpdated.getDataHandle();
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  RealVector flux0(0.0, nbEqs);
  RealVector flux1(0.0, nbEqs);
  RealVector flux2(0.0, nbEqs);
  RealVector normal(0.0, dim);
  const CFreal oEminusAlpha = 1. - _alpha;
  const CFreal third = 1./3.;

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
    // compute the normal fluxes corrections for both the states
    // of this cell
    computeNormalFlux(*state0, normal, flux0);
    computeNormalFlux(*state1, normal, flux1);
    computeNormalFlux(*state2, normal, flux2);

    const CFuint stateID0 = state0->getLocalID();
    const CFuint stateID1 = state1->getLocalID();
    const CFuint stateID2 = state2->getLocalID();
    // distribute contributions to the nodes
    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {

      if (!isUpdated[stateID0]) {
        rhs(stateID0, iEq, nbEqs) -=
          third*(_alpha*flux0[iEq] + oEminusAlpha*0.5*
                   (flux1[iEq] + flux2[iEq]));

        _flagState[stateID0] = true;
      }

      if (!isUpdated[stateID1]) {
        rhs(stateID1, iEq, nbEqs) -=
          third*(_alpha*flux1[iEq] + oEminusAlpha*0.5*
          (flux0[iEq] + flux2[iEq]));

        _flagState[stateID1] = true;
      }

      if (!isUpdated[stateID2]) {
        rhs(stateID2, iEq, nbEqs) -=
          third*(_alpha*flux2[iEq] + oEminusAlpha*0.5*
          (flux0[iEq] + flux1[iEq]));

        _flagState[stateID2] = true;
      }
    }

    // release the face
    geoBuilder->releaseGE();
  }

  // this could be done in a more efficient way ...
  const CFuint nbStates = _flagState.size();
  for (CFuint i = 0; i < nbStates; ++i) {
    if (_flagState[i] == true) {
      isUpdated[i] = true;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void WeakSlipWallMHD3DProjection::computeNormalFlux(const State& state,
                                            const RealVector& normal,
                                            RealVector& flux)
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
//   const CFreal _refPhi = _physicalData[MHDProjectionTerm::PHI];

//cout << "phi " << _refPhi << endl;
   flux[0] = rho*un;
   flux[1] = un*rho*u - Bn*Bx;
   flux[2] = un*rho*v - Bn*By;
   flux[3] = un*rho*w - Bn*Bz;
   flux[4] = un*Bx - Bn*u + nx*_refPhi;
   flux[5] = un*By - Bn*v + ny*_refPhi;
   flux[6] = un*Bz - Bn*w + nz*_refPhi;
   flux[7] = un*(gamma/(gamma - 1.)*_physicalData[MHDProjectionTerm::P] + 0.5*rho*V2 + B2) - Bn*VB;
   flux[8] = Bn*beta*beta;
}

//////////////////////////////////////////////////////////////////////////////

void WeakSlipWallMHD3DProjection::setFaceNormal(const CFuint faceID,
                                        RealVector& normal)
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
