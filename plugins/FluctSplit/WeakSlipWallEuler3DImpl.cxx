#include "FluctSplit/FluctSplitNavierStokes.hh"
#include "WeakSlipWallEuler3DImpl.hh"
#include "InwardNormalsData.hh"
#include "Framework/CFL.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/MeshData.hh"
#include "Framework/NamespaceSwitcher.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<WeakSlipWallEuler3DImpl, FluctuationSplitData, FluctSplitNavierStokesModule> 
weakSlipWallEuler3DImplProvider("WeakSlipWallEuler3DImpl");

MethodCommandProvider<WeakSlipWallEuler3DImpl, FluctuationSplitData, FluctSplitNavierStokesModule> 
weakSlipWallEuler3DConsImplProvider("WeakSlipWallEuler3DConsImpl");

//////////////////////////////////////////////////////////////////////////////

void WeakSlipWallEuler3DImpl::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("alpha","distribution coefficient");
}

//////////////////////////////////////////////////////////////////////////////

WeakSlipWallEuler3DImpl::WeakSlipWallEuler3DImpl(const std::string& name) :
  FluctuationSplitCom(name),
  socket_rhs("rhs"),
  socket_states("states"),
  socket_isUpdated("isUpdated"),
  socket_normals("normals"),
  socket_faceNeighCell("faceNeighCell"),
  _varSet(CFNULL),
  _im0(),
  _in0(),
  _im1(),
  _in1(),
  _accumulators(),
  _fluxes(),
  _fluxJacobs(),
  _physicalData(),
  _tJacob()
{
   addConfigOptionsTo(this);
  _alpha = 1.0;
   setParameter("alpha",&_alpha);
}

//////////////////////////////////////////////////////////////////////////////

WeakSlipWallEuler3DImpl::~WeakSlipWallEuler3DImpl()
 {
  for (CFuint i = 0; i < _accumulators.size(); ++i) {
    deletePtr(_accumulators[i]);
  }
 }

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
WeakSlipWallEuler3DImpl::needsSockets()
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

void WeakSlipWallEuler3DImpl::setup()
{

  FluctuationSplitCom::setup();
  
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

  _im0.resize(nbEqs);
  _in0.resize(nbEqs);
  _im1.resize(nbEqs);
  _in1.resize(nbEqs);

  _accumulators.resize(2);
  _accumulators[0] = getMethodData().getLinearSystemSolver()[0]->
    createBlockAccumulator(3,3,nbEqs);
  // this second one is needed if you deal with hybrid 3D meshes
  _accumulators[1] = getMethodData().getLinearSystemSolver()[0]->
    createBlockAccumulator(4,4,nbEqs);

  // size the fluxes
  _fluxes.resize(MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell()-1);
  for (CFuint i = 0; i < _fluxes.size(); ++i) {
    _fluxes[i].resize(nbEqs);
  }

  // size the flux jacobians
  _fluxJacobs.resize(MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell()-1);
  for (CFuint i = 0; i < _fluxJacobs.size(); ++i) {
    _fluxJacobs[i].resize(nbEqs, nbEqs);
  }
  
  _varSet = getMethodData().getUpdateVar().d_castTo<Euler3DVarSet>();
  
  // set up the physical data
  _varSet->getModel()->resizePhysicalData(_physicalData);
  
  _tJacob.resize(PhysicalModelStack::getActive()->getNbEq(), 
		 PhysicalModelStack::getActive()->getNbEq());
  _tJacob = 0.;
}

//////////////////////////////////////////////////////////////////////////////

void WeakSlipWallEuler3DImpl::executeOnTrs()
{
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFreal third = 1./3;
  const CFreal halfOEminAlphaThird = 0.5*(1. - _alpha)*third;
  const CFreal alphaThird = _alpha*third;

  DataHandle<bool> isUpdated = socket_isUpdated.getDataHandle();
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();

  RealVector normal(dim);

  SafePtr<LinearSystemSolver> lss =
    getMethodData().getLinearSystemSolver()[0];

  SafePtr<LSSMatrix> jacobMatrix = lss->getMatrix();

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

    const vector<State*> *const statesInFace =
      currFace->getStates();
    const CFuint nbStatesInFace = statesInFace->size();

    for (CFuint is = 0; is < nbStatesInFace; ++is) {
      // compute the normal fluxes corrections for both the states
      // of this cell

      /// WARNING the use of std::valarray<RealMatrix> could give a
      /// problem with icc ... CHECK IT
      computeNormalFluxAndJacob(*(*statesInFace)[is],
				normal,
				_fluxes[is],
                                  _fluxJacobs[is]);
    }

    // block accumulator 2*2
    BlockAccumulator *const acc = getAccumulator(nbStatesInFace);

    /// @todo implement case when the face is quadrilateral
    if (nbStatesInFace == 3) {
      const CFuint stateID0 = (*statesInFace)[0]->getLocalID();
      const CFuint stateID1 = (*statesInFace)[1]->getLocalID();
      const CFuint stateID2 = (*statesInFace)[2]->getLocalID();

      // set the row - column
      acc->setRowColIndex(0, stateID0);
      acc->setRowColIndex(1, stateID1);
      acc->setRowColIndex(2, stateID2);

      // distribute contributions to the nodes
      for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
	rhs(stateID0, iEq, nbEqs) -=
	  alphaThird*_fluxes[0][iEq] + halfOEminAlphaThird*
	  (_fluxes[1][iEq] + _fluxes[2][iEq]);

	rhs(stateID1, iEq, nbEqs) -=
            alphaThird*_fluxes[1][iEq] + halfOEminAlphaThird*
	  (_fluxes[0][iEq] + _fluxes[2][iEq]);

	rhs(stateID2, iEq, nbEqs) -=
	  alphaThird*_fluxes[2][iEq] + halfOEminAlphaThird*
	  (_fluxes[0][iEq] + _fluxes[1][iEq]);
      }

      // compute the jacobian contributions for the involved states
      for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
	for (CFuint jEq = 0; jEq < nbEqs; ++jEq) {

	  if ((*statesInFace)[0]->isParUpdatable()) {
	    acc->setValue(0, 0, iEq, jEq,
			  alphaThird*_fluxJacobs[0](iEq, jEq));

	    acc->setValue(0, 1, iEq, jEq,
			  halfOEminAlphaThird*_fluxJacobs[1](iEq, jEq));

	    acc->setValue(0, 2, iEq, jEq,
			  halfOEminAlphaThird*_fluxJacobs[2](iEq, jEq));
	  }


	  if ((*statesInFace)[1]->isParUpdatable()) {
	    acc->setValue(1, 0, iEq, jEq,
			  halfOEminAlphaThird*_fluxJacobs[0](iEq, jEq));

	    acc->setValue(1, 1, iEq, jEq,
			  alphaThird*_fluxJacobs[1](iEq, jEq));

	    acc->setValue(1, 2, iEq, jEq,
			  halfOEminAlphaThird*_fluxJacobs[2](iEq, jEq));
	  }

	  if ((*statesInFace)[2]->isParUpdatable()) {
	    acc->setValue(2, 0, iEq, jEq,
			  halfOEminAlphaThird*_fluxJacobs[0](iEq, jEq));

	    acc->setValue(2, 1, iEq, jEq,
			  halfOEminAlphaThird*_fluxJacobs[1](iEq, jEq));

              acc->setValue(2, 2, iEq, jEq,
                            alphaThird*_fluxJacobs[2](iEq, jEq));
	  }
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

void WeakSlipWallEuler3DImpl::computeNormalFluxAndJacob(const State& state,
							const RealVector& normal,
							RealVector& flux,
							RealMatrix& fluxJacob)
{
   State& ss = *(const_cast<State*>(&state));
   _varSet->computePhysicalData(ss, _physicalData);
   
   const CFreal nx = normal[XX];
   const CFreal ny = normal[YY];
   const CFreal nz = normal[ZZ];
   const CFreal rho = _physicalData[EulerTerm::RHO];
   const CFreal u = _physicalData[EulerTerm::VX];
   const CFreal v = _physicalData[EulerTerm::VY];
   const CFreal w = _physicalData[EulerTerm::VZ];
   const CFreal un = u*nx + v*ny + w*nz;
   const CFreal gamma = _varSet->getModel()->getGamma();
   const CFreal gammaMinus1 = gamma - 1.;
   const CFreal halfGammaMinus1 = 0.5*gammaMinus1;
   const CFreal vel2 = u*u + v*v +w*w;
// unused // const CFreal rhoV2 = rho*vel2;
   const CFreal rhoH = rho*_physicalData[EulerTerm::H];
   
   flux[0] = rho*un;
   flux[1] = un*rho*u;
   flux[2] = un*rho*v;
   flux[3] = un*rho*w;
   flux[4] = un*rhoH;
   
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
     
     fluxJacob(4,0) = un*(halfGammaMinus1*vel2 - rhoH/rho);
     fluxJacob(4,1) = rhoH/rho*nx - un*u*gammaMinus1;
     fluxJacob(4,2) = rhoH/rho*ny - un*v*gammaMinus1;
     fluxJacob(4,3) = rhoH/rho*nz - un*w*gammaMinus1;
     fluxJacob(4,4) = gamma*un;
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
     
     _tJacob(4,0) = un*(halfGammaMinus1*vel2 - rhoH/rho);
     _tJacob(4,1) = rhoH/rho*nx - un*u*gammaMinus1;
     _tJacob(4,2) = rhoH/rho*ny - un*v*gammaMinus1;
     _tJacob(4,3) = rhoH/rho*nz - un*w*gammaMinus1;
     _tJacob(4,4) = gamma*un;
     
     // set the transformation from update to solution in update
     SafePtr<VarSetMatrixTransformer> updateToSolInUpdate =
       getMethodData().getUpdateToSolutionInUpdateMatTrans();
     
     updateToSolInUpdate->setMatrix(state);
     const RealMatrix& tMatrix = *updateToSolInUpdate->getMatrix();
     
     fluxJacob = _tJacob*tMatrix;
   }
}
   
//////////////////////////////////////////////////////////////////////////////

void WeakSlipWallEuler3DImpl::setFaceNormal(const CFuint faceID,
						RealVector& normal)
{

  DataHandle< InwardNormalsData*> normals = socket_normals.getDataHandle();
  DataHandle<Common::Trio<CFuint, CFuint, Common::SafePtr<Framework::TopologicalRegionSet> > >
    faceNeighCell = socket_faceNeighCell.getDataHandle();

  const CFuint cellTrsID = faceNeighCell[faceID].first;
  const CFuint iFaceLocal = faceNeighCell[faceID].second;
  Common::SafePtr<TopologicalRegionSet> cellTrs = faceNeighCell[faceID].third;
  const CFuint cellLocalID = cellTrs->getLocalGeoID(cellTrsID);

  normal[0] = -normals[cellLocalID]->getFaceNormComp(iFaceLocal, 0);
  normal[1] = -normals[cellLocalID]->getFaceNormComp(iFaceLocal, 1);
  normal[2] = -normals[cellLocalID]->getFaceNormComp(iFaceLocal, 2);
}

//////////////////////////////////////////////////////////////////////////////

void WeakSlipWallEuler3DImpl::configure ( Config::ConfigArgs& args )
{
  FluctuationSplitCom::configure(args);
}
      
//////////////////////////////////////////////////////////////////////////////

  } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
