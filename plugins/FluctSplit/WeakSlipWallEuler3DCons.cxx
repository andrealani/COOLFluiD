#include "FluctSplit/FluctSplitNavierStokes.hh"
#include "WeakSlipWallEuler3DCons.hh"
#include "InwardNormalsData.hh"
#include "Framework/CFL.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/NamespaceSwitcher.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<WeakSlipWallEuler3DCons, FluctuationSplitData, FluctSplitNavierStokesModule> weakSlipWallEuler3DConsProvider("WeakSlipWallEuler3DCons");

//////////////////////////////////////////////////////////////////////////////

void WeakSlipWallEuler3DCons::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("alpha","distribution coefficient");
}

//////////////////////////////////////////////////////////////////////////////

WeakSlipWallEuler3DCons::WeakSlipWallEuler3DCons(const std::string& name) :
  FluctuationSplitCom(name),
  socket_rhs("rhs"),
  socket_states("states"),
  socket_isUpdated("isUpdated"),
  socket_normals("normals"),
  socket_faceNeighCell("faceNeighCell"),
  _varSet(),
  _fluxes(),
  _flagState()
{
   addConfigOptionsTo(this);
  _alpha = 1.0;
   setParameter("alpha",&_alpha);
}

//////////////////////////////////////////////////////////////////////////////

WeakSlipWallEuler3DCons::~WeakSlipWallEuler3DCons()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
WeakSlipWallEuler3DCons::needsSockets()
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

void WeakSlipWallEuler3DCons::setup()
{

  _varSet->setup();
  // size the fluxes
  _fluxes.resize(MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell()-1);
  for (CFuint i = 0; i < _fluxes.size(); ++i) {
    _fluxes[i].resize(PhysicalModelStack::getActive()->getNbEq());
  }
 
  // get the data handle for the states
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

  _flagState.resize(states.size());
  _flagState = false;
  
 }

//////////////////////////////////////////////////////////////////////////////

void WeakSlipWallEuler3DCons::executeOnTrs()
{
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  RealVector normal(0.0, dim);
  const CFreal oEminusAlpha = 1. - _alpha;
  const CFreal third = 1./3;

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<bool> isUpdated = socket_isUpdated.getDataHandle();
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();

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
      computeNormalFlux(*(*statesInFace)[is],
			normal,
			_fluxes[is]);
    }

    /// @todo implement case when the face is quadrilateral
    if (nbStatesInFace == 3) {
      const CFuint stateID0 = (*statesInFace)[0]->getLocalID();
      const CFuint stateID1 = (*statesInFace)[1]->getLocalID();
      const CFuint stateID2 = (*statesInFace)[2]->getLocalID();

      // distribute contributions to the nodes
      for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
	if (!isUpdated[stateID0]) {
	  rhs(stateID0, iEq, nbEqs) -=
	    third*(_alpha*_fluxes[0][iEq] + oEminusAlpha*0.5*
		   (_fluxes[1][iEq] + _fluxes[2][iEq]));
	  _flagState[stateID0] = true;
	}

	if (!isUpdated[stateID1]) {
	  rhs(stateID1, iEq, nbEqs) -=
	    third*(_alpha*_fluxes[1][iEq] + oEminusAlpha*0.5*
		   (_fluxes[0][iEq] + _fluxes[2][iEq]));
	  _flagState[stateID1] = true;
	}

	if (!isUpdated[stateID2]) {
	  rhs(stateID2, iEq, nbEqs) -=
	    third*(_alpha*_fluxes[2][iEq] + oEminusAlpha*0.5*
		   (_fluxes[0][iEq] + _fluxes[1][iEq]));

	  _flagState[stateID2] = true;
	}
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

void WeakSlipWallEuler3DCons::computeNormalFlux(const RealVector& state,
                                                const RealVector& normal,
                                                RealVector& flux) const
 {
   const CFreal nx = normal[0];
   const CFreal ny = normal[1];
   const CFreal nz = normal[2];
   const CFreal rho = state[0];
   const CFreal u = state[1]/rho;
   const CFreal v = state[2]/rho;
   const CFreal w = state[3]/rho;
   const CFreal un = u*nx + v*ny + w*nz;
   const CFreal gamma = _varSet->getModel()->getGamma();
   const CFreal halfGammaMinus1 = 0.5*(gamma - 1);
   const CFreal rhoV2 = rho*(u*u + v*v + w*w);

   flux[0] = rho*un;
   flux[1] = un*state[1];
   flux[2] = un*state[2];
   flux[3] = un*state[3];
   flux[4] = un*(gamma*state[4] - halfGammaMinus1*rhoV2);
 }

//////////////////////////////////////////////////////////////////////////////

void WeakSlipWallEuler3DCons::setFaceNormal(const CFuint faceID,
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

void WeakSlipWallEuler3DCons::configure ( Config::ConfigArgs& args )
{
  FluctuationSplitCom::configure(args);

  std::string name = getMethodData().getNamespace();
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance
    (SubSystemStatusStack::getCurrentName()).getNamespace(name);
  Common::SafePtr<PhysicalModel> physModel = PhysicalModelStack::getInstance().getEntryByNamespace(nsp);
  
  std::string varSetName = "Euler3DCons";
  _varSet.reset((Environment::Factory<ConvectiveVarSet>::getInstance().getProvider(varSetName)->
    create(physModel->getImplementor()->getConvectiveTerm())).d_castTo<Physics::NavierStokes::Euler3DCons>());

  cf_assert(_varSet.isNotNull());

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
