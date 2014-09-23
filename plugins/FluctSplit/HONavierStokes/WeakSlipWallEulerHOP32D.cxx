#include "FluctSplit/FluctSplitNavierStokes.hh"
#include "WeakSlipWallEulerHOP32D.hh"
#include "FluctSplit/InwardNormalsData.hh"
#include "Framework/CFL.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MethodCommandProvider.hh"
#include "NavierStokes/Euler2DVarSet.hh"
#include "Framework/MeshData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<WeakSlipWallEulerHOP32D,
		      FluctuationSplitData,
		      FluctSplitNavierStokesModule>
WeakSlipWallEulerHOP32DProvider("WeakSlipWallEulerHOP32D");

//////////////////////////////////////////////////////////////////////////////

void WeakSlipWallEulerHOP32D::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("alpha","distribution coefficient");
   options.addConfigOption< bool >("ExactNormNaca0012","Use the exact normals for a NACA0012");
}

//////////////////////////////////////////////////////////////////////////////

WeakSlipWallEulerHOP32D::WeakSlipWallEulerHOP32D(const std::string& name) :
  FluctuationSplitCom(name),
  socket_rhs("rhs"),
  socket_states("states"),
  socket_isUpdated("isUpdated"),
  socket_normals("normals"),
  socket_faceNeighCell("faceNeighCell"),
  _flagState(),
  _varSet(CFNULL),
  _physicalData(),
  _tmpFlux(),
  m_qdNormal()
{
   addConfigOptionsTo(this);
  _alpha = 1.0;
   setParameter("alpha",&_alpha);


  m_exact_norm = false;
  setParameter("ExactNormNaca0012",&m_exact_norm);
}

//////////////////////////////////////////////////////////////////////////////

WeakSlipWallEulerHOP32D::~WeakSlipWallEulerHOP32D()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
WeakSlipWallEulerHOP32D::needsSockets()
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

void WeakSlipWallEulerHOP32D::setup()
{
  _varSet = getMethodData().getUpdateVar().d_castTo<Euler2DVarSet>();

  // get the data handle for the states
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  _flagState.resize(states.size());
  _flagState = false;

  // set up the physical data
  _varSet->getModel()->resizePhysicalData(_physicalData);

  _tmpFlux.resize(PhysicalModelStack::getActive()->getNbEq());
 CFuint max_nb_states = 0;
 Common::SafePtr<GeometricEntityPool<StdTrsGeoBuilder> >
  geoBuilder = getMethodData().getStdTrsGeoBuilder();

  StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();

  std::vector< Common::SafePtr<TopologicalRegionSet> >& trss = getTrsList();

  std::vector< SafePtr<TopologicalRegionSet> >::iterator itrs = trss.begin();


  for (; itrs != trss.end(); ++itrs)
  {
    SafePtr<TopologicalRegionSet> faces = *itrs;

    geoData.trs = faces;

    const CFuint nbFaces = faces->getLocalNbGeoEnts();

    for (CFuint iFace = 0; iFace < nbFaces; ++iFace)
    {
      // build the GeometricEntity
      geoData.idx = iFace;

      GeometricEntity& currFace = *geoBuilder->buildGE();

      const CFuint nb_state = currFace.nbStates();

      max_nb_states = std::max(max_nb_states, nb_state);

    // release the face
    geoBuilder->releaseGE();
    }
  }

// vector of fluxes of each state
    m_flux.resize(max_nb_states);
const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
   for ( CFuint iState = 0; iState < max_nb_states; ++iState)
      m_flux[iState].resize(nbEqs);
  // vector of states of the face
    m_states.resize(max_nb_states);
}

//////////////////////////////////////////////////////////////////////////////

void WeakSlipWallEulerHOP32D::executeOnTrs()
{
  CFAUTOTRACE;
  DataHandle<bool> isUpdated = socket_isUpdated.getDataHandle();
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();

  RealVector normal(0.0, dim);
  const CFreal oEminusAlpha = 1. - _alpha;


 Common::SafePtr<GeometricEntityPool<StdTrsGeoBuilder> >
    geoBuilder = getMethodData().getStdTrsGeoBuilder();

  StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.trs = getCurrentTRS();

  const CFuint nbFaces = getCurrentTRS()->getLocalNbGeoEnts();
  for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {

    geoData.idx = iFace;
    GeometricEntity& face = *geoBuilder->buildGE();

    // set number of states on one face
    const CFuint nb_state = face.nbStates();

    // set the face normal
    setFaceNormal(face.getID(), normal);

    m_states = *face.getStates();
    m_qdNormal.resize(dim);
    // compute the normal fluxes corrections for both the states
    // of this cell
    for (CFuint iState = 0; iState < nb_state; ++iState){ 

    if (m_exact_norm){
      Node& node0 = (*m_states[iState]).getCoordinates();

      // Coordinate of the vertex
      const CFreal x0 = node0[XX];
      const CFreal y0 = node0[YY];

      if ((y0 >= 0.0) && (x0 != 0.0) && (x0 < 1.0)){
        m_qdNormal[XX] = -0.6*(( 0.2969/(2.0*sqrt(x0))) - 0.1260 - 0.3516*2.0*x0 + 3.0*0.2843*x0*x0 - 0.1015*4.0*x0*x0*x0 );
        m_qdNormal[YY] = 1.0;
      }
      else if ((y0 <= 0.0) && (x0 != 0.0) && (x0 < 1.0)){
        m_qdNormal[XX] = -0.6*(( 0.2969/(2.0*sqrt(x0))) - 0.1260 - 0.3516*2.0*x0 + 3.0*0.2843*x0*x0 - 0.1015*4.0*x0*x0*x0 );
        m_qdNormal[YY] = -1.0;
      }
      else if ((x0 == 0.0)){
        m_qdNormal[XX] = -1.0;
        m_qdNormal[YY] = 0.0;
      }
      else {m_qdNormal = normal;}
      m_qdNormal /= m_qdNormal.norm2();
      m_qdNormal *= normal.norm2();
      }
      
      else {m_qdNormal = normal; }
         computeNormalFlux(*m_states[iState], m_qdNormal, m_flux[iState]);

}

    // distribute contributions to the two nodes
    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      if (!isUpdated[m_states[0]->getLocalID()]){
        rhs(m_states[0]->getLocalID(), iEq, nbEqs) -= _alpha*(1.0/8.0)*m_flux[0][iEq] + (oEminusAlpha/3.0)*((1.0/8.0)*m_flux[1][iEq] + (3.0/8.0)*m_flux[2][iEq] + (3.0/8.0)*m_flux[3][iEq]);

      }

      if (!isUpdated[m_states[1]->getLocalID()]){
        rhs(m_states[1]->getLocalID(), iEq, nbEqs) -= _alpha*(1.0/8.0)*m_flux[1][iEq] + (oEminusAlpha/3.0)*((1.0/8.0)*m_flux[0][iEq] + (3.0/8.0)*m_flux[2][iEq] + (3.0/8.0)*m_flux[3][iEq]);

      }

      if (!isUpdated[m_states[2]->getLocalID()]){
        rhs(m_states[2]->getLocalID(), iEq, nbEqs) -= _alpha*(3.0/8.0)*m_flux[2][iEq] + (oEminusAlpha/3.0)*((1.0/8.0)*m_flux[1][iEq] + (1.0/8.0)*m_flux[0][iEq] + (3.0/8.0)*m_flux[3][iEq]);

      }

  if (!isUpdated[m_states[3]->getLocalID()]){
        rhs(m_states[3]->getLocalID(), iEq, nbEqs) -= _alpha*(3.0/8.0)*m_flux[3][iEq] + (oEminusAlpha/3.0)*((1.0/8.0)*m_flux[1][iEq] + (1.0/8.0)*m_flux[0][iEq] + (3.0/8.0)*m_flux[2][iEq]);

      }
}
   geoBuilder->releaseGE();
  }

}

//////////////////////////////////////////////////////////////////////////////

void WeakSlipWallEulerHOP32D::computeNormalFlux(const State& state,
                                            const RealVector& normal,
                                            RealVector& flux)
{
  State& ss = *(const_cast<State*>(&state));
  _varSet->computePhysicalData(ss, _physicalData);

  const CFreal un = _physicalData[EulerTerm::VX]*normal[XX] +
    _physicalData[EulerTerm::VY]*normal[YY];

  const CFreal rho = _physicalData[EulerTerm::RHO];

  //if (!getMethodData().isResidualTransformationNeeded()) {
  flux[0] = rho*un;
  flux[1] = un*rho*_physicalData[EulerTerm::VX];
  flux[2] = un*rho*_physicalData[EulerTerm::VY];
  flux[3] = un*rho*_physicalData[EulerTerm::H];
  //}
  //   else {
  //     _tmpFlux[0] = rho*un;
  //     _tmpFlux[1] = un*rho*_physicalData[EulerTerm::VX];
  //     _tmpFlux[2] = un*rho*_physicalData[EulerTerm::VY];
  //     _tmpFlux[3] = un*rho*_physicalData[EulerTerm::H];

  //     // set the transformation from update to solution in update
  //     SafePtr<VarSetMatrixTransformer> solToUpdateInUpdate =
  //       getMethodData().getSolToUpdateInUpdateMatTrans();

  //     solToUpdateInUpdate->setMatrix(state);
  //     const RealMatrix& tMatrix = *solToUpdateInUpdate->getMatrix();
  //     flux = tMatrix*_tmpFlux;
  //   }
}

//////////////////////////////////////////////////////////////////////////////

void WeakSlipWallEulerHOP32D::setFaceNormal(const CFuint faceID,
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
}

//////////////////////////////////////////////////////////////////////////////

} // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
