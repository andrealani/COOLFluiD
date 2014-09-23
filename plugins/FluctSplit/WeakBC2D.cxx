#include "WeakBC2D.hh"
#include "InwardNormalsData.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MeshData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

WeakBC2D::WeakBC2D(const std::string& name) :
  WeakBC(name),
  socket_pastNormals("pastNormals", false),
  socket_pastStates("pastStates", false),
  socket_isUpdated("isUpdated"),
  m_gstate(CFNULL),
  m_state(CFNULL)
{
}

//////////////////////////////////////////////////////////////////////////////

WeakBC2D::~WeakBC2D()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
WeakBC2D::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result = WeakBC::needsSockets();

  result.push_back(&socket_pastNormals);
  result.push_back(&socket_pastStates);
  result.push_back(&socket_isUpdated);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void WeakBC2D::setup()
{
  WeakBC::setup();

  m_gstate = new State();
  m_state = new State();

Common::SafePtr<GeometricEntityPool<StdTrsGeoBuilder> >
  geoBuilder = getMethodData().getStdTrsGeoBuilder();

  StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();

  std::vector< Common::SafePtr<TopologicalRegionSet> >& trss = getTrsList();

  vector< SafePtr<TopologicalRegionSet> >::iterator itrs = trss.begin();
  CFuint max_nb_states = 0;
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
      geoBuilder->releaseGE();
    }
  }
  m_flux.resize(max_nb_states);
const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
 for ( CFuint iState = 0; iState < max_nb_states; ++iState)
         m_flux[iState].resize(nbEqs);
}

//////////////////////////////////////////////////////////////////////////////

void WeakBC2D::unsetup()
{
  deletePtr(m_gstate);
  deletePtr(m_state);

  WeakBC::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void WeakBC2D::executeOnTrs()
{
  CFAUTOTRACE;
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
//  const CFuint nbTRs = getCurrentTRS()->getNbTRs();
  const CFreal oEminusAlpha = 1. - m_alpha;
  // states forming the ghost cell
  vector<State*> states(2);

  bool isUnsteady = false;

  CFreal dt = SubSystemStatusStack::getActive()->getDT();
  if(dt < 0.) dt = 1.;
  else isUnsteady = true;

  // get the data handle for the states
  DataHandle<State*,GLOBAL> allStates = socket_states.getDataHandle();
  // get the data handle for the is Updated flags
  DataHandle<bool> isUpdated = socket_isUpdated.getDataHandle();
  // get the data handle for the rhs
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();

  Common::SafePtr<GeometricEntityPool<StdTrsGeoBuilder> >
  geoBuilder = getMethodData().getStdTrsGeoBuilder();

  StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.trs = getCurrentTRS();
  const CFuint nbFaces = getCurrentTRS()->getLocalNbGeoEnts();

  for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
    CFLogDebugMed( "Computing iFace = " << iFace << "\n");



     // build the GeometricEntity
    geoData.idx = iFace;
    GeometricEntity& currFace = *geoBuilder->buildGE();

    const CFuint nb_states = currFace.nbStates();



    // set the face normal
    if (SubSystemStatusStack::getActive()->isMovingMesh()){
      setMovingFaceNormal(currFace.getID());
    }
    else{
      setFaceNormal(currFace.getID());
    }

    for (CFuint iState = 0 ; iState < nb_states; ++iState) {
    // UNSTEADY

        CFuint stateID = currFace.getState(iState)->getLocalID();

        if (SubSystemStatusStack::getActive()->isMovingMesh()){
           DataHandle<State*> pastStates = socket_pastStates.getDataHandle();

             // Unsteady, state=0.5*(pastState+interState)

               for (CFuint j=0; j<nbEqs ; ++j){

               (*m_state)[j] = 0.5*((*allStates[stateID])[j] + (*pastStates[stateID])[j]);
           }

           RealVector stateCoord = ((*allStates[stateID]).getCoordinates()) ;
           stateCoord += ((*pastStates[stateID]).getCoordinates());
           stateCoord *= 0.5;

           bool isOnMesh = false;
           Node* coord = new Node(stateCoord, isOnMesh);

           m_state->setSpaceCoordinates(coord);
        }
        else {
           // clone the states in order to use the indexes

          m_state->clone(*allStates[stateID]);
	  m_gstate->clone(*allStates[stateID]);	
          }

           // set the appropriate values in the ghost states
          setGhostState((*m_state), (*m_gstate));

         // set ghost state == first state
         // corresponding B-state == second state

         states[0] = m_gstate;
         states[1] = m_state;

        // compute the corresponding fluctuation and the update coeff
        computeFlux(states, m_flux[iState]);

        m_flux[iState] *= dt;
    }

    for (CFuint iEq = 0; iEq < nbEqs; ++iEq)
    {
      // distribute corrections to the boundary nodes
      for (CFuint iState = 0 ; iState < nb_states; ++iState) {
          CFuint stateID = currFace.getState(iState)->getLocalID();
          CFreal sum_flux_not_iState = 0.;
	 for (CFuint jState = 0; jState < nb_states ; ++jState){
                       if (jState != iState)
                           sum_flux_not_iState += m_flux[jState][iEq];
                  }
          if (!isUpdated[stateID]) {
             rhs(stateID, iEq, nbEqs) +=
                 (m_alpha*m_flux[iState][iEq] + oEminusAlpha*sum_flux_not_iState);

             m_flagState[stateID] = true;
          }
      }

    }
//This is the true implementation for High order
//This is put like that for the moment to see if it is really usefull.
//If it appears to be usefull, it will be trully implemented, doing a template

//   for (CFuint iEq = 0; iEq < nbEqs; ++iEq)
//     {
//
//            CFuint stateID = currFace.getState(0)->getLocalID();
//
//            if (!isUpdated[stateID]) {
//               rhs(stateID, iEq, nbEqs) +=  (5.0/12.0)*m_flux[0][iEq];
//
//               m_flagState[stateID] = true;
//            }
//
//           stateID = currFace.getState(2)->getLocalID();
//
//            if (!isUpdated[stateID]) {
//               rhs(stateID, iEq, nbEqs) +=  (1.0/3.0)*m_flux[2][iEq];
//
//               m_flagState[stateID] = true;
//            }
//
//           stateID = currFace.getState(1)->getLocalID();
//
//            if (!isUpdated[stateID]) {
//               rhs(stateID, iEq, nbEqs) +=  (5.0/12.0)*m_flux[1][iEq];
//
//               m_flagState[stateID] = true;
//            }
//
//        }


  //Release the geometric entity
  geoBuilder->releaseGE();
  }

  // this could be done in a more efficient way ...
  const CFuint nbStates = m_flagState.size();
  for (CFuint i = 0; i < nbStates; ++i) {
    if (m_flagState[i] == true) {
      isUpdated[i] = true;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void WeakBC2D::setMovingFaceNormal(const CFuint faceID)
{
  DataHandle< InwardNormalsData*> normals = socket_normals.getDataHandle();
  DataHandle< InwardNormalsData*> pastNormals = socket_pastNormals.getDataHandle();
  DataHandle<Common::Trio<CFuint, CFuint, Common::SafePtr<Framework::TopologicalRegionSet> > >
    faceNeighCell = socket_faceNeighCell.getDataHandle();

  const CFuint cellTrsID = faceNeighCell[faceID].first;
  const CFuint iFaceLocal = faceNeighCell[faceID].second;
  Common::SafePtr<TopologicalRegionSet> cellTrs = faceNeighCell[faceID].third;
  const CFuint cellLocalID = cellTrs->getLocalGeoID(cellTrsID);

  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  for (CFuint iDim = 0; iDim < dim; ++iDim) {
    // external normal is already inward for ghost state
    // The normal is the average of the past and current
    // normals...
    m_faceNormal[iDim] = -0.5*(normals[cellLocalID]->
			      getFaceNormComp(iFaceLocal, iDim)
			      +pastNormals[cellLocalID]->
			      getFaceNormComp(iFaceLocal, iDim));
  }
}

//////////////////////////////////////////////////////////////////////////////

} // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
