#include "FVMCC_PolyRecLin.hh"
#include "Framework/VolumeIntegrator.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MeshData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

FVMCC_PolyRecLin::FVMCC_PolyRecLin(const std::string& name) :
  PolyReconstructorLin(name),
  _cells(CFNULL),
  _cellBuilder(),
  socket_nodes("nodes"),
  socket_states("states"),
  socket_gstates("gstates")
{
}

//////////////////////////////////////////////////////////////////////////////

FVMCC_PolyRecLin::~FVMCC_PolyRecLin()
{
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_PolyRecLin::updateWeights()
{
  //nothing to do here
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_PolyRecLin::setup()
{
  PolyReconstructorLin::setup();

  _cells = MeshDataStack::getActive()->getTrs("InnerCells");

  SafePtr<CellTrsGeoBuilder> geoBuilderPtr = _cellBuilder.getGeoBuilder();
  geoBuilderPtr->setDataSockets(socket_states, socket_gstates, socket_nodes);
  _cellBuilder.setup();
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_PolyRecLin::computeLimiters()
{
  DataHandle<CFreal> newLimiter = socket_limiter.getDataHandle();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  RealVector limiter(1.0, nbEqs);
  RealVector limiterTmp(1.0, nbEqs);
  
  const CFreal residual =
    SubSystemStatusStack::getActive()->getResidual();
  
  Common::SafePtr<FVMCC_VolumeIntegrator> integrator = getMethodData().getVolumeIntegrator();

  CellTrsGeoBuilder::GeoData& geoData = _cellBuilder.getDataGE();
  geoData.trs = _cells;
  
  const CFuint nbCells = _cells->getLocalNbGeoEnts();
  for (CFuint iCell = 0; iCell < nbCells; ++iCell) {
    geoData.idx = iCell;
    GeometricEntity *const currCell = _cellBuilder.buildGE();

    const State *const state = currCell->getState(0);
    const GeomEntList *const faces = currCell->getNeighborGeos();
    const CFuint nbFaces = faces->size();
    
    for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
      // compute by shape function-based interpolation the
      // coordinates of all the quadrature points in the face
      integrator->computeCoordinatesAtQuadraturePointsOnGeoEnt
        ((*faces)[iFace],_quadPointCoord[iFace]);
    }
    
    const CFuint stateID = state->getLocalID();
    if (residual > _limitRes) {
      const CFuint startID = stateID*nbEqs;
      for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
	newLimiter[startID + iEq] = 1.0;
      }
      _limiter->limit(_quadPointCoord, currCell, &newLimiter[startID]);
    }
    else {
      if (!_freezeLimiter) {
	// historical modification of the limiter
	_limiter->limit(_quadPointCoord, currCell, &limiter[0]);
	CFuint currID = stateID*nbEqs;
	for (CFuint iVar = 0; iVar < nbEqs; ++iVar, ++currID) {
	  newLimiter[currID] = min(limiter[iVar],newLimiter[currID]);
	}
      }
    }
  }
  
  // release the current GeometricEntity
  _cellBuilder.releaseGE();
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
