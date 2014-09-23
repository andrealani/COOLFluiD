#include "FiniteVolumeMHD/FiniteVolumeMHD.hh"
#include "MirrorMHD3DProjectionTanakaLinearDensityExtrapolation.hh"
#include "Framework/MethodCommandProvider.hh"
#include "MHD/MHD3DProjectionVarSet.hh"
#include "Framework/MeshData.hh"

//////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::MHD;

//////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////

MethodCommandProvider<MirrorMHD3DProjectionTanakaLinearDensityExtrapolation, CellCenterFVMData, FiniteVolumeMHDModule>
mirrorMHD3DProjectionTanakaLinearDensityExtrapolationFVMCCProvider("MirrorMHD3DProjectionTanakaLinearDensityExtrapolationFVMCC");

//////////////////////////////////////////////////////////////////////

MirrorMHD3DProjectionTanakaLinearDensityExtrapolation::MirrorMHD3DProjectionTanakaLinearDensityExtrapolation(const std::string& name) :
  FVMCC_BC(name),
  socket_stencil("stencil"),
  _alignedNgbourStateID(CFNULL),
  _rNgbourState(CFNULL),
  _placeInStencil(CFNULL),
  _varSet(CFNULL),
  _dataInnerState(),
  _dataGhostState()
{
}

//////////////////////////////////////////////////////////////////////

MirrorMHD3DProjectionTanakaLinearDensityExtrapolation::~MirrorMHD3DProjectionTanakaLinearDensityExtrapolation()
{
}

//////////////////////////////////////////////////////////////////////

void MirrorMHD3DProjectionTanakaLinearDensityExtrapolation::setup()
{
 FVMCC_BC::setup();

 SafePtr<TopologicalRegionSet> trs = getCurrentTRS();
 CFLogDebugMin( "FVMCC_BC::execute() called for TRS: "
		       << trs->getName() << "\n");

 Common::SafePtr<GeometricEntityPool<FaceTrsGeoBuilder> >
	 geoBuilder = getMethodData().getFaceTrsGeoBuilder();

 SafePtr<FaceTrsGeoBuilder> geoBuilderPtr = geoBuilder->getGeoBuilder();
 geoBuilderPtr->setDataSockets(socket_states, socket_gstates,
			 socket_nodes);

 FaceTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
 geoData.isBFace = true;
 geoData.trs = trs;

 DataHandle<vector<State*> > stencil = socket_stencil.getDataHandle();

 const CFuint nbTrsFaces = trs->getLocalNbGeoEnts();
 const CFuint nbStates = MeshDataStack::getActive()->getNbStates();

 _alignedNgbourStateID.resize(nbStates);
 _rNgbourState.resize(nbStates);
 _placeInStencil.resize(nbStates);

 CFreal rStateTimesrNgbourStaemin = 0.;
 CFreal rNgbourState = 0.;
 CFuint alignedNgbourStateID = 0;
 CFuint placeInStencil = 0;

 for (CFuint iFace = 0; iFace < nbTrsFaces; ++iFace) {
   CFLogDebugMed( "iFace = " << iFace << "\n");

   // build the GeometricEntity
   geoData.idx = iFace;

   GeometricEntity *const face = geoBuilder->buildGE();

   const CFuint stateID = face->getState(LEFT)->getLocalID();
   RealVector stateCoord = face->getState(LEFT)->getCoordinates();

   const CFuint stencilSize = stencil[stateID].size();

   // finding the cell in the stencil of the boundary cell
   // whose cell center is the closest to the line passing
   // through the boundary cell center and the origin
   // (the center of the planet assumed to be located at the origin)

   rStateTimesrNgbourStaemin = 1000000.0;
   // loop over the neighbor cells belonging to the chosen stencil
   for (CFuint in = 0; in < stencilSize; ++in) {
     const State* const ngbourState = stencil[stateID][in];
     RealVector ngbourStateCoord = ngbourState->getCoordinates();
     const CFreal rStateTimesrNgbourState =
	     sqrt((stateCoord[1]*ngbourStateCoord[2]-stateCoord[2]*ngbourStateCoord[1])*
		  (stateCoord[1]*ngbourStateCoord[2]-stateCoord[2]*ngbourStateCoord[1])+
		  (stateCoord[2]*ngbourStateCoord[0]-stateCoord[0]*ngbourStateCoord[2])*
		  (stateCoord[2]*ngbourStateCoord[0]-stateCoord[0]*ngbourStateCoord[2])+
		  (stateCoord[0]*ngbourStateCoord[1]-stateCoord[1]*ngbourStateCoord[0])*
		  (stateCoord[0]*ngbourStateCoord[1]-stateCoord[1]*ngbourStateCoord[0]));
     if (rStateTimesrNgbourState<rStateTimesrNgbourStaemin) {
 	rStateTimesrNgbourStaemin = rStateTimesrNgbourState;
	alignedNgbourStateID = ngbourState->getLocalID();
	rNgbourState = sqrt(ngbourStateCoord[0]*ngbourStateCoord[0]+
			ngbourStateCoord[1]*ngbourStateCoord[1]+
			ngbourStateCoord[2]*ngbourStateCoord[2]);
	placeInStencil = in;
     }
   }

   _alignedNgbourStateID[stateID] = alignedNgbourStateID;
   _rNgbourState[alignedNgbourStateID] = rNgbourState;
   _placeInStencil[alignedNgbourStateID] = placeInStencil;

   // release the GeometricEntity
   geoBuilder->releaseGE();
 }

 _varSet = getMethodData().getUpdateVar().d_castTo<MHD3DProjectionVarSet>();
 _varSet->getModel()->resizePhysicalData(_dataInnerState);
 _varSet->getModel()->resizePhysicalData(_dataGhostState);
}

//////////////////////////////////////////////////////////////////////

void MirrorMHD3DProjectionTanakaLinearDensityExtrapolation::setGhostState(GeometricEntity *const face)
{
  State *const innerState = face->getState(0);
  State *const ghostState = face->getState(1);
  const CFuint faceID = face->getID();
  const CFuint startID = faceID*PhysicalModelStack::getActive()->getDim();

  DataHandle<CFreal> normals = socket_normals.getDataHandle();

  CFreal nx = normals[startID];
  CFreal ny = normals[startID + 1];
  CFreal nz = normals[startID + 2];
  const CFreal invFaceLength = 1./sqrt(nx*nx + ny*ny + nz*nz);
  nx *= invFaceLength;
  ny *= invFaceLength;
  nz *= invFaceLength;

  // set the physical data starting from the inner state
  _varSet->computePhysicalData(*innerState, _dataInnerState);

  const CFreal bn = _dataInnerState[MHDTerm::BX]*nx +
    _dataInnerState[MHDTerm::BY]*ny +
    _dataInnerState[MHDTerm::BZ]*nz;

  RealVector innerStateCoord = innerState->getCoordinates();
  RealVector ghostStateCoord = ghostState->getCoordinates();

  // the center of the planet assumed to be located at the origin

  const CFreal rInnerState = sqrt(innerStateCoord[0]*innerStateCoord[0]+
		  innerStateCoord[1]*innerStateCoord[1]+
		  innerStateCoord[2]*innerStateCoord[2]);
  const CFreal rGhostState = sqrt(ghostStateCoord[0]*ghostStateCoord[0]+
		  ghostStateCoord[1]*ghostStateCoord[1]+
		  ghostStateCoord[2]*ghostStateCoord[2]);

  const CFuint innerStateID = innerState->getLocalID();

  DataHandle<vector<State*> > stencil = socket_stencil.getDataHandle();

  // all boundary conditions except for density are according to Gabor Toth's suggestions

  _dataGhostState[MHDTerm::RHO] = (*stencil[innerStateID][_placeInStencil[_alignedNgbourStateID[innerStateID]]])[0] +
	  ((rGhostState-_rNgbourState[_alignedNgbourStateID[innerStateID]])/
	  (rInnerState-_rNgbourState[_alignedNgbourStateID[innerStateID]]))*
	  (_dataInnerState[MHDTerm::RHO]-(*stencil[innerStateID][_placeInStencil[_alignedNgbourStateID[innerStateID]]])[0]);
  _dataGhostState[MHDTerm::VX] = -_dataInnerState[MHDTerm::VX];
  _dataGhostState[MHDTerm::VY] = -_dataInnerState[MHDTerm::VY];
  _dataGhostState[MHDTerm::VZ] = -_dataInnerState[MHDTerm::VZ];
  _dataGhostState[MHDTerm::BX] = _dataInnerState[MHDTerm::BX] - 2.0*bn*nx;
  _dataGhostState[MHDTerm::BY] = _dataInnerState[MHDTerm::BY] - 2.0*bn*ny;
  _dataGhostState[MHDTerm::BZ] = _dataInnerState[MHDTerm::BZ] - 2.0*bn*nz;
  _dataGhostState[MHDTerm::V] = _dataInnerState[MHDTerm::V];
  _dataGhostState[MHDTerm::P] = _dataInnerState[MHDTerm::P];
  _dataGhostState[MHDTerm::A] = sqrt(_varSet->getModel()->getGamma()*_dataGhostState[MHDTerm::P]/_dataGhostState[MHDTerm::RHO]);
  _dataGhostState[MHDTerm::B] = _dataInnerState[MHDTerm::B];
  _dataGhostState[MHDProjectionTerm::PHI] = _dataInnerState[MHDProjectionTerm::PHI];

  _varSet->computeStateFromPhysicalData(_dataGhostState, *ghostState);

}

//////////////////////////////////////////////////////////////////////////////

std::vector<SafePtr<BaseDataSocketSink> > MirrorMHD3DProjectionTanakaLinearDensityExtrapolation::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result = FVMCC_BC::needsSockets();
  result.push_back(&socket_stencil);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
