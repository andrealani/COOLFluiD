#include "LeastSquareP1PolyRecNEQ2D.hh"
#include "Common/CFLog.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Framework/PhysicalChemicalLibrary.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"
#include "FiniteVolumeNEQ/FiniteVolumeNEQ.hh"
#include "FiniteVolume/CellCenterFVMData.hh"
#include "NavierStokes/EulerTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<LeastSquareP1PolyRecNEQ2D, CellCenterFVMData, 
		       PolyReconstructor<CellCenterFVMData>, 
		       FiniteVolumeNEQModule> 
leastSquareP1PolyRecNEQ2DProvider("LinearLSNEQ2D");

//////////////////////////////////////////////////////////////////////////////

void LeastSquareP1PolyRecNEQ2D::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::vector<string> >
    ("SubOutletTRS","TRSs corresponding to a subsonic outlet.");
}
    
//////////////////////////////////////////////////////////////////////////////

LeastSquareP1PolyRecNEQ2D::LeastSquareP1PolyRecNEQ2D(const std::string& name) :
  LeastSquareP1PolyRec2D(name),
  m_library(CFNULL),
  m_pdata(),
  m_subOutletStates()
{
  addConfigOptionsTo(this);

  m_subOutletTRS = vector<std::string>();
  setParameter("SubOutletTRS",&m_subOutletTRS);
}
      
//////////////////////////////////////////////////////////////////////////////

LeastSquareP1PolyRecNEQ2D::~LeastSquareP1PolyRecNEQ2D()
{
}

//////////////////////////////////////////////////////////////////////////////
      
void LeastSquareP1PolyRecNEQ2D::computeGradients()
{ 
  LeastSquareP1PolyRec2D::computeGradients();
  
  if (m_subOutletTRS.size() > 0) {
    SafePtr<ConvectiveVarSet> updateVarSet = getMethodData().getUpdateVar();
    
    // we override here the value of the partial density gradients for the states 
    // lying on the subsonic outlet boundary
    DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
    DataHandle<vector<State*> > stencil = socket_stencil.getDataHandle();
    DataHandle<CFreal> weights = socket_weights.getDataHandle();
    DataHandle<CFreal> uX = socket_uX.getDataHandle();
    DataHandle<CFreal> uY = socket_uY.getDataHandle();
    
    const CFuint nbStates = states.size();
    const CFuint nbEquations = PhysicalModelStack::getActive()->getNbEq();
    const CFuint nbSpecies = m_library->getNbSpecies();
    
    // special treatment for partial density gradients: they are computed from the pressure gradient
    _lf1 = 0.0;
    _lf2 = 0.0;
    CFuint count = 0;
    CFuint iEdge = 0;
    for(CFuint iState = 0; iState < nbStates; ++iState) {
      const State* const first = states[iState];
      cf_assert(!first->isGhost());
      const CFuint stencilSize = stencil[iState].size();
      
      if (m_subOutletStates[first->getLocalID()]) {
	count++;
	
	// loop over the neighbor cells belonging to the chosen stencil
	for(CFuint in = 0; in < stencilSize; ++in) {
	  const State* const last = stencil[iState][in];
	  const CFuint lastID = (!last->isGhost()) ? last->getLocalID() : numeric_limits<CFuint>::max();
	  const CFuint firstID = first->getLocalID();
	  cf_assert(firstID != lastID);
	  
	  if (lastID > firstID) {
	    // consider the next edge
	    const RealVector& nodeFirst = first->getCoordinates();
	    const RealVector& nodeLast = last->getCoordinates();
	    const CFreal weig = weights[iEdge];
	    const CFreal dx = weig*(nodeLast[0] - nodeFirst[0]);
	    const CFreal dy = weig*(nodeLast[1] - nodeFirst[1]);
	    
	    // compute the pressures here
	    updateVarSet->computePhysicalData((*last), m_pdata);
	    const CFreal pLast = m_pdata[EulerTerm::P];
	    
	    updateVarSet->computePhysicalData((*first), m_pdata);
	    const CFreal pFirst = m_pdata[EulerTerm::P];
	    
	    const CFreal du = weig*(pLast - pFirst);
	    const CFreal dxdu = dx*du;
	    const CFreal dydu = dy*du;
	    
	    _lf1[firstID] += dxdu;
	    _lf2[firstID] += dydu;
	    
	    if (!last->isGhost()) {
	      _lf1[lastID] += dxdu;
	      _lf2[lastID] += dydu;
	    }
	    
	    ++iEdge;
	  }
	}
      }
      else {
	// keep accountancy of the edges
	for(CFuint in = 0; in < stencilSize; ++in) {
	  const State* const last = stencil[iState][in];
	  const CFuint lastID = (!last->isGhost()) ? last->getLocalID() : numeric_limits<CFuint>::max();
	  const CFuint firstID = first->getLocalID();
	  cf_assert(firstID != lastID);
	  
	  if (lastID > firstID) {
	    ++iEdge;
	  }
	}
      }
    }
    
    for(CFuint iState = 0; iState < nbStates; ++iState) {
      if (m_subOutletStates[iState]) {
	const CFreal det = (_l11[iState]*_l22[iState] - _l12[iState]*_l12[iState]);
	cf_assert(std::abs(det) > 0.);
	const CFreal invDet = 1./det;
	const CFreal px = (_l22[iState]*_lf1[iState] - _l12[iState]*_lf2[iState])*invDet;
	const CFreal py = (_l11[iState]*_lf2[iState] - _l12[iState]*_lf1[iState])*invDet;
	updateVarSet->computePhysicalData(*states[iState], m_pdata);
	const CFreal p = m_pdata[EulerTerm::P];
	
	// \f$ \nabla \rho_i  = \frac{\rho_i}{p} \nabla p \f$
	for(CFuint iVar = 0; iVar < nbSpecies; ++iVar) {
	  uX(iState,iVar,nbEquations) = px*(*states[iState])[iVar]/p; 
	  uY(iState,iVar,nbEquations) = py*(*states[iState])[iVar]/p;
	}
      }
    }
  }
}
      
//////////////////////////////////////////////////////////////////////////////

void LeastSquareP1PolyRecNEQ2D::setup()
{
  LeastSquareP1PolyRec2D::setup();
  
  m_library = PhysicalModelStack::getActive()->getImplementor()->
    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
  assert(m_library.isNotNull());
  
  PhysicalModelStack::getActive()->getImplementor()->getConvectiveTerm()->
    resizePhysicalData(m_pdata);
  
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  const CFuint nbStates = states.size();
  m_subOutletStates.resize(nbStates, false);
  
  Common::SafePtr<GeometricEntityPool<FaceTrsGeoBuilder> > geoBuilder = 
    getMethodData().getFaceTrsGeoBuilder();
  
  SafePtr<FaceTrsGeoBuilder> geoBuilderPtr = geoBuilder->getGeoBuilder();
  geoBuilderPtr->setDataSockets(socket_states, socket_gstates, socket_nodes);
  
  FaceTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.isBFace = true;
  
  CFuint count = 0;
  for (CFuint iTRS = 0; iTRS < m_subOutletTRS.size(); ++iTRS) {
    SafePtr<TopologicalRegionSet> trs = MeshDataStack::getActive()->getTrs(m_subOutletTRS[iTRS]);
    geoData.trs = trs;
    
    const CFuint nbTrsFaces = trs->getLocalNbGeoEnts();
    for (CFuint iFace = 0; iFace < nbTrsFaces; ++iFace) {
      // build the GeometricEntity
      geoData.idx = iFace;
      GeometricEntity *const face = geoBuilder->buildGE();
      
      // flag the internal states on subsonic outlet boundaries 
      m_subOutletStates[face->getState(0)->getLocalID()] = true;
      count++;
      
      // release the GeometricEntity
      geoBuilder->releaseGE();
    }
  }
}
 
//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
