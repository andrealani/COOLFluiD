#include "WenoP1PolyRec2D.hh"
#include "Framework/VolumeIntegrator.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Framework/PhysicalModel.hh"
#include "Common/CFLog.hh"
#include "Framework/MeshData.hh"
#include "FiniteVolume/FiniteVolume.hh"
#include "FiniteVolume/CellCenterFVMData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<WenoP1PolyRec2D,
                       CellCenterFVMData,
                       PolyReconstructor<CellCenterFVMData>,
                       FiniteVolumeModule>
                       wenoP1PolyRec2DProvider("LinearWenoLS2D");

//////////////////////////////////////////////////////////////////////////////

WenoP1PolyRec2D::WenoP1PolyRec2D(const std::string& name) :
  LeastSquareP1PolyRec2D(name),
  socket_dr("dr")
{
}
      
//////////////////////////////////////////////////////////////////////////////

WenoP1PolyRec2D::~WenoP1PolyRec2D()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> > WenoP1PolyRec2D::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result = LeastSquareP1PolyRec2D::needsSockets();

  // Add the needed DataSocketSinks
  result.push_back(&socket_dr);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void WenoP1PolyRec2D::computeGradients()
{
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<vector<State*> > stencil = socket_stencil.getDataHandle();
  DataHandle<CFreal> dr = socket_dr.getDataHandle();
  DataHandle<CFreal> uX = socket_uX.getDataHandle();
  DataHandle<CFreal> uY = socket_uY.getDataHandle();

  const CFuint nbStates = states.size();
  const CFuint nbEquations = PhysicalModelStack::getActive()->getNbEq();

  for(CFuint iVar = 0; iVar < nbEquations; ++iVar) {
    _l11 = 0.0;
    _l12 = 0.0;
    _l22 = 0.0;
    _lf1 = 0.0;
    _lf2 = 0.0;
    CFuint iEdge = 0;
    
    for(CFuint iState = 0; iState < nbStates; ++iState) {
      const State* const first = states[iState];
      const CFuint stencilSize = stencil[iState].size();
      // loop over the neighbor cells belonging to the chosen stencil
      for(CFuint in = 0; in < stencilSize; ++in) {
	const State* const last = stencil[iState][in];
	const CFuint lastID = (!last->isGhost()) ? last->getLocalID() : 
	  numeric_limits<CFuint>::max();
	const CFuint firstID = first->getLocalID();
	cf_assert(firstID != lastID);
	
	if (lastID > firstID) {
	  // consider the next edge
	  const RealVector& nodeFirst = first->getCoordinates();
	  const RealVector& nodeLast = last->getCoordinates();
	  
	  const CFreal diffU = (*last)[iVar] - (*first)[iVar];
	  const CFreal weight = 1./sqrt(diffU*diffU + dr[iEdge]*dr[iEdge]);
	  const CFreal dx = weight*(nodeLast[XX] - nodeFirst[XX]);
	  const CFreal dy = weight*(nodeLast[YY] - nodeFirst[YY]);
	  const CFreal du = weight*diffU;
	  CFLogDebugMin( "du = " << du << "\n");
	  
	  _lf1[firstID] += dx*du;
	  _lf2[firstID] += dy*du;
	  
	  if (!last->isGhost()) {
	    _lf1[lastID] += dx*du;
	    _lf2[lastID] += dy*du;
	  }
	  
	  _l11[firstID] += dx*dx;
	  _l12[firstID] += dx*dy;
	  _l22[firstID] += dy*dy;
	  _l11[lastID] += dx*dx;
	  _l12[lastID] += dx*dy;
	  _l22[lastID] += dy*dy;
	  ++iEdge;
	}
      }
    }
    
    for(CFuint iState = 0; iState < nbStates; ++iState) {
      const CFreal invDet = 1./(_l11[iState]*_l22[iState] - _l12[iState]*_l12[iState]);
      uX(iState,iVar,nbEquations) = (_l22[iState]*_lf1[iState] - _l12[iState]*_lf2[iState])*invDet;
      uY(iState,iVar,nbEquations) = (_l11[iState]*_lf2[iState] - _l12[iState]*_lf1[iState])*invDet;
      
      CFLogDebugMax( "invDet = " << invDet
		     << ", l11 = " << _l11[iState]
		     << ", l12 = " << _l12[iState]
		     << ", l22 = " << _l22[iState]
		     <<", lf1 = " << _lf1[iState]
		     <<", lf2 = " << _lf2[iState]
		     << ", uX =" << uX(iState,iVar,nbEquations)
		     << ", uY =" << uY(iState,iVar,nbEquations)
		     << "\n");
    }
  }
}
      
//////////////////////////////////////////////////////////////////////////////
      
void WenoP1PolyRec2D::updateWeights()
{
  
  FVMCC_PolyRec::updateWeights();

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<vector<State*> > stencil = socket_stencil.getDataHandle();
  DataHandle<CFreal> dr = socket_dr.getDataHandle();

  const CFuint nbStates = states.size();

  _l11 = 0.0;
  _l12 = 0.0;
  _l22 = 0.0;
  _lf1 = 0.0;
  _lf2 = 0.0;

  // weight coefficients are calculated
  CFuint iEdge = 0;
 for(CFuint iState = 0; iState < nbStates; ++iState) {
   const State* const first = states[iState];
   const CFuint stencilSize = stencil[iState].size();

   // loop over the neighbor cells belonging to the chosen stencil
   for(CFuint in = 0; in < stencilSize; ++in) {
     const State* const last = stencil[iState][in];
     const CFuint lastID = last->getLocalID();
     const CFuint firstID = first->getLocalID();
     cf_assert(firstID != lastID);

     if (lastID > firstID) {
       const RealVector& nodeFirst = first->getCoordinates();
       const RealVector& nodeLast = last->getCoordinates();
       const CFreal deltaR = MathFunctions::getDistance(nodeFirst,nodeLast);

       dr[iEdge] = deltaR;
       ++iEdge;
     }
   }
 }

}

//////////////////////////////////////////////////////////////////////////////

void WenoP1PolyRec2D::setup()
{
  FVMCC_PolyRec::setup();

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<vector<State*> > stencil = socket_stencil.getDataHandle();
  DataHandle<CFreal> dr = socket_dr.getDataHandle();

  const CFuint nbStates = states.size();

  _l11.resize(nbStates);
  _l12.resize(nbStates);
  _l22.resize(nbStates);
  _lf1.resize(nbStates);
  _lf2.resize(nbStates);

  // weight coefficients are calculated
  CFuint iEdge = 0;
 for(CFuint iState = 0; iState < nbStates; ++iState) {
   const State* const first = states[iState];
   const CFuint stencilSize = stencil[iState].size();

   // loop over the neighbor cells belonging to the chosen stencil
   for(CFuint in = 0; in < stencilSize; ++in) {
     const State* const last = stencil[iState][in];
     const CFuint lastID = last->getLocalID();
     const CFuint firstID = first->getLocalID();
     cf_assert(firstID != lastID);

     if (lastID > firstID) {
       const RealVector& nodeFirst = first->getCoordinates();
       const RealVector& nodeLast = last->getCoordinates();
       const CFreal deltaR = MathFunctions::getDistance(nodeFirst,nodeLast);

       dr[iEdge] = deltaR;
       ++iEdge;
     }
   }
 }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
