#include "LeastSquareP1PolyRec3D.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Common/CFLog.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"
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

MethodStrategyProvider<LeastSquareP1PolyRec3D, CellCenterFVMData, PolyReconstructor<CellCenterFVMData>, FiniteVolumeModule> leastSquareP1PolyRec3DProvider("LinearLS3D");

//////////////////////////////////////////////////////////////////////////////

LeastSquareP1PolyRec3D::LeastSquareP1PolyRec3D(const std::string& name) :
  FVMCC_PolyRec(name),
  socket_stencil("stencil"),
  socket_weights("weights"),
  socket_uX("uX"),
  socket_uY("uY"),
  socket_uZ("uZ"),
  _l11(),
  _l12(),
  _l13(),
  _l22(),
  _l23(),
  _l33(),
  _lf1(),
  _lf2(),
  _lf3()
{
}

//////////////////////////////////////////////////////////////////////////////

LeastSquareP1PolyRec3D::~LeastSquareP1PolyRec3D()
{
}

//////////////////////////////////////////////////////////////////////////////

void LeastSquareP1PolyRec3D::configure ( Config::ConfigArgs& args )
{
  FVMCC_PolyRec::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> > LeastSquareP1PolyRec3D::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result = FVMCC_PolyRec::needsSockets();

  // Add the needed DataSocketSinks
  result.push_back(&socket_stencil);
  result.push_back(&socket_weights);
  result.push_back(&socket_uX);
  result.push_back(&socket_uY);
  result.push_back(&socket_uZ);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void LeastSquareP1PolyRec3D::computeGradients()
{
  prepareReconstruction();

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<vector<State*> > stencil = socket_stencil.getDataHandle();
  DataHandle<CFreal> weights = socket_weights.getDataHandle();
  DataHandle<CFreal> uX = socket_uX.getDataHandle();
  DataHandle<CFreal> uY = socket_uY.getDataHandle();
  DataHandle<CFreal> uZ = socket_uZ.getDataHandle();

  const CFuint nbStates = states.size();
  const CFuint nbEquations = PhysicalModelStack::getActive()->getNbEq();

  for(CFuint iVar = 0; iVar < nbEquations; ++iVar) {
    _lf1 = 0.0;
    _lf2 = 0.0;
    _lf3 = 0.0;
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
	  const RealVector& nodeFirst = first->getCoordinates();
	  const RealVector& nodeLast = last->getCoordinates();

	  const CFreal dx = weights[iEdge]*(nodeLast[0]
					    - nodeFirst[0]);
	  const CFreal dy = weights[iEdge]*(nodeLast[1]
					    - nodeFirst[1]);
	  const CFreal dz = weights[iEdge]*(nodeLast[2]
					    - nodeFirst[2]);

	  const CFreal du = weights[iEdge]*((*last)[iVar] - (*first)[iVar]);

	  CFLogDebugMin( "du = " << du << "\n");

	    _lf1[firstID] += dx*du;
	    _lf2[firstID] += dy*du;
	    _lf3[firstID] += dz*du;
	    
	    if (!last->isGhost()) {
		_lf1[lastID]  += dx*du; 
		_lf2[lastID]  += dy*du;
		_lf3[lastID]  += dz*du;
	    }
	  ++iEdge;
	}
      }
    }

    for(CFuint iState = 0; iState < nbStates; ++iState) {

      const CFreal det = _l11[iState]*_l22[iState]*_l33[iState]
	- _l11[iState]*_l23[iState]*_l23[iState]
	- _l12[iState]*_l12[iState]*_l33[iState]
	+ _l12[iState]*_l13[iState]*_l23[iState]
	+ _l13[iState]*_l12[iState]*_l23[iState]
	- _l13[iState]*_l13[iState]*_l22[iState];

      if (!(std::abs(det) > MathTools::MathConsts::CFrealEps())) {
	CFout << "Det is zero."<<"\n";
      }

      const CFreal linv11 = _l22[iState]*_l33[iState] - _l23[iState]*_l23[iState];
      const CFreal linv22 = _l11[iState]*_l33[iState] - _l13[iState]*_l13[iState];
      const CFreal linv33 = _l11[iState]*_l22[iState] - _l12[iState]*_l12[iState];
      const CFreal linv12 = -(_l12[iState]*_l33[iState] - _l13[iState]*_l23[iState]);
      const CFreal linv13 = _l12[iState]*_l23[iState] - _l13[iState]*_l22[iState];
      const CFreal linv23 = -(_l11[iState]*_l23[iState] - _l13[iState]*_l12[iState]);

      // A cure to the singularites in calculating the determinant
      if (!MathChecks::isZero(det)) {
	uX(iState,iVar,nbEquations) = (linv11*_lf1[iState] +
				       linv12*_lf2[iState] +
				       linv13*_lf3[iState])/det;
	uY(iState,iVar,nbEquations) = (linv12*_lf1[iState] +
				       linv22*_lf2[iState] +
				       linv23*_lf3[iState])/det;
	uZ(iState,iVar,nbEquations) = (linv13*_lf1[iState] +
				       linv23*_lf2[iState] +
				       linv33*_lf3[iState])/det;
	
	// if (std::abs(_uX(iState,iVar,nbEquations)) > 0.) CFout << "ux = " << _uX(iState,iVar,nbEquations) << "\n";
	// if (std::abs(_uY(iState,iVar,nbEquations)) > 0.) CFout << "uy = " << _uY(iState,iVar,nbEquations) << "\n";
	// if (std::abs(_uZ(iState,iVar,nbEquations)) > 0.) CFout << "uz = " << _uZ(iState,iVar,nbEquations) << "\n";
      }
      else {
	uX(iState,iVar,nbEquations) = 0.0;
	uY(iState,iVar,nbEquations) = 0.0;
	uZ(iState,iVar,nbEquations) = 0.0;
      }

      CFLogDebugMed( "det = " << det
		     << ", l11 = " << _l11[iState]
		     << ", l12 = " << _l12[iState]
		     << ", l13 = " << _l13[iState]
		     << ", l22 = " << _l22[iState]
		     << ", l23 = " << _l23[iState]
		     << ", l33 = " << _l33[iState]
		     <<", lf1 = " << _lf1[iState]
		     <<", lf2 = " << _lf2[iState]
		     <<", lf3 = " << _lf3[iState]
		     << ", uX =" << uX(iState,iVar,nbEquations)
		     << ", uY =" << uY(iState,iVar,nbEquations)
		     << ", uZ =" << uZ(iState,iVar,nbEquations)
		     << "\n");
    }
  }
  
  
  ///// TEST
  // CFuint idxr = 0;
  // for (CFuint cellID = 0; cellID < 5000; ++cellID) {
  //   for (CFuint i = 0; i < 9; ++i, ++idxr) {
  //     cout << "cellID["<< cellID << "], "<< i << " => UX (";
  //     cout.precision(12); cout << uX[idxr] << ", " << uY[idxr] << ", " << uZ[idxr] << ")\n";
  //   }
  // }
  // abort();
  ////
  
}

//////////////////////////////////////////////////////////////////////////////

void LeastSquareP1PolyRec3D::extrapolateImpl(GeometricEntity* const face)
{
  FVMCC_PolyRec::baseExtrapolateImpl(face);
  
  // first compute the limiter
  computeFaceLimiter(face);
  
  DataHandle<CFreal> uX = socket_uX.getDataHandle();
  DataHandle<CFreal> uY = socket_uY.getDataHandle();
  DataHandle<CFreal> uZ = socket_uZ.getDataHandle();
  DataHandle<CFreal> newLimiter = socket_limiter.getDataHandle();

  // please note that you work with references !!!
  // don't forget the "&" !!!
  const vector<Node*>& coord = getCoord();
  vector<State*>& valuesLR  = getExtrapolatedValues();
  const State *const state = face->getState(LEFT);
  const CFuint stateID = state->getLocalID();
  const RealVector& stateCoord = state->getCoordinates();
  const State *const neighState = face->getState(RIGHT);
  const RealVector& neighStateCoord = neighState->getCoordinates();
  const CFuint neighStateIDgrad = (!isBoundaryFace()) ? neighState->getLocalID() : stateID;
  
  // number of quadrature points associated with this face
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint startID = stateID*nbEqs;
  const CFuint neighStartID = neighStateIDgrad*nbEqs;
  
  const CFreal xq = (*coord[0])[XX];
  const CFreal yq = (*coord[0])[YY];
  const CFreal zq = (*coord[0])[ZZ];
  
  if (_vars.size() > 1) {
    _vFunction.evaluate((*coord[0]),_gradientCoeff);
  }
  
  for(CFuint iVar = 0; iVar < nbEqs; ++iVar) {
    // AL: gradientCoeff prevails on the interactive gradientFactor
    if (_vars.size() > 1) {_gradientFactor = _gradientCoeff[iVar];}
    
    const CFreal gradientCoeffState = (!(*_zeroGradient)[iVar]) ?
      _gradientFactor*(uX(stateID,iVar,nbEqs)*(xq - stateCoord[XX]) + 
		       uY(stateID,iVar,nbEqs)*(yq - stateCoord[YY]) + 
		       uZ(stateID,iVar,nbEqs)*(zq - stateCoord[ZZ])) : 0.0;
    
    const CFreal gradientCoeffNeighbor = (!(*_zeroGradient)[iVar]) ?
      _gradientFactor*(uX(neighStateIDgrad,iVar,nbEqs)*(xq - neighStateCoord[XX]) + 
		       uY(neighStateIDgrad,iVar,nbEqs)*(yq - neighStateCoord[YY]) + 
		       uZ(neighStateIDgrad,iVar,nbEqs)*(zq - neighStateCoord[ZZ])) : 0.0;
    
    // left reconstructed value (the one inside the current cell)
    const CFreal lValue = (_stopLimiting == 0) ? newLimiter[startID + iVar] : 1.;
    (*valuesLR[0])[iVar] = (*state)[iVar] + lValue*gradientCoeffState;
    
    // right reconstructed value (the one inside the neighbor cell)
    const CFreal rValue = (_stopLimiting == 0) ? newLimiter[neighStartID + iVar] : 1.;
    (*valuesLR[1])[iVar] = (*neighState)[iVar] + rValue*gradientCoeffNeighbor;
        
    getBackupValues(0)[iVar] =  (*valuesLR[0])[iVar];
    getBackupValues(1)[iVar] =  (*valuesLR[1])[iVar];
  }
}

//////////////////////////////////////////////////////////////////////////////

void LeastSquareP1PolyRec3D::extrapolateImpl(GeometricEntity* const face,
					     CFuint iVar, CFuint leftOrRight)
{
  DataHandle<CFreal> uX = socket_uX.getDataHandle();
  DataHandle<CFreal> uY = socket_uY.getDataHandle();
  DataHandle<CFreal> uZ = socket_uZ.getDataHandle();
  DataHandle<CFreal> newLimiter = socket_limiter.getDataHandle();

  const State *const state = face->getState(leftOrRight);
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  CFuint stateID = 0;
  if (leftOrRight == LEFT) {
    stateID = state->getLocalID();
  }
  else {
    stateID = (!isBoundaryFace()) ? state->getLocalID() :
      face->getState(LEFT)->getLocalID();
  }
  const RealVector& stateCoord = state->getCoordinates();
  
  // number of quadrature points associated with this face
  const vector<Node*>& coord = getCoord();
  const CFuint startID = stateID*PhysicalModelStack::getActive()->getNbEq();

  const CFreal xq = (*coord[0])[XX];
  const CFreal yq = (*coord[0])[YY];
  const CFreal zq = (*coord[0])[ZZ];
  copyBackupValues();
   
  if (_vars.size() > 1) {
    _vFunction.evaluate(iVar,(*coord[0]), _gradientFactor);
  }
  
  // AL: gradientCoeff prevails on the interactive gradientFactor
  const CFreal gradientCoeffState = (!(*_zeroGradient)[iVar]) ?
    _gradientFactor*(uX(stateID,iVar,nbEqs)*(xq - stateCoord[XX]) + 
		     uY(stateID,iVar,nbEqs)*(yq - stateCoord[YY]) + 
		     uZ(stateID,iVar,nbEqs)*(zq - stateCoord[ZZ])) : 0.0;
  
  // reconstructed value
  const CFreal lValue = (_stopLimiting == 0) ? newLimiter[startID + iVar] : 1.;
  getValues(leftOrRight)[iVar] = (*state)[iVar] + lValue*gradientCoeffState;
}

//////////////////////////////////////////////////////////////////////////////

void LeastSquareP1PolyRec3D::updateWeights()
{

  FVMCC_PolyRec::updateWeights();

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<vector<State*> > stencil = socket_stencil.getDataHandle();
  DataHandle<CFreal> weights = socket_weights.getDataHandle();

  const CFuint nbStates = states.size();

  _l11 = 0.0;
  _l12 = 0.0;
  _l13 = 0.0;
  _l22 = 0.0;
  _l23 = 0.0;
  _l33 = 0.0;

  _lf1 = 0.0;
  _lf2 = 0.0;
  _lf3 = 0.0;

  // weight coefficients are calculated
  // weight coefficients are calculated
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
	const RealVector& nodeFirst = first->getCoordinates();
	const RealVector& nodeLast = last->getCoordinates();
	const CFreal deltaR = MathFunctions::getDistance(first->getCoordinates(),last->getCoordinates());
	
	weights[iEdge] = 1.0/deltaR;

	
	// weights always != 0
	const CFreal dx = weights[iEdge]*(nodeLast[0]
					  - nodeFirst[0]);
	const CFreal dy = weights[iEdge]*(nodeLast[1]
					  - nodeFirst[1]);
	const CFreal dz = weights[iEdge]*(nodeLast[2]
					  - nodeFirst[2]);
	
	CFLogDebugMax( "weights = " << weights[iEdge]
		       << "dx = " << dx
		       << "dy = " << dy
		       << "dz = " << dz << "\n");
	
	_l11[firstID] += dx*dx;
	_l12[firstID] += dx*dy;
	_l13[firstID] += dx*dz;
	_l22[firstID] += dy*dy;
	_l23[firstID] += dy*dz;
	_l33[firstID] += dz*dz;
	
	if (!last->isGhost()) {
	  _l11[lastID] += dx*dx;
	  _l12[lastID] += dx*dy;
	  _l13[lastID] += dx*dz;
	  _l22[lastID] += dy*dy;
	  _l23[lastID] += dy*dz;
	  _l33[lastID] += dz*dz;
	}
	++iEdge;
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void LeastSquareP1PolyRec3D::setup()
{
  FVMCC_PolyRec::setup();

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<vector<State*> > stencil = socket_stencil.getDataHandle();
  DataHandle<CFreal> weights = socket_weights.getDataHandle();

  const CFuint nbStates = states.size();

  _l11.resize(nbStates);
  _l12.resize(nbStates);
  _l13.resize(nbStates);
  _l22.resize(nbStates);
  _l23.resize(nbStates);
  _l33.resize(nbStates);

  _lf1.resize(nbStates);
  _lf2.resize(nbStates);
  _lf3.resize(nbStates);

  _l11 = 0.0;
  _l12 = 0.0;
  _l13 = 0.0;
  _l22 = 0.0;
  _l23 = 0.0;
  _l33 = 0.0;

  _lf1 = 0.0;
  _lf2 = 0.0;
  _lf3 = 0.0;

  // weight coefficients are calculated
  // weight coefficients are calculated
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
	const RealVector& nodeFirst = first->getCoordinates();
	const RealVector& nodeLast = last->getCoordinates();
	const CFreal deltaR = MathFunctions::getDistance(first->getCoordinates(),last->getCoordinates());
	
	
	weights[iEdge] = 1.0/deltaR;
	//  weights[iEdge] = 1.0;
	
	// weights always != 0
	const CFreal dx = weights[iEdge]*(nodeLast[0] - nodeFirst[0]);
	const CFreal dy = weights[iEdge]*(nodeLast[1] - nodeFirst[1]);
	const CFreal dz = weights[iEdge]*(nodeLast[2] - nodeFirst[2]);
	
	CFLogDebugMax( "weights = " << weights[iEdge]
		       << "dx = " << dx
		       << "dy = " << dy
		       << "dz = " << dz << "\n");
	
	_l11[firstID] += dx*dx;
	_l12[firstID] += dx*dy;
	_l13[firstID] += dx*dz;
	_l22[firstID] += dy*dy;
	_l23[firstID] += dy*dz;
	_l33[firstID] += dz*dz;
	
	if (!last->isGhost()) {
	  _l11[lastID] += dx*dx;
	  _l12[lastID] += dx*dy;
	  _l13[lastID] += dx*dz;
	  _l22[lastID] += dy*dy;
	  _l23[lastID] += dy*dz;
	  _l33[lastID] += dz*dz;
	}
	++iEdge;
      }
    }
  }
}
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

 } // namespace COOLFluiD
