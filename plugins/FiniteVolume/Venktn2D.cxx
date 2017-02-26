#include "Venktn2D.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Framework/MeshData.hh"
#include "Common/CFLog.hh"
#include "MathTools/MathFunctions.hh"
#include "FiniteVolume/FiniteVolume.hh"

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

MethodStrategyProvider<Venktn2D,CellCenterFVMData,Limiter<CellCenterFVMData>,FiniteVolumeModule>
venktn2DProvider("Venktn2D");

//////////////////////////////////////////////////////////////////////////////

void Venktn2D::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< vector<CFreal> >
    ("magnitudeValues","State values with order of magnitude of solution.");
  options.addConfigOption< CFreal >
    ("coeffEps","Coefficient for the epsilon.");
  options.addConfigOption< CFreal >
    ("length","Characteristic solution length in the smooth flow region.");
  options.addConfigOption< bool >
    ("isMFMHD","Fix for smooth flow region.");
}

//////////////////////////////////////////////////////////////////////////////

Venktn2D::Venktn2D(const std::string& name) :
  Limiter<CellCenterFVMData>(name),
  socket_stencil("stencil"),
  socket_uX("uX"),
  socket_uY("uY"),
  _deltaMin(0.)
{
  addConfigOptionsTo(this);
  
  _magnitudeValues = vector<CFreal>();
  setParameter("magnitudeValues",&_magnitudeValues);
  
  _coeffEps = 1.0;
  setParameter("coeffEps",&_coeffEps);
  
  // default negative value for convenience  
  _length = -1.0;
  setParameter("length",&_length);

  _isMFMHD = false;
  setParameter("isMFMHD",&_isMFMHD);
}

//////////////////////////////////////////////////////////////////////////////

Venktn2D::~Venktn2D()
{
}

//////////////////////////////////////////////////////////////////////////////

void Venktn2D::configure ( Config::ConfigArgs& args )
{
  Limiter<CellCenterFVMData>::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void Venktn2D::setup()
{
  Limiter<CellCenterFVMData>::setup();
  
  if (_length <= 0.) {
    CFLog(VERBOSE, "WARNING: Venktn2D::configure() => using default reference length\n");
    _length = PhysicalModelStack::getActive()->getImplementor()->getRefLength();
  }
  
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const RealVector& refValues =
    PhysicalModelStack::getActive()->getImplementor()->getRefStateValues();
  
  if (_magnitudeValues.size() < nbEqs) {
    CFLog(VERBOSE, "WARNING: Venktn2D::configure() => using default reference state values\n");
    _magnitudeValues.resize(nbEqs);
    for (CFuint i = 0; i< nbEqs; ++i) {
      _magnitudeValues[i] = refValues[i]; 
    }
  }
}
      
//////////////////////////////////////////////////////////////////////////////

void Venktn2D::limit(const vector<vector<Node*> >& coord,
                     GeometricEntity* const cell,
                     CFreal* limiterValue)
{
  const GeomEntList *const faces = cell->getNeighborGeos();
  const CFuint nbFaces = faces->size();
  const State *const state = cell->getState(0);
  const CFuint stateID = state->getLocalID();
  const CFuint nbEquations = PhysicalModelStack::getActive()->getNbEq();
  DataHandle< vector<State*> > stencil = socket_stencil.getDataHandle();
  
  for(CFuint iVar = 0; iVar < nbEquations; ++iVar) {
    //std::cout<<"iVar = "<< iVar <<"\n";
    CFreal min0 = (*state)[iVar];
    CFreal max0 = (*state)[iVar];
    CFreal avgDistance = 0.0;
    
    if (this->m_useFullStencil) {
      CFuint nbNeighbors = 0;
      if (!this->m_useNodalExtrapolationStencil) {
	const vector<State*>& s = stencil[stateID];
	nbNeighbors = s.size();
	for (CFuint is = 0; is < nbNeighbors; ++is) {
	  State *const currState = s[is];
	  if(state != currState) {
	    min0 = min(min0,(*currState)[iVar]);
	    max0 = max(max0,(*currState)[iVar]);
	    avgDistance += MathFunctions::getDistance
	      (state->getCoordinates(), currState->getCoordinates());
	  }
	}
      }
      else {
	// old method kept for backward compatibility: it may escludes some ghost states at corners
	const CFuint nbNodesInCell = cell->nbNodes();
	for (CFuint iNode = 0; iNode < nbNodesInCell; ++iNode) {
	  const CFuint nodeID = cell->getNode(iNode)->getLocalID();
	  const vector<State*>& s = getMethodData().
	    getNodalStatesExtrapolator()->getNodalStateNeighbors(nodeID);
	  for (CFuint is = 0; is < s.size(); ++is) {
	    State *const currState = s[is];
	    if(state != currState) {
	      min0 = min(min0,(*currState)[iVar]);
	      max0 = max(max0,(*currState)[iVar]);
	      avgDistance += MathFunctions::getDistance
		(state->getCoordinates(), currState->getCoordinates());
	      nbNeighbors++;
	    }
	  }
	}
      }
      
      // only uniquely counted neighbors should be considered ...
      avgDistance /= nbNeighbors;
    }
    else {
      for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
	const GeometricEntity *const neighborFace = (*faces)[iFace];
        State *const neighborState = (state == neighborFace->getState(0)) ?
	  neighborFace->getState(1) : neighborFace->getState(0);
	
	min0 = min(min0,(*neighborState)[iVar]);
	max0 = max(max0,(*neighborState)[iVar]);
	
	// the average distance between the centroids of the cell and
	// each of its face neighbors is calculated
	avgDistance += MathFunctions::getDistance
	  (state->getCoordinates(), neighborState->getCoordinates());
      }
      
      avgDistance /= nbFaces;
    }
    
    const CFreal deltaPlusMax = max0 - (*state)[iVar];
    const CFreal deltaPlusMin = min0 - (*state)[iVar];
    
    CFreal psi = 1.0;
    CFreal psimin = 1.1;
    if(_isMFMHD){ // IMPLEMENTATION FROM VENKATAKRISHNAN PAPER 
      const CFreal epsilon2 = pow(_coeffEps*avgDistance/_length, 3.0); //_magnitudeValues[iVar]*_magnitudeValues[iVar]*
    
      for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
        // number of quadrature points associated with this face
        const CFuint nbPoints = 1;
      
        for (CFuint ip = 0; ip < nbPoints; ++ip) {
	  computeDeltaMin((*coord[iFace][ip]), *state, iVar);
	
          if (_deltaMin > 0.0) {
	    const CFreal dMinStar = MathTools::MathFunctions::sign(_deltaMin)*(std::abs(_deltaMin)/_magnitudeValues[iVar] + 1e-24);
	    const CFreal dPlus = deltaPlusMax/_magnitudeValues[iVar];
	    const CFreal dPlus2   = dPlus*dPlus;
	    const CFreal dPlusMin = dPlus*dMinStar;
	    const CFreal Num = (dPlus2 + epsilon2)*dMinStar + 2*dMinStar*dMinStar*dPlus;
	    const CFreal Den = (dPlus2 + 2*dMinStar*dMinStar + dPlusMin + epsilon2);
	    psi = 1./dMinStar*(Num/Den);
	  }
          if (_deltaMin < 0.0) {
	    const CFreal dMinStar = MathTools::MathFunctions::sign(_deltaMin)*(std::abs(_deltaMin)/_magnitudeValues[iVar] + 1e-24);
            const CFreal dPlus = deltaPlusMin/_magnitudeValues[iVar];
            const CFreal dPlus2   = dPlus*dPlus;
            const CFreal dPlusMin = dPlus*dMinStar;
            const CFreal Num = (dPlus2 + epsilon2)*dMinStar + 2*dMinStar*dMinStar*dPlus;
            const CFreal Den = (dPlus2 + 2*dMinStar*dMinStar + dPlusMin + epsilon2);
            psi = 1./dMinStar*(Num/Den);
	  }
          psimin = min(psi, psimin);
        }
      }
    }
    else{ //OLD IMPLEMENTATION
      const CFreal epsilon = _coeffEps*_magnitudeValues[iVar]*_magnitudeValues[iVar]*
        pow(avgDistance/_length, 3.0);

      for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
        const CFuint nbPoints = 1;

        for (CFuint ip = 0; ip < nbPoints; ++ip) {
          computeDeltaMin((*coord[iFace][ip]), *state, iVar);
 
          if (_deltaMin > 0.0) {
            const CFreal dPlusMax2   = deltaPlusMax*deltaPlusMax;
            const CFreal dPlusMaxMin = deltaPlusMax*_deltaMin;
            psi = (dPlusMax2 + 2.0*dPlusMaxMin + epsilon)/
              (dPlusMax2 + dPlusMaxMin + 2.0*_deltaMin*_deltaMin + epsilon);
          }
          if (_deltaMin < 0.0) {
            const CFreal y = deltaPlusMin/_deltaMin;
            psi = (y*y + 2*y)/(y*y + y + 2.);
          }
          psimin = min(psi, psimin);
        }
      }
    }
    limiterValue[iVar] = psimin;
    
    const CFreal maxAllowableLimiterFunctionValue = 1.094;
    
    if (limiterValue[iVar] > maxAllowableLimiterFunctionValue) {
      CFout << "wrong limiterValue = " << limiterValue[iVar] << "\n";
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

