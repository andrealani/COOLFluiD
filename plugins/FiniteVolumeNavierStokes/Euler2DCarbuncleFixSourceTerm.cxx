#include "FiniteVolumeNavierStokes/Euler2DCarbuncleFixSourceTerm.hh"
#include "NavierStokes/Euler2DVarSet.hh"
#include "Common/CFLog.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/GeometricEntity.hh"
#include "FiniteVolumeNavierStokes/FiniteVolumeNavierStokes.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "FiniteVolume/CellCenterFVMData.hh"
#include "FiniteVolume/DerivativeComputer.hh"
#include "NavierStokes/NavierStokes2DVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<Euler2DCarbuncleFixSourceTerm,
		       CellCenterFVMData,
		       ComputeSourceTerm<CellCenterFVMData>,
		       FiniteVolumeNavierStokesModule>
euler2DCarbunleFixSTFVMCCProvider("Euler2DCarbuncleFixST");

//////////////////////////////////////////////////////////////////////////////
      
void Euler2DCarbuncleFixSourceTerm::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal, Config::DynamicOption<> >
    ("DissipEps", "Coefficient to control the dissipation.");
  
  options.addConfigOption< CFuint, Config::DynamicOption<> >
    ("FreezeArtificialViscosityCoeff", "The artificial viscous coefficient mu_s is not recomputed anymore");
		
  options.addConfigOption< CFuint >("Variant", "Selects a particular implementation of the Carbuncle Fix");
}
      
//////////////////////////////////////////////////////////////////////////////
      
Euler2DCarbuncleFixSourceTerm::Euler2DCarbuncleFixSourceTerm(const std::string& name) :
  ComputeSourceTermFVMCC(name), 
  socket_fix_active("fix_active"),
  socket_ArtViscCoeff("ArtViscCoeff"),
  socket_faceAreas("faceAreas"),
  _varSet(CFNULL),
  _diffVarSet(CFNULL),
  _temp(),
  _physicalData(),
  _dX(2),
  _speed(),
  _soundSpeed(),
  _values(),
  _states()
{ 
  addConfigOptionsTo(this);
 
  _dissipEps = 0.1;
  setParameter("DissipEps",&_dissipEps);
  
  _freeze_mu_s = 0;
  setParameter("FreezeArtificialViscosityCoeff",&_freeze_mu_s);
	
  _variantChosen = 0;
  setParameter("Variant", &_variantChosen);
}
      
//////////////////////////////////////////////////////////////////////////////

Euler2DCarbuncleFixSourceTerm::~Euler2DCarbuncleFixSourceTerm()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DCarbuncleFixSourceTerm::setup()
{
  ComputeSourceTermFVMCC::setup();
  
  const CFuint nbCells =
    MeshDataStack::getActive()->getTrs("InnerCells")->getLocalNbGeoEnts();
  
  DataHandle<CFreal> fix_active = socket_fix_active.getDataHandle();
  fix_active.resize(nbCells);
  
  DataHandle<CFreal> artViscCoeff = socket_ArtViscCoeff.getDataHandle();
  artViscCoeff.resize(nbCells);
  
  _temp.resize(PhysicalModelStack::getActive()->getNbEq());
  
  _varSet = getMethodData().getUpdateVar().d_castTo<Euler2DVarSet>();
  _varSet->getModel()->resizePhysicalData(_physicalData);
	
  _diffVarSet = getMethodData().getDiffusiveVar().d_castTo<Physics::NavierStokes::NavierStokesVarSet>();
  
  const CFuint maxNbFacesIn2D = 4;
  _speed.resize(maxNbFacesIn2D);
  _soundSpeed.resize(maxNbFacesIn2D);
  
  const CFuint maxNbNodesIn2DCV = 4;
  _values.resize(PhysicalModelStack::getActive()->getNbEq(), maxNbNodesIn2DCV);
  _states.reserve(maxNbNodesIn2DCV);
  
  // default setting for _uvID array
  if (_uvID.size() != 2) {
    _uvID.resize(2);
    _uvID[0] = 1;
    _uvID[1] = 2;
  }
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DCarbuncleFixSourceTerm::computeSource(Framework::GeometricEntity *const element,
						  RealVector& source,
						  RealMatrix& jacobian)
{
  CFLogDebugMin( "Euler2DCarbuncleFixSourceTerm::computeSource()" << "\n");
  
  // ALWAYS reset the cell to 0
  source = 0.0;
  
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
  DataHandle<CFreal> normals = socket_normals.getDataHandle();
  
  // remove the following
  DataHandle<CFint> isOutward =  this->socket_isOutward.getDataHandle();
  DataHandle<RealVector> nstates = _sockets.getSocketSink<RealVector>("nstates")->getDataHandle();
	
  cf_assert(_varSet.isNotNull());
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  
  // this source term has to be used in combination with the axisymmetric one
  State *const currState = element->getState(0);
  const vector<GeometricEntity*>& faces = *element->getNeighborGeos();
  const CFuint nbFaces = faces.size();
  
  // first check if the cell is in the shock layer
  _varSet->computePhysicalData(*currState, _physicalData);
  
  const CFreal rhoCell  = _physicalData[EulerTerm::RHO];
  const CFreal velCell  = _physicalData[EulerTerm::V];
  const CFreal uCell    = _physicalData[EulerTerm::VX];
  const CFreal vCell    = _physicalData[EulerTerm::VY];
  const CFreal aCell    = _physicalData[EulerTerm::A];
  const CFreal machCell = velCell/aCell;
  const bool isSupersonic = (std::abs(machCell) >= 1.);
  bool isCrossedByShock = false;
  
  _speed = 0;
  _soundSpeed = 0.;
	
	
	
	
	
	
	
	
  
  // detect if the current cell is crossed by a shock:
  // if the currect cell is supersonic and at least one neighbor is subsonic
  // if the currect cell is subsonic and at least one neighbor is supersonic
  // and gradient of mach in the speed direction is < 0 (Mach decreases through shock) 
  const CFuint elemID = element->getID();
  for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
    const GeometricEntity *const currFace = faces[iFace];
    State* const neighbor = (currFace->getState(1) != currState) ?
      currFace->getState(1) : currFace->getState(0);
    
    // this could be expensive
    _varSet->computePhysicalData(*neighbor, _physicalData);
    
    _speed[iFace]      = _physicalData[EulerTerm::V];
    _soundSpeed[iFace] = _physicalData[EulerTerm::A];
    const CFreal mach = _physicalData[EulerTerm::V]/_physicalData[EulerTerm::A];
    const bool isNeighborSuperSonic = (std::abs(mach) >= 1.);
    
    if ((!isCrossedByShock) && (isNeighborSuperSonic != isSupersonic)) {
      isCrossedByShock = true;
    }
  }
  
  bool isCompression = false;
  if (isCrossedByShock) {
    _states.clear();
    const vector<Node*>* const nodes = element->getNodes();
    const CFuint nbNodesInElem = nodes->size();
    for (CFuint i = 0; i < nbNodesInElem; ++i) {
      _states.push_back(&nstates[(*nodes)[i]->getLocalID()]);
    }
    
    setGradientVars(_states, _values, nbNodesInElem);
    cf_assert(faces.size() == nbNodesInElem);
    
    // compute the gradient of Mach number by applying Green Gauss in the cell
    const CFuint machID = 2;
    CFreal dMachDx = 0.0;
    CFreal dMachDy = 0.0;
    for (CFuint i = 0; i < nbNodesInElem; ++i) {
      // get the face normal
      const CFuint faceID = faces[i]->getID();
      const CFuint startID = faceID*dim;
      CFreal nx = normals[startID];
      CFreal ny = normals[startID + 1];
      if (static_cast<CFuint>(isOutward[faceID]) != elemID) {
        nx *= -1.;
        ny *= -1.;
      }
      
      const CFreal machSum = (i < (nbNodesInElem - 1)) ? 	(_values(machID, i) + _values(machID, i+1)) : (_values(machID, i) +_values(machID,0));
      dMachDx += nx*machSum;
      dMachDy += ny*machSum;
    }
    
    // the following extra computation is irrelevant if we are just interested in the sign of the Mach gradient in velocity direction
    // dMachDx *= 0.5/volumes[elemID]; dMachDy *= 0.5/volumes[elemID];
    
    if (dMachDx*uCell + dMachDy*vCell < 0.) {
      isCompression = true;
    }
  }
  
  // compute source term only if the cell is crossed by a compression shock
  if (isCompression) {
    DataHandle<CFreal> fix_active = socket_fix_active.getDataHandle();
    DataHandle<CFreal> artViscCoeff = socket_ArtViscCoeff.getDataHandle();
    static CFuint nbIter = 0;
    if (SubSystemStatusStack::getActive()->getNbIter() > nbIter) {
      nbIter = SubSystemStatusStack::getActive()->getNbIter();
      fix_active = 0.0;
      if (_freeze_mu_s == 0) {
				artViscCoeff = 0.0;
      }
    }
    
    if (!this->getMethodData().isPerturb()) {
      fix_active[elemID] = 1.;
    }
    
    const CFreal cosD = uCell/velCell;
    const CFreal sinD = vCell/velCell;
    const CFreal lambdaMax = velCell + aCell; 
    
    CFreal csiTerm = 0.; 
    CFreal nDotCsi = 0.; 
    CFreal h = 0.;
    for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
      const GeometricEntity *const currFace = faces[iFace];
      const State* const neighbor = (currFace->getState(1) != currState) ? currFace->getState(1) : currFace->getState(0);
      
      // project the L/R velocity gradient onto the eta direction 
      _dX = neighbor->getCoordinates() - currState->getCoordinates();
      const CFreal dl = _dX.norm2();
      const CFreal lDotEta = -_dX[XX]*sinD + _dX[YY]*cosD;
      const CFreal dVdEta = (_speed[iFace] - velCell)/dl*lDotEta;
      
      // get the outward normal
      const CFuint startID = currFace->getID()*dim;
      CFreal nx = normals[startID];
      CFreal ny = normals[startID + 1];
      if (static_cast<CFuint>(isOutward[currFace->getID()]) != elemID) {
				nx *= -1.;
				ny *= -1.;
      }
      
      const CFreal newNDotCsi = std::abs(nx*cosD + ny*sinD);
      if (newNDotCsi > nDotCsi) {
				// choose as characteristic dimension the size of the face which is more aligned with the shock
				// i.e. for which the norm of the dot product between face normal and velocity direction is maximum 
				// const CFreal h = sqrt(volumes[elemID]);
						nDotCsi = newNDotCsi;
				h = std::max(h, socket_faceAreas.getDataHandle()[currFace->getID()]);
      }
      csiTerm += dVdEta*(-nx*sinD + ny*cosD);
    }
    
    // if the mu_s is not frozen and we are not computing jacobian terms, the mu_s must be recalculated
    // otherwise we use the value previously stored 
    if (_freeze_mu_s == 0 && (!this->getMethodData().isPerturb())) {
      artViscCoeff[elemID] = _dissipEps*lambdaMax*h*rhoCell;
    }
    csiTerm *= artViscCoeff[elemID];
    
// 		std::cout << "Epsilon " << _dissipEps << "\n";
		switch (_variantChosen){
			case 0 :{ 
// 				std::cout << "Case 0 \n";
				source[1] += csiTerm*cosD;
				source[2] += csiTerm*sinD;

				break;
			} // end - case 0
			
			case 11 :{ 
// 				std::cout << "Case 11 \n";
				source[1] += csiTerm*cosD;
				source[2] += csiTerm*sinD;
				source[3] += csiTerm*velCell;
				break;
			} // end - case 11
			
			default :{ 
				std::cout << "Variant not implemented \n";
				abort();
				break;
			} // end - default
			
		}// end - switch

  }
  
  // // // // // // // // // // // // // // // // // // // 	
	
// 	// compute the gradients by applying Green Gauss in the
//     // cell volume
//     CFreal dUdX = 0.0;
//     CFreal dVdR = 0.0;
//     
//     const CFuint uID = _uvID[0];
//     const CFuint vID = _uvID[1];
// 		
// 		const CFuint totalNbEqs = 4;
//     
//     CFLog(DEBUG_MED, "Euler2DCarbuncleFixSourceTerm::computeSource() => uID = " << uID << ", vID = " << vID << "\n"); 
//     
//     if (this->m_useGradientLS && this->m_gradientsExist) {
//       const CFuint start = elemID*totalNbEqs;
//       dUdX = this->m_ux[start+uID];
//       dVdR = this->m_uy[start+vID];
//       CFLog(DEBUG_MED, "Euler2DCarbuncleFixSourceTerm::computeSource() => LS gradient in cell [" << 
// 	    elemID << " ] => dUdX = [" << dUdX << "], dVdR = [" << dVdR << "]\n");
//     }
//     else {
//       _states.clear();
//       const vector<Node*>* const nodes = element->getNodes();
//       const CFuint nbNodesInElem = nodes->size();
//       for (CFuint i = 0; i < nbNodesInElem; ++i) {
// 				_states.push_back(&nstates[(*nodes)[i]->getLocalID()]);
//       }
//       
//       const vector<GeometricEntity*>& faces = *element->getNeighborGeos();
//       cf_assert(faces.size() == nbNodesInElem);
//       
//       _diffVarSet->setGradientVars(_states, _values, nbNodesInElem);
// 			for (CFuint i = 0; i < nbNodesInElem; ++i) {
// 				// get the face normal
// 				const CFuint faceID = faces[i]->getID();
// 				const CFuint startID = faceID*PhysicalModelStack::getActive()->getDim();
// 				CFreal nx = normals[startID];
// 				CFreal ny = normals[startID + 1];
// 				if (static_cast<CFuint>(isOutward[faceID]) != elemID) {
// 					nx *= -1.;
// 					ny *= -1.;
// 				}
// 				
// 				if (i < (nbNodesInElem - 1)) {
// 					dUdX += nx*(_values(uID, i) + _values(uID, i+1));
// 					dVdR += ny*(_values(vID, i) + _values(vID, i+1));
// 				}
// 				else {
// 					dUdX += nx*(_values(uID, i) + _values(uID, 0));
// 					dVdR += ny*(_values(vID, i) + _values(vID, 0));
// 				}
//       }
// 		}
// // // // // // // // // // // // // // // // // // // // //       

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
}
      
//////////////////////////////////////////////////////////////////////////////
    
void Euler2DCarbuncleFixSourceTerm::setGradientVars(const vector<RealVector*>& states,
						    RealMatrix& values,
						    CFuint stateSize)
{  
  // We suppose to have u T variables here 
  // gradients of V, T, M
  const CFreal R = 287.046;
  const CFreal gamma = 1.4;
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  
  for (CFuint i = 0; i < stateSize; ++i) {
    const RealVector& state = *states[i];
    CFreal V2 = 0.;
    for (CFuint iDim = 0; iDim < dim; ++iDim) {
      V2 += state[1+iDim]*state[1+iDim];
    }
    values(0,i) = sqrt(V2);                               // velocity magnitude
    values(1,i) = state[dim + 1];                        // temperature
    values(2,i) = values(0,i)/sqrt(gamma*R*values(1,i)); // Mach number 
  }
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<Framework::BaseDataSocketSink> >
Euler2DCarbuncleFixSourceTerm::needsSockets()
{
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result =
    ComputeSourceTermFVMCC::needsSockets();
  result.push_back(&socket_faceAreas);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > 
Euler2DCarbuncleFixSourceTerm::providesSockets()
{
  std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > result = 
    ComputeSourceTermFVMCC::providesSockets();
  result.push_back(&socket_fix_active); 
  result.push_back(&socket_ArtViscCoeff);
  return result;
}
      
//////////////////////////////////////////////////////////////////////////////

  } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

