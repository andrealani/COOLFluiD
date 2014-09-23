#include "NavierStokesKOmega3DSourceTerm.hh"
#include "KOmega/NavierStokesKOmegaVarSetTypes.hh"
#include "Common/CFLog.hh"
#include "Framework/GeometricEntity.hh"
#include "Framework/MeshData.hh"
#include "FiniteVolumeTurb/FiniteVolumeKOmega.hh"
#include "Framework/SubSystemStatus.hh"

#include "Framework/MethodStrategyProvider.hh"
#include "FiniteVolume/CellCenterFVMData.hh"

#include "FiniteVolume/DerivativeComputer.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::Physics::KOmega;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<NavierStokesKOmega3DSourceTerm,
		       CellCenterFVMData,
		       ComputeSourceTerm<CellCenterFVMData>,
		       FiniteVolumeKOmegaModule>
NavierStokesKOmega3DSourceTermFVMCCProvider("NavierStokesKOmega3DSourceTerm");

//////////////////////////////////////////////////////////////////////////////

NavierStokesKOmega3DSourceTerm::NavierStokesKOmega3DSourceTerm(const std::string& name) :
  ComputeSourceTermFVMCC(name),
  _varSet(CFNULL),
  _diffVarSet(CFNULL),
  _temp(),
  _avState(),
  _physicalData(),
  _nstates(CFNULL),
  _wallDistance(CFNULL),
  _values(),
  _states(),
  _unperturbedPositivePart(),
  _unperturbedNegativePart(),
  _gradients(),
  _averageFaceValue()
{
}

//////////////////////////////////////////////////////////////////////////////

NavierStokesKOmega3DSourceTerm::~NavierStokesKOmega3DSourceTerm()
{
  for(CFuint iGrad = 0; iGrad < _gradients.size(); iGrad++)
  {
    deletePtr(_gradients[iGrad]);
  }

}

//////////////////////////////////////////////////////////////////////////////

void NavierStokesKOmega3DSourceTerm::setup()
{
  CFAUTOTRACE;

  ComputeSourceTermFVMCC::setup();
  
  _varSet = getMethodData().getUpdateVar().d_castTo<EulerVarSet>();
  _diffVarSet = getMethodData().getDiffusiveVar().d_castTo<NavierStokes3DKOmega>();
  
  _temp.resize(PhysicalModelStack::getActive()->getNbEq());
  _avState.resize(PhysicalModelStack::getActive()->getNbEq());
  
  cf_assert(_varSet.isNotNull());
  _varSet->getModel()->resizePhysicalData(_physicalData);

  _nstates = _sockets.getSocketSink<RealVector>("nstates")->getDataHandle();

  _wallDistance = _sockets.getSocketSink<CFreal>("wallDistance")->getDataHandle();

  SafePtr<DerivativeComputer> derComput =
    this->getMethodData().getDerivativeComputer();
  const CFuint nbNodesInControlVolume =
    derComput->getMaxNbVerticesInControlVolume();

  _values.resize(PhysicalModelStack::getActive()->getNbEq(), nbNodesInControlVolume);
  _states.reserve(PhysicalModelStack::getActive()->getNbEq());

  const CFuint nbScalarVars = _varSet->getModel()->getNbScalarVars(0);
  _unperturbedPositivePart.resize(nbScalarVars);
  _unperturbedNegativePart.resize(nbScalarVars);

  _gradients.resize(PhysicalModelStack::getActive()->getNbEq());
  for(CFuint iGrad = 0; iGrad < _gradients.size(); iGrad++)
  {
    _gradients[iGrad] = new RealVector(DIM_3D);
  }

  _averageFaceValue.resize(PhysicalModelStack::getActive()->getNbEq());
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokesKOmega3DSourceTerm::computeSource(Framework::GeometricEntity *const element,
						   RealVector& source,
						   RealMatrix& jacobian)
{
  DataHandle<CFreal> normals = this->socket_normals.getDataHandle();
  DataHandle<CFint> isOutward = this->socket_isOutward.getDataHandle();
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();

  cf_assert(_varSet.isNotNull());

  // Set the physical data for the cell considered
  State *const currState = element->getState(0);
  _varSet->computePhysicalData(*currState, _physicalData);

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

  // fill in the nodal states
  const vector<Node*>* const nodes = element->getNodes();
  const CFuint nbNodesInElem = nodes->size();
  _states.clear();
  for (CFuint i = 0; i < nbNodesInElem; ++i) {
    _states.push_back(&_nstates[(*nodes)[i]->getLocalID()]);
  }

  //From now on, we will use the gradient vars
  _diffVarSet->setGradientVars(_states, _values, _states.size());

  const vector<GeometricEntity*>& faces = *element->getNeighborGeos();
//   cf_assert(faces.size() == nbNodesInElem);
  const CFuint elemID = element->getID();

  // compute the gradients by applying Green Gauss in the
  // cell d's
  for(CFuint iGrad = 0; iGrad < _gradients.size(); iGrad++)
  {
    *(_gradients[iGrad]) = 0.0;
  }

  ///@todo move this as a member data
  Common::CFMap<CFuint,CFuint> nodeID2cellNodeIDMap(nbNodesInElem);
  for (CFuint iNode = 0; iNode < nbNodesInElem; ++iNode) {
    const CFuint nodeID = element->getNode(iNode)->getLocalID();
      nodeID2cellNodeIDMap.insert(nodeID, iNode);
  }
  nodeID2cellNodeIDMap.sortKeys();

  const CFuint nbFaces = faces.size();
  for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
    // get the face normal
    const CFuint faceID = faces[iFace]->getID();
    const CFuint startID = faceID*PhysicalModelStack::getActive()->getDim();
    CFreal nx = normals[startID];
    CFreal ny = normals[startID + 1];
    CFreal nz = normals[startID + 2];
    if (static_cast<CFuint>( isOutward[faceID]) != elemID) {
      nx *= -1.;
      ny *= -1.;
      nz *= -1.;
    }
    _averageFaceValue = 0.;
    const CFuint nbNodesInFace = faces[iFace]->nbNodes();
    for (CFuint iNode = 0; iNode < nbNodesInFace; ++iNode) {
      const CFuint nodeID = faces[iFace]->getNode(iNode)->getLocalID();
      const CFuint cellNodeID = nodeID2cellNodeIDMap.find(nodeID);
      for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
	_averageFaceValue[iEq] += _values(iEq, cellNodeID);
      }
    }
    _averageFaceValue /= nbNodesInFace;

    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      (*(_gradients[iEq]))[XX] += nx*_averageFaceValue[iEq];
      (*(_gradients[iEq]))[YY] += ny*_averageFaceValue[iEq];
      (*(_gradients[iEq]))[ZZ] += nz*_averageFaceValue[iEq];
    }
  }

  for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
    *(_gradients[iEq]) /= volumes[elemID];
  }

// compute PUVTKOmega by averaging the nodes
// NO!!! If we do it like that we nearly certainly
// get negative values!!!
// So we just take the state value
  const CFuint iK = _varSet->getModel()->getFirstScalarVar(0);

  _avState[0] = _physicalData[EulerTerm::P];
  _avState[1] = _physicalData[EulerTerm::VX];
  _avState[2] = _physicalData[EulerTerm::VY];
  _avState[3] = _physicalData[EulerTerm::VZ];
  _avState[4] = _physicalData[EulerTerm::T];
  _avState[5] = _physicalData[iK];
  _avState[6] = _physicalData[iK+1];

  CFreal avK = _physicalData[iK];
  CFreal avOmega = _physicalData[iK+1];

  const CFreal rho = _diffVarSet->getDensity(_avState);

  ///Get the wall distance
  const CFreal avDist = _wallDistance[currState->getLocalID()];

  ///Set the wall distance before computing the turbulent viscosity
  _diffVarSet->setWallDistance(avDist);
  CFreal mut = _diffVarSet->getTurbDynViscosityFromGradientVars(_avState, _gradients);
  _diffVarSet->computeBlendingCoefFromGradientVars(_avState, *(_gradients[5]), *(_gradients[6]));

  const CFreal blendingCoefF1 = _diffVarSet->getBlendingCoefficientF1();
  const CFreal sigmaOmega2 = _diffVarSet->getSigmaOmega2();

  ///Compute Reynolds stress tensor
  const CFreal twoThirdRhoK = (2./3.)*(avK * rho);
  const CFreal twoThirdDivV = (2./3.)*((*(_gradients[1]))[XX] + (*(_gradients[2]))[YY] + (*(_gradients[3]))[ZZ]);

  const CFreal coeffTauMu = _diffVarSet->getModel().getCoeffTau();
  CFreal tauXX = coeffTauMu * (mut * (2.*(*(_gradients[1]))[XX] - twoThirdDivV)) - twoThirdRhoK;
  CFreal tauYY = coeffTauMu * (mut * (2.*(*(_gradients[2]))[YY] - twoThirdDivV)) - twoThirdRhoK;
  CFreal tauZZ = coeffTauMu * (mut * (2.*(*(_gradients[3]))[ZZ] - twoThirdDivV)) - twoThirdRhoK;

  CFreal tauXY = coeffTauMu * (mut * ((*(_gradients[1]))[YY] + (*(_gradients[2]))[XX]));
  CFreal tauYX = tauXY;
  CFreal tauYZ = coeffTauMu * (mut * ((*(_gradients[2]))[ZZ] + (*(_gradients[3]))[YY]));
  CFreal tauZY = tauYZ;
  CFreal tauXZ = coeffTauMu * (mut * ((*(_gradients[3]))[XX] + (*(_gradients[1]))[ZZ]));
  CFreal tauZX = tauXZ;

  ///Production term: k
  CFreal prodTerm_k = tauXX*(*(_gradients[1]))[XX] + tauXY*(*(_gradients[1]))[YY] + tauXZ*(*(_gradients[1]))[ZZ] +
                      tauYX*(*(_gradients[2]))[XX] + tauYY*(*(_gradients[2]))[YY] + tauYZ*(*(_gradients[2]))[ZZ] +
                      tauZX*(*(_gradients[3]))[XX] + tauZY*(*(_gradients[3]))[YY] + tauZZ*(*(_gradients[3]))[ZZ];

  ///Production term: Omega
  CFreal prodTerm_Omega = (_diffVarSet->getGammaCoef()*rho/mut) *
                          (tauXX*(*(_gradients[1]))[XX] + tauXY*(*(_gradients[1]))[YY] + tauXZ *(*(_gradients[1]))[ZZ] +
                           tauYX*(*(_gradients[2]))[XX] + tauYY*(*(_gradients[2]))[YY] + tauYZ *(*(_gradients[2]))[ZZ] +
                           tauZX*(*(_gradients[3]))[XX] + tauZY*(*(_gradients[3]))[YY] + tauZZ *(*(_gradients[3]))[ZZ]);

  ///This is used in (BSL,SST), not for normal kOmega
  const CFreal overOmega = 1./avOmega;
  prodTerm_Omega += (1. - blendingCoefF1) * 2. * rho * overOmega * sigmaOmega2
                    * ((*(_gradients[5]))[XX]*(*(_gradients[6]))[XX] + (*(_gradients[5]))[YY]*(*(_gradients[6]))[YY] + (*(_gradients[5]))[ZZ]*(*(_gradients[6]))[ZZ]);

  ///Destruction term: k
  CFreal destructionTerm_k =(-1.) * rho * avOmega * avK * _diffVarSet->getBetaStar(_avState);

  ///Destruction term: Omega
  CFreal destructionTerm_Omega = (-1.) * rho * avOmega * avOmega * _diffVarSet->getBeta(_avState);

  //Make sure negative values dont propagate...
  destructionTerm_k     = min(0., destructionTerm_k);
  destructionTerm_Omega = min(0., destructionTerm_Omega);
  prodTerm_k            = max(0., prodTerm_k);
  prodTerm_Omega        = max(0., prodTerm_Omega);

  // This trick was used by W. Dieudonne in Euphoria
  prodTerm_k = min(20.*fabs(destructionTerm_k), prodTerm_k);
  prodTerm_Omega = min(20.*fabs(destructionTerm_Omega), prodTerm_Omega);

  ///Computation of the source term
  source[0] = 0.0;
  source[1] = 0.0;
  source[2] = 0.0;
  source[3] = 0.0;
  source[4] = 0.0;

  //What we do with the source term depends if
  //we are computing the jacobian or not
  const bool isPerturb = this->getMethodData().isPerturb();
  const CFuint iPerturbVar = this->getMethodData().iPerturbVar();
  
  if(isPerturb)
  {
    /// Compute the jacobian contribution
    // only perturb the negative part of the source term
    if(iPerturbVar == 5)
    {
      source[5]  = destructionTerm_k;
      source[5] += _unperturbedPositivePart[0];
    }
    else
    {
      source[5]  = _unperturbedNegativePart[0];
      source[5] += _unperturbedPositivePart[0];
    }

    if(iPerturbVar == 6)
    {
      source[6]  = destructionTerm_Omega;
      source[6] += _unperturbedPositivePart[1];
    }
    else
    {
      source[6]  = _unperturbedNegativePart[1];
      source[6] += _unperturbedPositivePart[1];
    }
  }
  else
  {
    /// Compute the rhs contribution
    // and Store the unperturbed source terms
    source[5] = prodTerm_k;
    source[5] += destructionTerm_k;
    _unperturbedPositivePart[0] = prodTerm_k;
    _unperturbedNegativePart[0] = destructionTerm_k;

    source[6] = prodTerm_Omega;
    source[6] += destructionTerm_Omega;
    _unperturbedPositivePart[1] = prodTerm_Omega;
    _unperturbedNegativePart[1] = destructionTerm_Omega;
  }

  ///Finally multiply by the cell volume
  source *= volumes[elemID];

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
