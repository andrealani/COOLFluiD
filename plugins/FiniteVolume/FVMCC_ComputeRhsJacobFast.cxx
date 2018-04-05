#include "FiniteVolume/FiniteVolume.hh"
#include "FVMCC_ComputeRhsJacobFast.hh"
#include "FiniteVolume/FVMCC_BC.hh"
#include "Framework/CFL.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/BlockAccumulator.hh"
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

MethodCommandProvider<FVMCC_ComputeRhsJacobFast,
		      CellCenterFVMData,
		      FiniteVolumeModule>
fvmcc_computeRhsJacobFast("NumJacobFast");

//////////////////////////////////////////////////////////////////////////////

FVMCC_ComputeRhsJacobFast::FVMCC_ComputeRhsJacobFast
(const std::string& name) :
  FVMCC_ComputeRhsJacob(name)
{
}

//////////////////////////////////////////////////////////////////////////////

FVMCC_ComputeRhsJacobFast::~FVMCC_ComputeRhsJacobFast()
{
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacobFast::computeConvDiffFluxes(CFuint iVar, CFuint iCell)
{  
  // extrapolate (and LIMIT, if the reconstruction is linear or more)
  // the solution in the quadrature points
  _polyRec->extrapolate(_currFace, iVar, iCell);
  
  // compute the physical data for each left and right reconstructed
  // state and in the left and right cell centers
  computePerturbedStatesData(iVar, iCell);

  _pertFlux = 0.;
  
  // linearization will be done in the flux splitter if needed
  _fluxSplitter->computeFlux(_pertFlux);
  
  if (_hasDiffusiveTerm && _isDiffusionActive) {
    //_nodalExtrapolator->extrapolateInNodes(*_currFace->getNodes());
    _diffusiveFlux->computeFlux(_dFlux);
    _pertFlux -= _dFlux;
  }
  
  // compute the finite difference derivative of the flux
  _numericalJacob->computeDerivative(getJacobianFlux(),_pertFlux,_fluxDiff);
}

//////////////////////////////////////////////////////////////////////////////

} // namespace FiniteVolume

} // namespace Numerics

} // namespace COOLFluiD
