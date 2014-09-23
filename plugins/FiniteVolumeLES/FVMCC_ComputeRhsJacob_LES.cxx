#include "Framework/CFL.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/BlockAccumulator.hh"
#include "Framework/MeshData.hh"

#include "FiniteVolume/FVMCC_BC.hh"

#include "FiniteVolumeLES/FiniteVolumeLES.hh"
#include "FiniteVolumeLES/FVMCC_ComputeRhsJacob_LES.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<FVMCC_ComputeRhsJacob_LES,
                      CellCenterFVMData,
                      FiniteVolumeLESModule>
fvmcc_computeRhsJacob_LES("NumJacobLES");

//////////////////////////////////////////////////////////////////////////////

FVMCC_ComputeRhsJacob_LES::FVMCC_ComputeRhsJacob_LES(const std::string& name) :
  FVMCC_ComputeRhsJacob(name),
  m_lesVarSet()
{
  addConfigOptionsTo(this);
}

//////////////////////////////////////////////////////////////////////////////

FVMCC_ComputeRhsJacob_LES::~FVMCC_ComputeRhsJacob_LES()
{
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacob_LES::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacob_LES::setup()
{
  FVMCC_ComputeRhsJacob::setup();

  // dynamic cast the diffusive varset to LESVarSet
  m_lesVarSet = _diffVar.d_castTo< LES::LESVarSet >();
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacob_LES::configure ( Config::ConfigArgs& args )
{
  FVMCC_ComputeRhsJacob::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacob_LES::computeConvDiffFluxes(CFuint iVar, 
						  CFuint iCell)
{  
  // get the volumes
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
  
  // extrapolate (and LIMIT, if the reconstruction is linear or more)
  // the solution in the quadrature points
  _polyRec->extrapolate(_currFace, iVar, iCell);
  
  // compute the physical data for each left and right reconstructed
  // state and in the left and right cell centers
  computeStatesData();
  
  // set the perturbed variable
  getMethodData().setIPerturbVar(iVar);
  
  // linearization will be done in the flux splitter if needed
  _fluxSplitter->computeFlux(_pertFlux);
  
  if (_hasDiffusiveTerm) {
    //_nodalExtrapolator->extrapolateInNodes(*_currFace->getNodes());
	  // set the volume in the LESVarSet
    CFreal volume;
    volume  = volumes[_currFace->getState(0)->getLocalID()];
    volume += volumes[_currFace->getState(1)->getLocalID()];
    volume *= 0.5;
    m_lesVarSet->setVolume(volume);
    
    _diffusiveFlux->computeFlux(_dFlux);
    _pertFlux -= _dFlux;
  }
    
  // compute the finite difference derivative of the flux
  _numericalJacob->computeDerivative(getJacobianFlux(),_pertFlux,_fluxDiff);
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacob_LES::computeBoundaryJacobianTerm()
{
  // get the volumes
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();

  getMethodData().setIsPerturb(true);

  // set the index of the block corresponding to the current
  // state in the jacobian matrix
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  State& currState = *_currFace->getState(0);
  State& ghostState = *_currFace->getState(1);
  cf_assert(ghostState.isGhost());

  if (currState.isParUpdatable()) {
    const bool isAxi = getMethodData().isAxisymmetric();
    _upFactor[LEFT] = (!isAxi) ? getResFactor() : getResFactor()*(_rMid*_invr[0]);
    _upStFactor[LEFT] = (!isAxi) ? -getResFactor() :-getResFactor()*_invr[0];
    
    // copy the original value of the ghost state
    _origState = ghostState;
    
    _bAcc->setRowColIndex(0, currState.getLocalID());
    
    for (CFuint iVar = 0; iVar < nbEqs; ++iVar) {
      CFLogDebugMax( "Perturbing iVar = " << iVar << "\n");

      // perturb the given component of the state vector
      _numericalJacob->perturb(iVar, currState[iVar]);

      // compute the ghost state in the perturbed inner state
      _currBC->setGhostState(_currFace);

      // extrapolate (and LIMIT, if the reconstruction is linear or more)
      // the solution in the quadrature points
      _polyRec->extrapolate(_currFace);

      // compute the physical data for each left and right reconstructed
      // state and in the left and right cell centers
      computeStatesData();

      // set the perturbed variable
      getMethodData().setIPerturbVar(iVar);
      
      _currBC->computeFlux(_pertFlux);

      if (_hasDiffusiveTerm) {
	      //_nodalExtrapolator->extrapolateInNodes(*_currFace->getNodes());
    	  // set the volume in the LESVarSet
        CFreal volume;
        volume  = volumes[_currFace->getState(0)->getLocalID()];
        volume += volumes[_currFace->getState(1)->getLocalID()];
        volume *= 0.5;
        m_lesVarSet->setVolume(volume);
        _diffusiveFlux->computeFlux(_dFlux);
        _pertFlux -= _dFlux;
      }

      // compute the finite difference derivative of the flux
      _numericalJacob->computeDerivative(getJacobianFlux(),_pertFlux,_fluxDiff);

      addJacobTerm(0, iVar, 0, _bAcc.get());

      // restore the unperturbed value
      _numericalJacob->restore(currState[iVar]);

      // restore the original ghost state
      ghostState = _origState;
    }

    // compute analytical jacobian for source term 
     if (computeSourceTermJacob(LEFT,_stAnJacobIDs)) {
       addAnalyticSourceTermJacob(LEFT, _bAcc.get());
     }
    
    // add the values in the jacobian matrix
    _lss->getMatrix()->addValues(*_bAcc);

    // reset to zero the entries in the block accumulator
    _bAcc->reset();
    _sourceJacobOnCell[LEFT] = false;
  }
}

//////////////////////////////////////////////////////////////////////////////

} // namespace FiniteVolume

} // namespace Numerics

} // namespace COOLFluiD
