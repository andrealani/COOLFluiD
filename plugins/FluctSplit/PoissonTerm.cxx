#include "Common/BadValueException.hh"
#include "Framework/GeometricEntity.hh"
#include "Framework/MeshData.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Poisson/PoissonConvVarSet.hh"
#include "Poisson/PoissonDiffVarSet.hh"
#include "FluctSplit/PoissonTerm.hh"
#include "FluctSplit/InwardNormalsData.hh"
#include "FluctSplit/FluctSplitPoisson.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::Poisson;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<PoissonTerm,
		       FluctuationSplitData,
		       ComputeDiffusiveTerm,
		       FluctSplitPoissonModule>
poissonDiffusiveTermProvider("Poisson");

//////////////////////////////////////////////////////////////////////////////

PoissonTerm::PoissonTerm(const std::string& name) :
  ComputeDiffusiveTerm(name),
  _upVar(CFNULL),
  _diffVar(CFNULL),
  _states(),
  _values(),
  _avValues(CFNULL),
  _normal()
{
}

//////////////////////////////////////////////////////////////////////////////

PoissonTerm::~PoissonTerm()
{
  deletePtr(_avValues);
}

//////////////////////////////////////////////////////////////////////////////

void PoissonTerm::setDiffusiveVarSet(Common::SafePtr<Framework::DiffusiveVarSet> diffVar)
{
  _diffVar = diffVar.d_castTo<Physics::Poisson::PoissonDiffVarSet>();
}

//////////////////////////////////////////////////////////////////////////////

void PoissonTerm::setUpdateVarSet(Common::SafePtr<Framework::ConvectiveVarSet> updateVar)
{
  if (ComputeDiffusiveTerm::addToDerivedTerm()) {	
    ComputeDiffusiveTerm::setUpdateVarSet(updateVar);
  }
  _upVar = updateVar.d_castTo<Physics::Poisson::PoissonConvVarSet>();
}

//////////////////////////////////////////////////////////////////////////////

void PoissonTerm::computeDiffusiveTerm(Framework::GeometricEntity *const geo, 
				       std::vector<RealVector>& result,
				       bool updateCoeffFlag)
{
  if (ComputeDiffusiveTerm::addToDerivedTerm()) {			
    ComputeDiffusiveTerm::computeDiffusiveTerm(geo, result, updateCoeffFlag);
  }
  
  // gradient and average values computation
  const RealVector& edata = _upVar->getModel()->getPhysicalData();
  computeCellGradientsAndAverageState(geo, edata);
  
  DataHandle< InwardNormalsData*> normals = this->socket_normals.getDataHandle();
  DataHandle< CFreal> updateCoeff = this->socket_updateCoeff.getDataHandle();
  const CFuint nbCellStates = geo->getStates()->size();
  
  DistributionData& dd = this->getMethodData().getDistributionData();
  
  const CFuint cellID = geo->getID();
  cellVolume = this->socket_volumes.getDataHandle()[cellID];
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFreal dimCoeff = 1./dim;
  const CFreal ovDimCoeff2 = 1./(dimCoeff*dimCoeff);
  
  vector<RealVector*>& gradients = dd.gradients;
  RealVector& avState = dd.avState;

  const CFreal radius = 0;

  // set the diffusive term
  for (CFuint i = 0; i < nbCellStates; ++i) {
    // this is not the unit normal !!
    for (CFuint iDim = 0; iDim < dim; ++iDim) {
      _normal[iDim] = normals[cellID]->getNodalNormComp(i,iDim);
    }
    
    const RealVector& flux = _diffVar->getFlux(avState, gradients, _normal, radius);
    if (ComputeDiffusiveTerm::addToDerivedTerm()) {
      result[i] += (-dimCoeff)*flux;
    }
    else { 
      result[i] = (-dimCoeff)*flux;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void PoissonTerm::setup()
{  
  ComputeDiffusiveTerm::setup();	
  
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbNodesInControlVolume =
    MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell();
  
  _states.resize(nbNodesInControlVolume);
  _values.resize(nbEqs, nbNodesInControlVolume);
  _normal.resize(PhysicalModelStack::getActive()->getDim());
  _avValues = new State();
}

//////////////////////////////////////////////////////////////////////////////

void PoissonTerm::configure ( Config::ConfigArgs& args )
{
  ComputeDiffusiveTerm::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void PoissonTerm::computeCellGradientsAndAverageState
(Framework::GeometricEntity *const geo, const RealVector& pdata) 
{
  DataHandle< InwardNormalsData*> normals = this->socket_normals.getDataHandle();
  vector<State*> *const cellStates = geo->getStates();
  const CFuint nbCellStates = cellStates->size();
  
  // store the pointers to state in another array (of RealVector*)
  for (CFuint i = 0; i < nbCellStates; ++i) {
    _states[i] = (*cellStates)[i];
  }
  
  // compute vars that will be used to compute the gradients
  _diffVar->setGradientVars(_states, _values, geo->nbNodes());
  
  const CFuint cellID = geo->getID();
  const CFreal cellVolume = this->socket_volumes.getDataHandle()[cellID];
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFreal dimCoeff = 1./dim;
  const CFreal coeffGrad = dimCoeff/cellVolume;
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  
  vector<RealVector*>& gradients = this->getMethodData().getDistributionData().gradients;
  
  for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
    RealVector& grad = *gradients[iEq];
    grad = 0.0;
    
    for (CFuint is = 0; is < nbCellStates; ++is) {
      for (CFuint iDim = 0; iDim < dim; ++iDim) {
	grad[iDim] += _values(iEq,is)*normals[cellID]->getNodalNormComp(is,iDim);
      }
    }
    grad *= coeffGrad;
  }
  
  // avesrage state is not needed to be filled in 
  this->getMethodData().getDistributionData().avState = *_avValues;
}
      
//////////////////////////////////////////////////////////////////////////////
      
void PoissonTerm::computePicardDiffJacob(Framework::GeometricEntity *const geo, 
					 std::vector<RealMatrix*>& jacob)
{
}

//////////////////////////////////////////////////////////////////////////////


} // namespace FluctSplit

} // namespace COOLFluiD

///////////////////////////////////////////////////////////////////////////////
