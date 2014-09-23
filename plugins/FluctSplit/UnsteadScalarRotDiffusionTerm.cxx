#include "Environment/ObjectProvider.hh"

#include "Common/BadValueException.hh"
#include "Framework/GeometricEntity.hh"
#include "Framework/MeshData.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "RotationAdv/RotationDiffusionVarSet.hh"
#include "RotationAdv/RotationAdv2DVarSet.hh"
#include "Framework/SubSystemStatus.hh"

#include "FluctSplit/UnsteadScalarRotDiffusionTerm.hh"
#include "FluctSplit/InwardNormalsData.hh"
#include "FluctSplit/FluctSplitRotationDiffusion.hh"
//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::RotationAdv;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<UnsteadScalarRotDiffusionTerm,
		       FluctuationSplitData,
		       ComputeDiffusiveTerm,
		       FluctSplitRotationDiffusionModule>
unsrdDiffusiveTermProvider("UnsteadScalarRotDiffusion");

//////////////////////////////////////////////////////////////////////////////

UnsteadScalarRotDiffusionTerm::UnsteadScalarRotDiffusionTerm(const std::string& name) :
  ComputeDiffusiveTerm(name),
  _diffVar(CFNULL),
  _updateVar(CFNULL),
  _states(),
  _values(),
  _gradients(),
  _avValues(CFNULL),
  _normal(),
  socket_pastStates("pastStates"),
  _pastStates(0)
{
}

//////////////////////////////////////////////////////////////////////////////

UnsteadScalarRotDiffusionTerm::~UnsteadScalarRotDiffusionTerm()
{
  for (CFuint i = 0; i< _values.size(); ++i) {
    deletePtr(_values[i]);
  }

  for (CFuint i = 0; i< _gradients.size(); ++i) {
    deletePtr(_gradients[i]);
  }
  
  deletePtr(_avValues); 

// deallocating data for the temporary local residual
  for (CFuint i = 0; i < _pastStates.size(); ++i) {
    deletePtr(_pastStates[i]);
  }
}

//////////////////////////////////////////////////////////////////////////////

void UnsteadScalarRotDiffusionTerm::setDiffusiveVarSet(SafePtr<DiffusiveVarSet> diffVar)
{
  _diffVar = diffVar.d_castTo<RotationDiffusionVarSet>();
 }

//////////////////////////////////////////////////////////////////////////////

void UnsteadScalarRotDiffusionTerm::setUpdateVarSet(SafePtr<ConvectiveVarSet> updateVar)
{
  _updateVar = updateVar.d_castTo<RotationAdv2DVarSet>();
}

//////////////////////////////////////////////////////////////////////////////

void UnsteadScalarRotDiffusionTerm::computeDiffusiveTerm
(GeometricEntity *const geo, vector<RealVector>& result, bool updateCoeffFlag)
{
  DistributionData& ddata = getMethodData().getDistributionData();
  DataHandle< InwardNormalsData*> normals = socket_normals.getDataHandle();
  DataHandle< CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
  
  // Start with the part from the past 
  DataHandle<State*> pastStatesStorage = socket_pastStates.getDataHandle();
  const CFuint nbStatesInCell = ddata.states->size();
  const CFreal cellVolume = geo->computeVolume();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFreal dimCoeff = 1./(dim);
  const CFreal coeffGrad = 1.0/(cellVolume*(dim));
  const CFreal dt = SubSystemStatusStack::getActive()->getDT()/dim;
  const CFuint cellID = geo->getID();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
 RDTerm& model = _diffVar->getModel();
  const CFreal ovDimCoeff2 = 1.0/(dimCoeff*dimCoeff);


  // get the paststates in this cell
  for (CFuint i = 0; i < nbStatesInCell; ++i) {

    const CFuint stateID = (*ddata.states)[i]->getLocalID();

    *_pastStates[i] = *pastStatesStorage[stateID];
  }

  // store the pointers to state in another array (of RealVector*)
  for (CFuint i = 0; i < nbStatesInCell; ++i) {
    _states[i] = (_pastStates)[i];
  }

  // compute vars that will be used to compute the gradients
  _diffVar->setGradientVars(_states, _values, geo->nbNodes());

 
  // compute the radius (axysimmetric computations)
  CFreal radius = 0.0;
  if (getMethodData().isAxisymmetric()) {
    for (CFuint i = 0; i < nbStatesInCell; ++i) {
      const Node& node = *geo->getNode(i);
      radius += node[YY];
    }
    radius /= nbStatesInCell;
  }

  // gradient and average values computation


  for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
    RealVector& grad = *_gradients[iEq];
    const RealVector& values = *_values[iEq];

    grad = 0.0;

    for (CFuint is = 0; is < nbStatesInCell; ++is) {
      for (CFuint iDim = 0; iDim < dim; ++iDim) {
	grad[iDim] += values[is]*normals[cellID]->getNodalNormComp(is,iDim);
      }
    }
    grad *= coeffGrad;

 }

  // set the diffusive term
  for (CFuint i = 0; i < nbStatesInCell; ++i) {
    // this is not the unit normal !!
    for (CFuint iDim = 0; iDim < dim; ++iDim) {
      _normal[iDim] = normals[cellID]->getNodalNormComp(i,iDim);
    }
    result[i] = dt*(-dimCoeff)*_diffVar->getFlux
      (*_avValues, _gradients, _normal, radius);
 
  }

  // Now contribution from the present
  // store the pointers to state in another array (of RealVector*)
  vector<State*> *const cellStates = geo->getStates();

  for (CFuint i = 0; i < nbStatesInCell; ++i) {
    _states[i] = (*cellStates)[i];
  }

  // compute vars that will be used to compute the gradients
  _diffVar->setGradientVars(_states, _values, geo->nbNodes());

  
  for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
    RealVector& grad = *_gradients[iEq];
    const RealVector& values = *_values[iEq];

    grad = 0.0;

    for (CFuint is = 0; is < nbStatesInCell; ++is) {
      for (CFuint iDim = 0; iDim < dim; ++iDim) {
	grad[iDim] += values[is]*normals[cellID]->getNodalNormComp(is,iDim);
      }
    }
    grad *= coeffGrad;

 }
 
  // set the diffusive term
  for (CFuint i = 0; i < nbStatesInCell; ++i) {
    // this is not the unit normal !!
    for (CFuint iDim = 0; iDim < dim; ++iDim) {
      _normal[iDim] = normals[cellID]->getNodalNormComp(i,iDim);
    }
    result[i] = dt*(-dimCoeff)*_diffVar->getFlux
      (*_avValues, _gradients, _normal, radius);
 
    if (updateCoeffFlag) {
      const CFreal faceArea = normals[cellID]->getAreaNode(i);
       const CFreal mu = (model.getPhysicalData())[RDTerm::NU];

     const CFreal diffCoeff = mu;
      updateCoeff[geo->getState(i)->getLocalID()] +=
	diffCoeff*faceArea*faceArea/(cellVolume*ovDimCoeff2);
    }
  }


}

//////////////////////////////////////////////////////////////////////////////

void UnsteadScalarRotDiffusionTerm::setup()
{  
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbNodesInControlVolume =
    MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell();
  _states.resize(nbNodesInControlVolume);
  _values.resize(nbEqs);
  
  for (CFuint i = 0; i< nbEqs; ++i) {
    _values[i] = new RealVector(nbNodesInControlVolume);
  }
  
  _gradients.resize(nbEqs);
  for (CFuint i = 0; i< nbEqs; ++i) {
    _gradients[i] = new RealVector(PhysicalModelStack::getActive()->getDim());
  }
  
  _normal.resize(PhysicalModelStack::getActive()->getDim());
  
  _avValues = new State();

const CFuint maxNbStatesInCell = MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell();

  _pastStates.resize(maxNbStatesInCell);
  // Resizing pastStates
  for (CFuint i = 0; i < maxNbStatesInCell; ++i) {
    _pastStates[i] = new State();
   }
}


//////////////////////////////////////////////////////////////////////////////

void UnsteadScalarRotDiffusionTerm::configure ( Config::ConfigArgs& args )
{
  ComputeDiffusiveTerm::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit


} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
