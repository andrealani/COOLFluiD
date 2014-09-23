#include "Environment/ObjectProvider.hh"

#include "Common/BadValueException.hh"
#include "Framework/GeometricEntity.hh"
#include "Framework/MeshData.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "LinearAdv/AdvectionDiffusionVarSet.hh"
#include "LinearAdv/LinearAdv2DVarSet.hh"

#include "FluctSplit/ScalarDiffusionTerm.hh"
#include "FluctSplit/InwardNormalsData.hh"
#include "FluctSplit/FluctSplitAdvectionDiffusion.hh"
//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::LinearAdv;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<ScalarDiffusionTerm,
           FluctuationSplitData,
           ComputeDiffusiveTerm,
           FluctSplitAdvectionDiffusionModule>
sdDiffusiveTermProvider("ScalarDiffusion");

//////////////////////////////////////////////////////////////////////////////

ScalarDiffusionTerm::ScalarDiffusionTerm(const std::string& name) :
  ComputeDiffusiveTerm(name),
  _diffVar(CFNULL),
  _updateVar(CFNULL),
  _states(),
  _values(),
  _gradients(),
  _avValues(CFNULL),
  _normal()
{
}

//////////////////////////////////////////////////////////////////////////////

ScalarDiffusionTerm::~ScalarDiffusionTerm()
{
  for (CFuint i = 0; i< _values.size(); ++i) {
    deletePtr(_values[i]);
  }

  for (CFuint i = 0; i< _gradients.size(); ++i) {
    deletePtr(_gradients[i]);
  }

  deletePtr(_avValues);
}

//////////////////////////////////////////////////////////////////////////////

void ScalarDiffusionTerm::setDiffusiveVarSet(SafePtr<DiffusiveVarSet> diffVar)
{
  _diffVar = diffVar.d_castTo<Physics::LinearAdv::AdvectionDiffusionVarSet>();
 }

//////////////////////////////////////////////////////////////////////////////

void ScalarDiffusionTerm::setUpdateVarSet(SafePtr<ConvectiveVarSet> updateVar)
{
  _updateVar = updateVar.d_castTo<Physics::LinearAdv::LinearAdv2DVarSet>();
}

//////////////////////////////////////////////////////////////////////////////

void ScalarDiffusionTerm::computeDiffusiveTerm
(GeometricEntity *const geo, vector<RealVector>& result, bool updateCoeffFlag)
{
  DataHandle< InwardNormalsData*> normals = socket_normals.getDataHandle();
  DataHandle< CFreal> updateCoeff = socket_updateCoeff.getDataHandle();

  vector<State*> *const cellStates = geo->getStates();
  const CFuint nbCellStates = cellStates->size();

  // store the pointers to state in another array (of RealVector*)
  for (CFuint i = 0; i < nbCellStates; ++i) {
    _states[i] = (*cellStates)[i];
  }

  // compute vars that will be used to compute the gradients
  _diffVar->setGradientVars(_states, _values, geo->nbNodes());

  // compute the radius (axysimmetric computations)
  CFreal radius = 0.0;
  if (getMethodData().isAxisymmetric()) {
    for (CFuint i = 0; i < nbCellStates; ++i) {
      const Node& node = *geo->getNode(i);
      radius += node[YY];
    }
    radius /= nbCellStates;
  }

  const CFreal cellVolume = geo->computeVolume();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFreal dimCoeff = 1./dim;
  const CFreal coeffGrad = dimCoeff/cellVolume;

  // gradient and average values computation
  const CFuint cellID = geo->getID();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

  for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
    RealVector& grad = *_gradients[iEq];
    const RealVector& values = *_values[iEq];

    grad = 0.0;

    for (CFuint is = 0; is < nbCellStates; ++is) {
      for (CFuint iDim = 0; iDim < dim; ++iDim) {
  grad[iDim] += values[is]*normals[cellID]->getNodalNormComp(is,iDim);
      }
    }
    grad *= coeffGrad;

 }


  ADTerm& model = _diffVar->getModel();
  const CFreal ovDimCoeff2 = (1./dimCoeff*dimCoeff);

  // set the diffusive term
  for (CFuint i = 0; i < nbCellStates; ++i) {
    // this is not the unit normal !!
    for (CFuint iDim = 0; iDim < dim; ++iDim) {
      _normal[iDim] = normals[cellID]->getNodalNormComp(i,iDim);
    }
    result[i] = (-dimCoeff)*_diffVar->getFlux
      (*_avValues, _gradients, _normal, radius);

    if (updateCoeffFlag) {
      const CFreal faceArea = normals[cellID]->getAreaNode(i);
      const CFreal mu = (model.getPhysicalData())[ADTerm::NU];

     const CFreal diffCoeff = mu;
      updateCoeff[geo->getState(i)->getLocalID()] +=
  diffCoeff*faceArea*faceArea/(cellVolume*ovDimCoeff2);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ScalarDiffusionTerm::setup()
{
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbNodesInControlVolume =
    MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell();
  _states.resize(nbNodesInControlVolume);
  _values.resize(nbEqs);
CF_DEBUG_POINT;
  for (CFuint i = 0; i< nbEqs; ++i) {
    _values[i] = new RealVector(nbNodesInControlVolume);
  }
CF_DEBUG_POINT;
  _gradients.resize(nbEqs);
  for (CFuint i = 0; i< nbEqs; ++i) {
    _gradients[i] = new RealVector(PhysicalModelStack::getActive()->getDim());
  }

  _normal.resize(PhysicalModelStack::getActive()->getDim());
CF_DEBUG_POINT;
  _avValues = new State();
  normaljState.resize(PhysicalModelStack::getActive()->getDim());
CF_DEBUG_POINT;
}

//////////////////////////////////////////////////////////////////////////////

void ScalarDiffusionTerm::configure ( Config::ConfigArgs& args )
{
  ComputeDiffusiveTerm::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void ScalarDiffusionTerm::computePicardDiffJacob(GeometricEntity *const geo, std::vector<RealMatrix*>& jacob){
  // carefull with the signs !!!
 DataHandle< InwardNormalsData*> normals = socket_normals.getDataHandle();
 vector<State*> *const cellStates = geo->getStates();
 const CFuint nbCellStates = cellStates->size();
 const CFuint cellID = geo->getID();
 const CFuint dim = PhysicalModelStack::getActive()->getDim();
 CFreal nxi;
 CFreal nyi;
  const CFreal cellVolume = geo->computeVolume();
ADTerm& model = _diffVar->getModel();
const CFreal mu = (model.getPhysicalData())[ADTerm::NU];

 for (CFuint iState = 0; iState < nbCellStates; ++iState) {
  nxi = normals[cellID]->getNodalNormComp(iState,0);
    nyi = normals[cellID]->getNodalNormComp(iState,1);

   const CFuint nStart = iState*nbCellStates;
    for (CFuint jState = 0; jState < nbCellStates; ++jState) {
      RealMatrix *const block = jacob[nStart + jState];
      normaljState[XX] =  normals[cellID]->getNodalNormComp(jState,0);
      normaljState[YY] =  normals[cellID]->getNodalNormComp(jState,1);
      (*block) = mu*(nxi*normaljState[XX]+nyi*normaljState[YY])/(4.0*cellVolume);
    }
  }
}
//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
