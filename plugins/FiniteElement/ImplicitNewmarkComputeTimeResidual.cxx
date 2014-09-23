#include "Framework/MeshData.hh"
#include "Framework/VolumeIntegrator.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/BlockAccumulator.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MethodCommandProvider.hh"

#include "FiniteElement/FiniteElement.hh"
#include "FiniteElement/ImplicitNewmarkComputeTimeResidual.hh"
#include "FiniteElement/ComputeInertiaTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<ImplicitNewmarkComputeTimeResidual, FiniteElementMethodData, FiniteElementModule> ImplicitNewmarkComputeTimeResidualProvider("ImplicitNewmarkComputeTimeResCom");

//////////////////////////////////////////////////////////////////////////////

void ImplicitNewmarkComputeTimeResidual::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("Alpha","Alpha coef");
   options.addConfigOption< CFreal >("Gamma","Gamma coef");
   options.addConfigOption< bool >("LumpedMassMatrix","Use or not a lumped Lumped Mass Matrix using sum row technique");
}

//////////////////////////////////////////////////////////////////////////////

ImplicitNewmarkComputeTimeResidual::ImplicitNewmarkComputeTimeResidual(const std::string& name) :
FiniteElementMethodCom(name),
  socket_rhs("rhs"),
  socket_pastStates("pastStates"),
  socket_pastStatesD("pastStatesD"),
  socket_pastStatesD2("pastStatesD2"),
  _pastStates(CFNULL),
  _pastStatesD(CFNULL),
  _pastStatesD2(CFNULL),
  _inertiaTerm(CFNULL)
{
   addConfigOptionsTo(this);
  _alpha = 0.5;
   setParameter("Alpha",&_alpha);

  _gamma = 0.5;
   setParameter("Gamma",&_gamma);

  _lumpMassMatrix = false;
   setParameter("LumpedMassMatrix",&_lumpMassMatrix);

}

//////////////////////////////////////////////////////////////////////////////

ImplicitNewmarkComputeTimeResidual::~ImplicitNewmarkComputeTimeResidual()
{
}

//////////////////////////////////////////////////////////////////////////////

void ImplicitNewmarkComputeTimeResidual::unsetup()
{
  for(CFuint i = 0; i < m_map_femdata.size(); ++i)
  {
    deletePtr(m_map_femdata[i].first);
    deletePtr(m_map_femdata[i].second);
    deletePtr(m_map_femdata[i].third);
    deletePtr(m_map_femdata[i].fourth);
  }

  // then call parent method
  FiniteElementMethodCom::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void ImplicitNewmarkComputeTimeResidual::setup()
{
  CFAUTOTRACE;

  // first call parent method
  FiniteElementMethodCom::setup();

  // rhs storage
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();

  SafePtr<vector<ElementTypeData> > elementType =
    MeshDataStack::getActive()->getElementTypeData();

  const CFuint nbElemTypes = elementType->size();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

  // settin correct size for integration result
  _integResultMat.resize(nbEqs,nbEqs);

  // loop over types since it can happen to deal with an hybrid mesh
  for (CFuint iType = 0; iType < nbElemTypes; ++iType) {
     // number of states in a cell of this type
     const CFuint nbStatesInType = (*elementType)[iType].getNbStates();
     // create a block matrix accumulator for this type of cell
     BlockAccumulator* ptr = getMethodData().getLinearSystemSolver()[0]->
       createBlockAccumulator(nbStatesInType,nbStatesInType,nbEqs);
     RealVector* vec = new RealVector(nbStatesInType*nbEqs);
     RealVector* vec2 = new RealVector(nbEqs);
     RealMatrix* mat = new RealMatrix(nbStatesInType*nbEqs,nbStatesInType*nbEqs);
     vector<RealVector>* residual = new vector<RealVector>(nbStatesInType,*vec2);

     FElemTypeData elemTypeData(ptr,mat,vec,residual);

    m_map_femdata.insert(nbStatesInType, elemTypeData);
  }

  m_map_femdata.sortKeys();


  // storage of the past States
  _pastStates  = socket_pastStates.getDataHandle();

  // storage of the past States derivatives
  _pastStatesD  = socket_pastStatesD.getDataHandle();

  // storage of the past States second derivatives
  _pastStatesD2  = socket_pastStatesD2.getDataHandle();

  _inertiaTerm = getMethodData().getInertiaTermComputer();


  _tempRes.resize(PhysicalModelStack::getActive()->getNbEq());
  _otherResidual.resize(MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell());
  for (CFuint i = 0; i < _otherResidual.size(); ++i) {
    _otherResidual[i].resize(PhysicalModelStack::getActive()->getNbEq());
    _otherResidual[i] = 0.;
  }


}

//////////////////////////////////////////////////////////////////////////////

void ImplicitNewmarkComputeTimeResidual::executeOnTrs()
{
  CFAUTOTRACE;

  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();

  SafePtr<LSSMatrix> jacobMatrix =
    getMethodData().getLinearSystemSolver()[0]->getMatrix();

  SafePtr<TopologicalRegionSet> trs = getCurrentTRS();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbGeos = trs->getLocalNbGeoEnts();
  FiniteElementMethodData& femdata  = getMethodData();
  LocalElementData& local_elem_data = femdata.getLocalElementData();
  local_elem_data.nbEqs = nbEqs;
  local_elem_data.trs   = getCurrentTRS();
  PhysicalModelStack::getActive()->setCurrentZone(getCurrentTRS()->getName());

  Common::SafePtr<GeometricEntityPool<StdTrsGeoBuilder> >
    geoBuilder = getMethodData().getStdTrsGeoBuilder();

  StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.trs = trs;

  const CFreal dt = SubSystemStatusStack::getActive()->getDT();
  _a3 = 2./(_gamma * dt * dt);
  _a4 = dt * _a3;
  _a5 = (1./_gamma) - 1.;

  for(CFuint iGeoEnt = 0; iGeoEnt < nbGeos; ++iGeoEnt) {

    CFLogDebugMax("Cell " << iGeoEnt << "\n");

    // build the GeometricEntity
    geoData.idx = iGeoEnt;
    GeometricEntity *const cell = geoBuilder->buildGE();

    vector<State*>& states = *(cell->getStates());
    const CFuint nbStatesInCell = states.size();

    local_elem_data.nbStates = nbStatesInCell;
    local_elem_data.cell = cell;

    getMethodData().getFEMVolumeIntegrator()->precomputeCellData(cell);

    /// @todo this must be improved. Finding in a map for every cell is
    /// efficiently speaking not acceptable.
    FElemTypeData elemData = m_map_femdata.find(nbStatesInCell);

    BlockAccumulator& acc = *elemData.first;
    acc.reset();

    RealMatrix& elemMat = *elemData.second;
    RealVector& elemVec = *elemData.third;
    vector<RealVector>& residual = *elemData.fourth;

  local_elem_data.blockacc = &acc;
  local_elem_data.stiff_mat = &elemMat;
  local_elem_data.load_vec = &elemVec;
  local_elem_data.residual = &residual;



    computeElementTimeResidual(*cell,elemMat,elemVec,residual);

    // write to the right hand side
    for (CFuint iState = 0; iState < nbStatesInCell; ++iState) {
      const CFuint stateID = states[iState]->getLocalID();
      for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
        rhs(stateID, iEq, nbEqs) -= residual[iState][iEq];
      }
    }

    // compute the jacobian matrix if it is not frozen
    if(!getMethodData().isSysMatrixFrozen())
    {
      computeJacobianTerm();

      // add the values in the jacobian matrix
      jacobMatrix->addValues(acc);
    }

    //release the GeometricEntity
    geoBuilder->releaseGE();
  }
}

//////////////////////////////////////////////////////////////////////////////

void ImplicitNewmarkComputeTimeResidual::computeElementTimeResidual(
                                      Framework::GeometricEntity& cell,
                                      RealMatrix& elemMat,
                                      RealVector& elemVec,
                                      vector<RealVector>& residual)
{
  CFAUTOTRACE;

  vector<State*>& states = *cell.getStates();
  const CFuint nbStatesInCell = states.size();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  FiniteElementMethodData& femdata  = getMethodData();
  LocalElementData& local_elem_data = femdata.getLocalElementData();

  cf_assert(elemMat.nbRows() == nbStatesInCell * nbEqs);
  cf_assert(elemMat.nbCols() == nbStatesInCell * nbEqs);
  cf_assert(elemVec.size() == nbStatesInCell * nbEqs);
  cf_assert(residual.size() <= nbStatesInCell);

      // Compute and  Insert result into the Block Accumulator
      std::vector<RealVector> inertiaRhs(nbStatesInCell);
      for (CFuint iState = 0; iState < nbStatesInCell; ++iState) {
        inertiaRhs[iState].resize(nbEqs);
      }

  fill(residual.begin(),residual.end(),0.0);
  elemMat = 0.0;
  // compute element matrix and vector
  for (CFuint iState = 0; iState < nbStatesInCell; ++iState) {
    const State& currState = *states[iState];
    const CFuint stateID = currState.getLocalID();
    local_elem_data.iState = iState;
    // Computing the part of the Inertia term going to the RHS (in the explicit case)
    // This term is: M *(a3*u_n + a4*up_n + a5*upp_n)
    // so we first compute: (a3*u_n + a4*up_n + a5*upp_n)
    inertiaRhs[iState] =  _a3 * (*_pastStates[stateID]);
    inertiaRhs[iState] += _a4 * (*_pastStatesD[stateID]);
    inertiaRhs[iState] += _a5 * (*_pastStatesD2[stateID]);

// CFout << "InertiaRhs["<<iState<<"]: "<<inertiaRhs[iState]<<"\n";
    if (currState.isParUpdatable()) {
      for (CFuint jState = 0; jState < nbStatesInCell; ++jState) {
        local_elem_data.jState = jState;
        _inertiaTerm->computeTerm(&cell, _integResultMat);

        for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
          for (CFuint jEq = 0; jEq < nbEqs; ++jEq) {
            if(!_lumpMassMatrix){
              ///This is for consistent mass matrix
              elemMat(iState*nbEqs + iEq, jState*nbEqs + jEq) += _integResultMat(iEq,jEq);
            }
            else{
              ///This is for lumped mass matrix
              elemMat(iState*nbEqs + iEq, iState*nbEqs + jEq) += _integResultMat(iEq,jEq);
            }
          }
        }
      } // for jState
    } // if isParUpdatable
  } // for iState


  // compute the first part of residual: (a3*M*u)
  // compute the second part of residual: -M*(a3*un + a4*upn + a5*uppn)
  /// @todo try to make this product better
  for (CFuint iState = 0; iState < nbStatesInCell; ++iState) {
    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      for (CFuint jState = 0; jState < nbStatesInCell; ++jState) {
        State& jStateRef = *(states[jState]);
        for (CFuint jEq = 0; jEq < nbEqs; ++jEq) {
          residual[iState][iEq] += _a3 * elemMat(iState*nbEqs + iEq, jState*nbEqs + jEq) * jStateRef[jEq];
          residual[iState][iEq] -= elemMat(iState*nbEqs + iEq, jState*nbEqs + jEq) * inertiaRhs[jState][jEq];
        }
      }
    }
  }

}


//////////////////////////////////////////////////////////////////////////////

void ImplicitNewmarkComputeTimeResidual::computeJacobianTerm()
{
  FiniteElementMethodData& femdata  = getMethodData();
  LocalElementData& local_elem_data = femdata.getLocalElementData();

  BlockAccumulator& acc = *local_elem_data.blockacc;
  RealMatrix& elemMat   = *local_elem_data.stiff_mat;
  RealVector& elemVec   = *local_elem_data.load_vec;
  const CFuint nbEqs    =  local_elem_data.nbEqs;
  const CFuint nbStates =  local_elem_data.nbStates;
  GeometricEntity& cell = *local_elem_data.cell;
  std::vector<RealVector>& residual = *local_elem_data.residual;

  NumericalJacobian& numericalJacob = femdata.getNumericalJacobian();

  vector<State*>& states = *cell.getStates();

  // compute the perturbed states for the evaluation of the
  // jacobian matrix looping over the state vectors in this cell
  for (CFuint iState = 0; iState < nbStates; ++iState) {

    State& currState = *states[iState];

    CFLogDebugMax( "Perturbing iState = " << iState << " with stateID = " << currState.getLocalID() << "\n");

    // set the index of the block corresponding to the current
    // state in the jacobian matrix
    acc.setRowColIndex(iState, currState.getLocalID());
    // loop over the variables in the state vector to perturb one
    // component at a time
    for (CFuint iVar = 0; iVar < nbEqs; ++iVar) {
      CFLogDebugMax( "Perturbing iVar = " << iVar << "\n");

      // reset the residual to 0.
      cleanOtherResidual();

      // perturb the given component of the state vector
      numericalJacob.perturb(iVar, currState[iVar]);

      // compute element matrix and vector and put in _otherResidual the perturbed residual
      computeElementTimeResidual(cell,elemMat,elemVec,_otherResidual);

      // compute and distribute the jacobian contributions
      for (CFuint jState = 0; jState < nbStates; ++jState) {
        CFLogDebugMax( "Perturbing jState = " << jState << "\n");

        if (states[jState]->isParUpdatable()) {

          // jacobian contribution (dR_jState/dU_iState)_k
          // compute (R_jState[U_iState + dU_k] - R_jState[U_iState])/eps
          numericalJacob.computeDerivative(residual[jState],
                                           _otherResidual[jState],
                                           _tempRes);
          acc.addValues(jState, iState, iVar, &_tempRes[0]);
        }
      }
      // restore the unperturbed value
      numericalJacob.restore(currState[iVar]);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
ImplicitNewmarkComputeTimeResidual::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_rhs);
  result.push_back(&socket_pastStates);
  result.push_back(&socket_pastStatesD);
  result.push_back(&socket_pastStatesD2);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement

  } // namespace Numerics

} // namespace COOLFluiD
