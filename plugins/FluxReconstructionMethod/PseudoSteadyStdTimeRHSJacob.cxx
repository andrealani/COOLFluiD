#include "Framework/CFL.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/LSSIdxMapping.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/MeshData.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/BlockAccumulator.hh"

#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"
#include "FluxReconstructionMethod/PseudoSteadyStdTimeRHSJacob.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< PseudoSteadyStdTimeRHSJacob,
                       FluxReconstructionSolverData,
                       FluxReconstructionModule>
PseudoSteadyRHSJacob("PseudoSteadyTimeRHS");

//////////////////////////////////////////////////////////////////////////////

PseudoSteadyStdTimeRHSJacob::PseudoSteadyStdTimeRHSJacob(const std::string& name) :
  FluxReconstructionSolverCom(name),
  socket_volumes("volumes",false),
  socket_rhs("rhs"),
  socket_updateCoeff("updateCoeff"),
  socket_pastStates("pastStates"),
  m_cellBuilder(CFNULL),
  m_lss(CFNULL),
  m_jacobMatrix(CFNULL),
  m_nbrEqs(),
  m_cellStates(CFNULL),
  m_solPntsLocalCoords(CFNULL),
  m_diagValues(),
  m_isUnsteady(),
  m_tempState(),
  m_updateToSolutionVecTrans(CFNULL),
  m_acc(CFNULL),
  m_numericalJacob(CFNULL),
  m_fluxDiff(),
  m_nbrSolPnts()
{
//  m_useGlobalDT = false;
//  setParameter("useGlobalDT",&m_useGlobalDT);
//
//  m_useAnalyticalMatrix = false;
//  setParameter("useAnalyticalMatrix",&m_useAnalyticalMatrix);
}

//////////////////////////////////////////////////////////////////////////////

PseudoSteadyStdTimeRHSJacob::~PseudoSteadyStdTimeRHSJacob()
{
}

//////////////////////////////////////////////////////////////////////////////

void PseudoSteadyStdTimeRHSJacob::defineConfigOptions(Config::OptionList& options)
{
//  options.addConfigOption< bool >("useGlobalDT", "Flag telling if to use global DT.");
//  options.addConfigOption< bool >("useAnalyticalMatrix", "Flag telling if to use analytical matrix.");  
}

//////////////////////////////////////////////////////////////////////////////

void PseudoSteadyStdTimeRHSJacob::setup()
{
  // setup parent class
  FluxReconstructionSolverCom::setup();

  // get number of equations in the physical model
  m_nbrEqs = PhysicalModelStack::getActive()->getNbEq();

  // get cell builder
  m_cellBuilder = getMethodData().getStdTrsGeoBuilder();

  // linear system solver
  m_lss = getMethodData().getLinearSystemSolver()[0];

  // get maximum number of solution points in a cell
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
  const CFuint nbrElemTypes = frLocalData.size();
  CFuint maxNbrSolPnts = 0;
  for (CFuint iElemType = 0; iElemType < nbrElemTypes; ++iElemType)
  {
    const CFuint nbrSolPnts = frLocalData[iElemType]->getNbrOfSolPnts();
    maxNbrSolPnts = maxNbrSolPnts > nbrSolPnts ? maxNbrSolPnts : nbrSolPnts;
  }
  m_nbrSolPnts = maxNbrSolPnts;

  // resize m_diagValues
  m_diagValues.resize(maxNbrSolPnts);
  
  m_tempState.resize(m_nbrEqs);
  
  m_updateToSolutionVecTrans = getMethodData().getUpdateToSolutionVecTrans();
  
  m_acc.reset(m_lss->createBlockAccumulator(maxNbrSolPnts,maxNbrSolPnts,m_nbrEqs));
  
  // numerical jacobian
  m_numericalJacob = getMethodData().getNumericalJacobian();
  
  m_fluxDiff.resize(m_nbrEqs);
}

//////////////////////////////////////////////////////////////////////////////

void PseudoSteadyStdTimeRHSJacob::unsetup()
{
  // unsetup parent class
  FluxReconstructionSolverCom::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void PseudoSteadyStdTimeRHSJacob::execute()
{
  CFAUTOTRACE;

  // get jacobian matrix
  m_jacobMatrix = m_lss->getMatrix();

  // get cfl number and time step
  const CFreal cfl = getMethodData().getCFL()->getCFLValue();
  const CFreal dt = SubSystemStatusStack::getActive()->getDT();

  // check if computation is unsteady (time accurate)
  m_isUnsteady = dt > 0.;
  
  m_acc.reset(m_lss->createBlockAccumulator(m_nbrSolPnts,m_nbrSolPnts,m_nbrEqs));

  // get datahandle of volumes if necessary
  DataHandle<CFreal> volumes(CFNULL);
  if(m_isUnsteady)
  {
    cf_assert(socket_volumes.isConnected());
    volumes = socket_volumes.getDataHandle();
  }

  // get update coefficients
  DataHandle< CFreal> updateCoeff = socket_updateCoeff.getDataHandle();

  // get the elementTypeData
  SafePtr< vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();

  // get InnerCells TopologicalRegionSet
  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");

  // get the geodata of the geometric entity builder and set the TRS
  StdTrsGeoBuilder::GeoData& geoData = m_cellBuilder->getDataGE();
  geoData.trs = cells;

  // loop over element types
  const CFuint nbrElemTypes = elemType->size();
  for (CFuint iElemType = 0; iElemType < nbrElemTypes; ++iElemType)
  {
    // get start and end indexes for this type of element
    const CFuint startIdx = (*elemType)[iElemType].getStartIdx();
    const CFuint endIdx   = (*elemType)[iElemType].getEndIdx();

    // get the local spectral FD data
    vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();

    // get solution point local coordinates
    m_solPntsLocalCoords = frLocalData[iElemType]->getSolPntsLocalCoords();

    // add the diagonal entries in the jacobian (updateCoeff/CFL)
    for (CFuint elemIdx = startIdx; elemIdx < endIdx; ++elemIdx)
    {
      // build the GeometricEntity
      geoData.idx = elemIdx;
      GeometricEntity* cell = m_cellBuilder->buildGE();

      // get the states in this cell
      m_cellStates = cell->getStates();

      // compute diagonal values
      if ((*m_cellStates)[0]->isParUpdatable())
      {
        const CFuint nbrSolPnts = m_cellStates->size();
        cf_assert(m_diagValues.size() >= nbrSolPnts);
        if (m_isUnsteady)
        {
          for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
          {
            m_diagValues[iSol] = volumes[(*m_cellStates)[iSol]->getLocalID()]/dt;
          }
        }
        else
        {
          // get the cell volume
          const CFreal invCellVolume = 1.0/cell->computeVolume();

          // get jacobian determinants at solution points
          const std::valarray<CFreal> jacobDet =
              cell->computeGeometricShapeFunctionJacobianDeterminant(*m_solPntsLocalCoords);

          const CFuint firstStateID = (*m_cellStates)[0]->getLocalID();
          const CFreal updateCoeffDivCFL = invCellVolume*updateCoeff[firstStateID]/cfl;// should be the same for all states
          for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
          {
            m_diagValues[iSol] = updateCoeffDivCFL*jacobDet[iSol];
          }
        }

        // add the contribution of the time residual to the rhs and the jacobian
        addTimeResidual();
    
      }

      //release the GeometricEntity
      m_cellBuilder->releaseGE();
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void PseudoSteadyStdTimeRHSJacob::addTimeResidual()
{
  // get factor for the residual
  const CFreal resFactor = getMethodData().getResFactor();

  // multiply diagonal values with the residual factor
  m_diagValues *= resFactor;

  // get rhs and past states
  DataHandle< CFreal> rhs = socket_rhs.getDataHandle();
  DataHandle<State*> pastStates = socket_pastStates.getDataHandle();

  // local (current processor) to global mapping
  const LSSIdxMapping& idxMapping =
      getMethodData().getLinearSystemSolver()[0]->getLocalToGlobalMapping();

  // add time residual contribution (Backward Euler)
  const CFuint nbrSolPnts = m_cellStates->size();
  
  for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
  {
    // now the contribution to the jacobian matrix is calculated
    m_acc->setRowColIndex(iSol, (*m_cellStates)[iSol]->getLocalID());
  }
    
  for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
  {
    // get state
    State *const currState = (*m_cellStates)[iSol];

    // get state ID
    const CFuint stateID = currState->getLocalID();

    // get past state
    const State& pastState = *pastStates[stateID];
    
    // this first transformed state HAS TO BE stored,
    // since the returned pointer will change pointee after
    // the second call to transform()
    m_tempState = static_cast<RealVector&>(*m_updateToSolutionVecTrans->transform(currState));

    const State *const tPastState = m_updateToSolutionVecTrans->transform(pastStates[stateID]);

    // add contribution to rhs and jacobian
    CFuint globalID = idxMapping.getColID(stateID)*m_nbrEqs;
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq, ++globalID)
    {
      rhs(stateID, iEq, m_nbrEqs) -= (m_tempState[iEq] - (*tPastState)[iEq])*m_diagValues[iSol];
//       CF_DEBUG_OBJ(rhs(stateID, iEq, m_nbrEqs));

      if((!getMethodData().isSysMatrixFrozen()) && getMethodData().doComputeJacobian())
      {
        // perturb the given component of the state vector
        m_numericalJacob->perturb(iEq, (*currState)[iEq]);

        const RealVector& tempPertState = static_cast<RealVector&>(*m_updateToSolutionVecTrans->transform(currState));

        // compute the finite difference derivative of the flux
        m_numericalJacob->computeDerivative(m_tempState,tempPertState,m_fluxDiff);

        // _fluxDiff corresponds to a column vector of the dU/dP matrix
        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq) 
        {
	  m_fluxDiff[iEq] *= m_diagValues[iSol];
        }

        m_acc->addValues(iSol, iSol, iEq, &m_fluxDiff[0]);

        // restore the unperturbed value
        m_numericalJacob->restore((*currState)[iEq]);

      }
    }

    if((!getMethodData().isSysMatrixFrozen()) && getMethodData().doComputeJacobian())
    {
      // add the values in the jacobian matrix
      //getMethodData().getLSSMatrix(0)->addValues(*m_acc);
      m_lss->getMatrix()->addValues(*m_acc);

      // reset to zero the entries in the block accumulator
      m_acc->reset();
    }

    //m_jacobMatrix->addValue(globalID, globalID, m_diagValues[iSol]);
  }
//   CF_DEBUG_POINT;
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > PseudoSteadyStdTimeRHSJacob::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_rhs);
  result.push_back(&socket_volumes);
  result.push_back(&socket_updateCoeff);
  result.push_back(&socket_pastStates);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD
