#include "Common/CFLog.hh"
#include "Framework/MethodCommandProvider.hh"

#include "Framework/BlockAccumulator.hh"
#include "Framework/LSSMatrix.hh"

#include "Framework/NamespaceSwitcher.hh"

#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"
#include "FluxReconstructionMethod/StdSourceTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

// MethodCommandProvider<StdSourceTerm, FluxReconstructionSolverData, FluxReconstructionModule>
// StdSourceTermProvider("StdSourceTerm");

//////////////////////////////////////////////////////////////////////////////

StdSourceTerm::StdSourceTerm(const std::string& name) :
  FluxReconstructionSolverCom(name),
  socket_rhs("rhs"),
  m_cell(CFNULL),
  m_cellStates(CFNULL),
  m_nbrEqs(),
  m_iElemType(),
  m_solPntsLocalCoords(CFNULL),
  m_solPntJacobDets(),
  m_lss(CFNULL),
  m_numJacob(CFNULL),
  m_acc(CFNULL),
  m_pertResUpdates(),
  m_resUpdates(),
  m_derivResUpdates(),
  m_nbrSolPnts(),
  m_pertSol(),
  m_isPerturbed(),
  m_useAnaJacob(),
  m_stateJacobian()
{
  addConfigOptionsTo(this);
  
  m_addJacob = false;
  setParameter("AddJacob",&m_addJacob);
  
  m_useAnaJacob = false;
  setParameter("AnalyticalJacob",&m_useAnaJacob); 
}

//////////////////////////////////////////////////////////////////////////////

StdSourceTerm::~StdSourceTerm()
{
}

//////////////////////////////////////////////////////////////////////////////

void StdSourceTerm::defineConfigOptions(Config::OptionList& options)
{
    options.addConfigOption< bool >("AddJacob","Flag telling whether to add the jacobian");
    
    options.addConfigOption< bool >("AnalyticalJacob","Flag telling whether to compute the jacobian analytically (not implemented for all ST)");
}

//////////////////////////////////////////////////////////////////////////////

void StdSourceTerm::setup()
{
  CFAUTOTRACE;
  FluxReconstructionSolverCom::setup();

  // get number of physical variables
  m_nbrEqs = PhysicalModelStack::getActive()->getNbEq();

  // get maximum number of solution points in a cell
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
  const CFuint nbrElemTypes = frLocalData.size();
  CFuint maxNbrSolPnts = 0;
  for (CFuint iElemType = 0; iElemType < nbrElemTypes; ++iElemType)
  {
    const CFuint nbrSolPnts = frLocalData[iElemType]->getNbrOfSolPnts();
    maxNbrSolPnts = maxNbrSolPnts > nbrSolPnts ? maxNbrSolPnts : nbrSolPnts;
  }

  // resize m_solPntJacobDets
  m_solPntJacobDets.resize(maxNbrSolPnts);
  
  m_nbrSolPnts = maxNbrSolPnts;
  
  const CFuint iter = SubSystemStatusStack::getActive()->getNbIter();
    
  const CFuint iterFreeze = getMethodData().getFreezeJacobIter();
    
  const CFuint interval = iter - iterFreeze;
          
  if (m_addJacob)
  { 
    // get the linear system solver
    m_lss = getMethodData().getLinearSystemSolver()[0];

    // get the numerical Jacobian computer
    m_numJacob = getMethodData().getNumericalJacobian();
  }
  
  const CFuint resSize = m_nbrSolPnts*m_nbrEqs;
  
  m_pertResUpdates.resize(resSize);
  m_resUpdates.resize(resSize);
  m_derivResUpdates.resize(resSize);
  
  m_stateJacobian.resize(m_nbrEqs);
  
  for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
  {
    m_stateJacobian[iEq].resize(m_nbrEqs);
  }
}

//////////////////////////////////////////////////////////////////////////////

void StdSourceTerm::unsetup()
{
  CFAUTOTRACE;
  FluxReconstructionSolverCom::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void StdSourceTerm::execute()
{
  CFAUTOTRACE;

  // get the ElementTypeData
  SafePtr< vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();
  const CFuint nbrElemTypes = elemType->size();

  // get inner cells TRS
  SafePtr<TopologicalRegionSet> trs = MeshDataStack::getActive()->getTrs("InnerCells");
  CFLogDebugMin("StdSourceTerm::executeOnTrs() called for TRS: " << trs->getName() << "\n");

  // prepares to loop over cells by getting the GeometricEntityPool
  SafePtr< GeometricEntityPool<StdTrsGeoBuilder> > geoBuilder = getMethodData().getStdTrsGeoBuilder();
  StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.trs = trs;
  
  const CFuint iter = SubSystemStatusStack::getActive()->getNbIter();
    
  const CFuint iterFreeze = getMethodData().getFreezeJacobIter();
    
  const CFuint interval = iter - iterFreeze;
  
  if (m_addJacob)
  {
    // create blockaccumulator
    m_acc.reset(m_lss->createBlockAccumulator(m_nbrSolPnts,m_nbrSolPnts,m_nbrEqs));
  }

  // loop over elements
  for (m_iElemType = 0; m_iElemType < nbrElemTypes; ++m_iElemType)
  {
    // get the number of elements
    const CFuint nbrElems = (*elemType)[m_iElemType].getNbElems();

    // get start index of this element type in global element list
    CFuint cellIdx = (*elemType)[m_iElemType].getStartIdx();

    // get data needed for source term computation
    getSourceTermData();

    // loop over elements
    for (CFuint iElem = 0; iElem < nbrElems; ++iElem, ++cellIdx)
    {
      // build the GeometricEntity
      geoData.idx = cellIdx;
      m_cell = geoBuilder->buildGE();

      // get the states
      m_cellStates = m_cell->getStates();

      // get jacobian determinants at solution points
      m_solPntJacobDets = m_cell->computeGeometricShapeFunctionJacobianDeterminant(*m_solPntsLocalCoords);

      // add the source term if the current cell is parallel updatable
      if ((*m_cellStates)[0]->isParUpdatable())
      {
	m_resUpdates = 0.;
	
        addSourceTerm(m_resUpdates);
	
	updateRHS();
	
	if (m_addJacob && (!getMethodData().freezeJacob() || iter < iterFreeze || interval % getMethodData().getFreezeJacobInterval() == 0))
	{
          m_isPerturbed = true;
            
          if (!m_useAnaJacob)
          {
	    addSrcTermJacob();
          }
          else
          {
            addSrcTermJacobAna(); 
          }
          
          m_isPerturbed = false;
	}
      }

      //release the GeometricEntity
      geoBuilder->releaseGE();
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void StdSourceTerm::getSourceTermData()
{
  // get the local FR data
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();

  // get solution point local coordinates
  m_solPntsLocalCoords = frLocalData[m_iElemType]->getSolPntsLocalCoords();
}

//////////////////////////////////////////////////////////////////////////////

void StdSourceTerm::updateRHS()
{
  // get the datahandle of the rhs
  DataHandle< CFreal > rhs = socket_rhs.getDataHandle();

  // get residual factor
  const CFreal resFactor = getMethodData().getResFactor();
  
  for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
  {
    CFuint resID = m_nbrEqs*( (*m_cellStates)[iState]->getLocalID() );
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      rhs[resID+iVar] += resFactor*m_solPntJacobDets[iState]*m_resUpdates[m_nbrEqs*iState+iVar];  //*m_divContFlx[iState][iVar];
    }
  }

//  // loop over solution points in this cell to add the source term
//  CFuint resID = m_nbrEqs*( (*m_cellStates)[0]->getLocalID() );
//
//  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
//  {
//    // loop over physical variables
//    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq, ++resID)
//    {
//      rhs[resID] += resFactor*m_solPntJacobDets[iSol]*m_resUpdates[m_nbrEqs*iSol+iEq];
//    }
//  }
}

//////////////////////////////////////////////////////////////////////////////

void StdSourceTerm::addSrcTermJacob()
{
  // get residual factor
  const CFreal resFactor = getMethodData().getResFactor();

  // dereference accumulator
  BlockAccumulator& acc = *m_acc;

  // set block row and column indices
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    acc.setRowColIndex(iSol,(*m_cellStates)[iSol]->getLocalID());
  }

  // loop over the states/solpnts in this cell to perturb the states
  for (m_pertSol = 0; m_pertSol < m_nbrSolPnts; ++m_pertSol)
  {
    // dereference state
    State& pertState = *(*m_cellStates)[m_pertSol];

    // loop over the variables in the state
    for (CFuint iEqPert = 0; iEqPert < m_nbrEqs; ++iEqPert)
    {
      // perturb physical variable in state
      m_numJacob->perturb(iEqPert,pertState[iEqPert]);

      m_pertResUpdates = 0.;

      // compute the perturbed residual updates (-divFD+divhFD)
      addSourceTerm(m_pertResUpdates);

      // compute the finite difference derivative of the volume term
      m_numJacob->computeDerivative(m_pertResUpdates,m_resUpdates,m_derivResUpdates);

      // multiply residual update derivatives with residual factor
      m_derivResUpdates *= resFactor;

      // add the derivative of the residual updates to the accumulator
      CFuint resUpdIdx = 0;
      for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol, resUpdIdx += m_nbrEqs)
      {
        for (CFuint iEq=0; iEq < m_nbrEqs; ++iEq)
        {
          m_derivResUpdates[resUpdIdx + iEq] *= m_solPntJacobDets[iSol];
        }
          
        acc.addValues(iSol,m_pertSol,iEqPert,&m_derivResUpdates[resUpdIdx]);
      }

      // restore physical variable in state
      m_numJacob->restore(pertState[iEqPert]);
    }
  }
//   if (m_cell->getID() == 49)
//   {
//   CFLog(VERBOSE,"accVol:\n");
//    acc.printToScreen();
//   }

  if (getMethodData().doComputeJacobian())
  {
    // add the values to the jacobian matrix
    m_lss->getMatrix()->addValues(acc);
  }

  // reset to zero the entries in the block accumulator
  acc.reset();
}

//////////////////////////////////////////////////////////////////////////////

void StdSourceTerm::addSrcTermJacobAna()
{
  // get residual factor
  const CFreal resFactor = getMethodData().getResFactor();

  // dereference accumulator
  BlockAccumulator& acc = *m_acc;

  // set block row and column indices
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    acc.setRowColIndex(iSol,(*m_cellStates)[iSol]->getLocalID());
  }

  // loop over the states/solpnts in this cell to perturb the states
  for (m_pertSol = 0; m_pertSol < m_nbrSolPnts; ++m_pertSol)
  {
    // compute the jacobian to state
    getSToStateJacobian(m_pertSol);

    // loop over the variables in the state
    for (CFuint iEqPert = 0; iEqPert < m_nbrEqs; ++iEqPert)
    {
//      // perturb physical variable in state
//      m_numJacob->perturb(iEqPert,pertState[iEqPert]);
//
//      m_pertResUpdates = 0.;
//
//      // compute the perturbed residual updates (-divFD+divhFD)
//      addSourceTerm(m_pertResUpdates);
//
//      // compute the finite difference derivative of the volume term
//      m_numJacob->computeDerivative(m_pertResUpdates,m_resUpdates,m_derivResUpdates);

      // multiply residual update derivatives with residual factor
      m_stateJacobian[iEqPert] *= resFactor;

      // add the derivative of the residual updates to the accumulator
      m_stateJacobian[iEqPert] *= m_solPntJacobDets[m_pertSol];
          
      acc.addValues(m_pertSol,m_pertSol,iEqPert,&m_stateJacobian[iEqPert][0]);

//      // restore physical variable in state
//      m_numJacob->restore(pertState[iEqPert]);
    }
  }
  
  if (isGradDependent())
  {
//    for (m_pertSol = 0; m_pertSol < m_nbrSolPnts; ++m_pertSol)
//    {
//      // compute the jacobian to gradient
//      getSToGradJacobian(m_pertSol);
//      
//      // loop over the variables in the state
//      for (CFuint iEqPert = 0; iEqPert < m_nbrEqs; ++iEqPert)
//      {
//        // multiply residual update derivatives with residual factor
//        m_stateJacobian[iEqPert] *= resFactor;
//
//        // add the derivative of the residual updates to the accumulator
//        m_stateJacobian[iEqPert] *= m_solPntJacobDets[m_pertSol];
//          
//        //acc.addValues(m_pertSol,m_pertSol,iEqPert,&m_stateJacobian[iEqPert][0]);
//      }  
//    } 
  }
  
//   if (m_cell->getID() == 49)
//   {
//   CFLog(VERBOSE,"accVol:\n");
//    acc.printToScreen();
//   }

  if (getMethodData().doComputeJacobian())
  {
    // add the values to the jacobian matrix
    m_lss->getMatrix()->addValues(acc);
  }

  // reset to zero the entries in the block accumulator
  acc.reset();
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > >
    StdSourceTerm::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result;

  result.push_back(&socket_rhs);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD
