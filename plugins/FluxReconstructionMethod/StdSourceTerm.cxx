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
  m_stateJacobian(),
  m_gradientStateJacobian()
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
  
  m_dim = PhysicalModelStack::getActive()->getDim();
  
  m_pertResUpdates.resize(resSize);
  m_resUpdates.resize(resSize);
  m_derivResUpdates.resize(resSize);
  
  m_stateJacobian.resize(m_nbrEqs);
  
  for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
  {
    m_stateJacobian[iEq].resize(m_nbrEqs);
  }
  
  m_gradientStateJacobian.resize(m_nbrSolPnts);
  
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  { 
    m_gradientStateJacobian[iSol].resize(2);
    m_gradientStateJacobian[iSol].resize(2);
    
    for (CFuint iSide = 0; iSide < 2; ++iSide)
    {
      m_gradientStateJacobian[iSol][iSide].resize(m_nbrSolPnts);
      m_gradientStateJacobian[iSol][iSide].resize(m_nbrSolPnts);
      
      for (CFuint jSol = 0; jSol < m_nbrSolPnts; ++jSol)
      {
        m_gradientStateJacobian[iSol][iSide][jSol].resize(m_dim);
        m_gradientStateJacobian[iSol][iSide][jSol].resize(m_dim);
      }
    }
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
//   if (m_cell->getID() == 1)
//   {
//   CFLog(VERBOSE,"accST:\n");
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
  const CFreal resFactor = -getMethodData().getResFactor();

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
      
      //if (m_cell->getID() == 1) CFLog(INFO,"solJacob: " << m_solPntJacobDets[m_pertSol] << "\n");
          
      acc.addValues(m_pertSol,m_pertSol,iEqPert,&m_stateJacobian[iEqPert][0]);

//      // restore physical variable in state
//      m_numJacob->restore(pertState[iEqPert]);
    }
  }
  
  if (isGradDependent())
  {
//    for (m_pertSol = 0; m_pertSol < m_nbrSolPnts; ++m_pertSol)
//    {
//      computeGradToStateJacobianAna();
//              
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
  
//   if (m_cell->getID() == 1)
//   {
//   CFLog(VERBOSE, "accST:\n");
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

void StdSourceTerm::computeGradToStateJacobianAna()
{
  CFLog(VERBOSE, "computeGradToStateJacobianAna\n");
  
  // reset the jacobian
  for (m_pertSol = 0; m_pertSol < m_nbrSolPnts; ++m_pertSol)
  {
    for (CFuint jSol = 0; jSol < m_nbrSolPnts; ++jSol)
    {
      for (CFuint iDim = 0; iDim < m_dim; ++iDim)
      { 
        m_gradientStateJacobian[m_pertSol][LEFT][jSol][iDim] = 0.0;
        m_gradientStateJacobian[m_pertSol][RIGHT][jSol][iDim] = 0.0;
      }
    }
  }

//  
//  // get the face flux point normals
//  DataHandle< CFreal > flxPntNormals = socket_flxPntNormals.getDataHandle();
//    
//  for (m_pertSide = 0; m_pertSide < 2; ++m_pertSide)
//  {
//    // variable for the other side
//    const CFuint iOtherSide = m_pertSide == LEFT ? RIGHT : LEFT;
//    
//    // loop over the sol pnts of which the gradients will be derived
//    for (m_pertSol = 0; m_pertSol < m_nbrSolPnts; ++m_pertSol)
//    {
//      // loop over the sol pnts to which to derive
//      for (CFuint jSol = 0; jSol < m_nbrSolSolDep; ++jSol)
//      {
//        const CFuint jSolIdx = (*m_solSolDep)[m_pertSol][jSol];  
//        
//        for (CFuint iDim = 0; iDim < m_dim; ++iDim)
//        {                    
//          for (CFuint jDim = 0; jDim < m_dim; ++jDim)
//          {
//            m_gradientStateJacobian[m_pertSide][jSolIdx][m_pertSide][m_pertSol][iDim] += m_neighbCellFluxProjVects[m_pertSide][jDim][m_pertSol][iDim] * (*m_solPolyDerivAtSolPnts)[jSolIdx][jDim][m_pertSol];
//          }
//        }
//      }
//    }
//    
//    for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
//    {
//      // local flux point indices in the left and right cell
//      const CFuint flxPntIdxThis = (*m_faceFlxPntConnPerOrient)[m_orient][m_pertSide][iFlx];
//      const CFuint flxPntIdxOther = (*m_faceFlxPntConnPerOrient)[m_orient][iOtherSide][iFlx];
//      
//      for (m_pertSol = 0; m_pertSol < m_nbrSolDep; ++m_pertSol)
//      {
//        const CFuint pertSolIdx = (*m_flxSolDep)[flxPntIdxThis][m_pertSol]; 
//    
//        for (CFuint jSol = 0; jSol < m_nbrSolDep; ++jSol)
//        {
//          const CFuint jSolIdx = (*m_flxSolDep)[flxPntIdxThis][jSol];  
//          const CFuint jSolIdxOther = (*m_flxSolDep)[flxPntIdxOther][jSol];  
//        
//          for (CFuint iDim = 0; iDim < m_dim; ++iDim)
//          {
//            m_gradientStateJacobian[m_pertSide][jSolIdx][m_pertSide][pertSolIdx][iDim] -= 0.5 * m_corrFctDiv[jSolIdx][flxPntIdxThis] * (*m_solPolyValsAtFlxPnts)[flxPntIdxThis][pertSolIdx] *
//                    (*m_faceMappedCoordDir)[m_orient][m_pertSide]*m_faceJacobVecs[iFlx][iDim];
//         
//            m_gradientStateJacobian[iOtherSide][jSolIdxOther][m_pertSide][pertSolIdx][iDim] += 0.5 * m_corrFctDiv[jSolIdxOther][flxPntIdxOther] * (*m_solPolyValsAtFlxPnts)[flxPntIdxThis][pertSolIdx] *
//                    (*m_faceMappedCoordDir)[m_orient][iOtherSide]*m_faceJacobVecs[iFlx][iDim];
//          }
//        }   
//      }
//    }
//  
//    // Add the contribution of the correction to the gradients for each face
//    // compute other face contributions to the gradients
//    const CFuint nbrOtherFaces = m_otherFaceLocalIdxs[m_pertSide].size();
//  
//    for (CFuint iFace = 0; iFace < nbrOtherFaces; ++iFace)
//    {
//      // get local face index
//      const CFuint faceIdx = m_otherFaceLocalIdxs[m_pertSide][iFace];
//      
//      const CFuint faceID = (*m_faces[m_pertSide])[faceIdx]->getID();
//      
//      if ((*m_isFaceOnBoundary[m_pertSide])[faceIdx])
//      {  
//        for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
//        {
//          const CFuint currFlxIdx = (*m_faceFlxPntConn)[faceIdx][iFlxPnt];
//          
//          for (m_pertSol = 0; m_pertSol < m_nbrSolDep; ++m_pertSol)
//          {
//            const CFuint pertSolIdx = (*m_flxSolDep)[currFlxIdx][m_pertSol];        
//        
//            for (CFuint jSol = 0; jSol < m_nbrSolDep; ++jSol)
//            {
//              const CFuint jSolIdx = (*m_flxSolDep)[currFlxIdx][jSol]; 
//    
//              for (CFuint iDim = 0; iDim < m_dim; ++iDim)
//              {
//                m_gradientStateJacobian[m_pertSide][jSolIdx][m_pertSide][pertSolIdx][iDim] -= 0.5 * m_corrFctDiv[jSolIdx][currFlxIdx] * (*m_solPolyValsAtFlxPnts)[currFlxIdx][pertSolIdx] *
//                          (*m_faceLocalDir)[faceIdx]*flxPntNormals[faceID*m_nbrFaceFlxPnts*m_dim+iFlxPnt*m_dim+iDim];
//              }
//            }
//          }
//        }
//      }
//      else
//      {
//        // Get orientation of face
//        const CFuint orient = (*m_faceOrients[m_pertSide])[faceIdx];
//        
//        // cell side with respect to this face
//        const CFuint cellSide = (*m_currCellSide[m_pertSide])[faceIdx];
//      
//        for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
//        {
//          const CFuint currFlxIdx = (*m_faceFlxPntConnPerOrient)[orient][cellSide][iFlxPnt];
//          
//          for (m_pertSol = 0; m_pertSol < m_nbrSolDep; ++m_pertSol)
//          {
//            const CFuint pertSolIdx = (*m_flxSolDep)[currFlxIdx][m_pertSol];        
//        
//            for (CFuint jSol = 0; jSol < m_nbrSolDep; ++jSol)
//            {
//              const CFuint jSolIdx = (*m_flxSolDep)[currFlxIdx][jSol]; 
//    
//              for (CFuint iDim = 0; iDim < m_dim; ++iDim)
//              {
//                m_gradientStateJacobian[m_pertSide][jSolIdx][m_pertSide][pertSolIdx][iDim] -= 0.5 * m_corrFctDiv[jSolIdx][currFlxIdx] * (*m_solPolyValsAtFlxPnts)[currFlxIdx][pertSolIdx] *
//                          ((*m_faceMappedCoordDir)[orient][cellSide])*flxPntNormals[faceID*m_nbrFaceFlxPnts*m_dim+iFlxPnt*m_dim+iDim];
//              }
//            }
//          }
//        } 
//      }
//    }
//  }
//  
//  for (m_pertSide = 0; m_pertSide < 2; ++m_pertSide)
//  { 
//    for (m_pertSol = 0; m_pertSol < m_nbrSolPnts; ++m_pertSol)
//    {
//      
//      
//      for (CFuint jSol = 0; jSol < m_nbrSolPnts; ++jSol)
//      {
//        const CFreal invJacobL = 1./m_solJacobDet[LEFT][jSol];
//        const CFreal invJacobR = 1./m_solJacobDet[RIGHT][jSol];
//          
//        for (CFuint iDim = 0; iDim < m_dim; ++iDim)
//        { 
//          m_gradientStateJacobian[LEFT][jSol][m_pertSide][m_pertSol][iDim] *= invJacobL;
//          m_gradientStateJacobian[RIGHT][jSol][m_pertSide][m_pertSol][iDim] *= invJacobR;
//        }
////        if (m_cells[LEFT]->getID()==1) CFLog(INFO,"side: " << m_pertSide << ", sol: " << m_pertSol << ", to side: 0, sol: " << jSol << ": " << m_gradientStateJacobian[LEFT][jSol][m_pertSide][m_pertSol] <<"\n");
////      if (m_cells[LEFT]->getID()==1) CFLog(INFO,"side: " << m_pertSide << ", sol: " << m_pertSol << ", to side: 1, sol: " << jSol << ": " << m_gradientStateJacobian[RIGHT][jSol][m_pertSide][m_pertSol] <<"\n");
//      }
//    }
//  }
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD
