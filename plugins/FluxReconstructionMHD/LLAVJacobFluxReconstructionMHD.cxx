// Copyright (C) 2022 KU Leuven CmPA, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/MethodCommandProvider.hh"

#include "Framework/CFSide.hh"
#include "Framework/MeshData.hh"
#include "Framework/BaseTerm.hh"

#include "MathTools/MathFunctions.hh"

#include "FluxReconstructionMHD/LLAVJacobFluxReconstructionMHD.hh"
#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"

#include "MHD/MHD3DProjectionVarSet.hh"

#include "FluxReconstructionMHD/FluxReconstructionMHD.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::MHD;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {
    
//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< LLAVJacobFluxReconstructionMHD,
		       FluxReconstructionSolverData,
		       FluxReconstructionMHDModule >
LLAVJacobFluxReconstructionMHDFluxReconstructionProvider("LLAVJacobMHD");

//////////////////////////////////////////////////////////////////////////////
  
LLAVJacobFluxReconstructionMHD::LLAVJacobFluxReconstructionMHD(const std::string& name) :
  LLAVJacobFluxReconstruction(name),
  m_varSet(CFNULL)
  {
  }

//////////////////////////////////////////////////////////////////////////////

void LLAVJacobFluxReconstructionMHD::configure ( Config::ConfigArgs& args )
{
  LLAVJacobFluxReconstruction::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSource > >
  LLAVJacobFluxReconstructionMHD::providesSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSource > > result;
  //result.push_back(&socket_artVisc);
  //result.push_back(&socket_monPhysVar);
  //result.push_back(&socket_smoothness);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void LLAVJacobFluxReconstructionMHD::setCellData(const CFuint side)
{ 
  m_cellNodes = m_cells[side]->getNodes();
  
  //DataHandle< CFreal > artVisc = socket_artVisc.getDataHandle();
  
  // loop over flx pnts to extrapolate the states to the flux points
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {   
//     m_solEpsilons[iSol] = m_cellEpsilons[m_cell->getID()];
    
    // reset the states in the flx pnts
    m_solEpsilons[side][iSol] = 0.0;

    // loop over the sol pnts to compute the states and grads in the flx pnts
    for (CFuint iNode = 0; iNode < m_nbrCornerNodes; ++iNode)
    {
      // get node local index
      //const CFuint nodeIdx = (*m_cellNodesConn)(m_elemIdx,iNode);
      
      const CFuint nodeIdx = (*m_cellNodes)[iNode]->getLocalID();
      
      m_solEpsilons[side][iSol] += m_nodePolyValsAtSolPnts[iSol][iNode]*m_nodeEpsilons[nodeIdx]/m_nbNodeNeighbors[nodeIdx];
    }
    
    //artVisc[(((*m_states[side])[iSol]))->getLocalID()] = m_solEpsilons[side][iSol];
  }
}

//////////////////////////////////////////////////////////////////////////////

void LLAVJacobFluxReconstructionMHD::computeSmoothness()
{ 
  CFreal sNum = 0.0;
  
  CFreal sDenom = 0.0;
  
  // get datahandle
  //DataHandle< CFreal > monPhysVar = socket_monPhysVar.getDataHandle();
  
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    CFreal stateP = 0.0;
    CFreal diffStatesPPMinOne = 0.0;
    
    if (m_monitoredPhysVar < m_pData.size())
    {
      RealVector statePMinOne = *((*m_cellStates)[iSol]->getData()) - m_statesPMinOne[iSol];
            
      m_varSet->computePhysicalData(*((*m_cellStates)[iSol]),m_pData);
      m_varSet->computePhysicalData(m_statesPMinOne[iSol],m_pData2);
      
      stateP = m_pData[m_monitoredPhysVar];
      diffStatesPPMinOne = stateP - m_pData2[m_monitoredPhysVar];
    
      //monPhysVar[(((*m_cellStates)[iSol]))->getLocalID()] = stateP;
    }
    else
    {
      stateP = (*((*m_cellStates)[iSol]))[m_monitoredVar];
      diffStatesPPMinOne = stateP - m_statesPMinOne[iSol][m_monitoredVar];
    }
    
    sNum += diffStatesPPMinOne*diffStatesPPMinOne;
    sDenom += stateP*stateP;
  }
  if (sNum <= MathTools::MathConsts::CFrealEps() || sDenom <= MathTools::MathConsts::CFrealEps())
  {
    m_s = -100.0;
  }
  else
  {
    m_s = log10(sNum/sDenom);
  }
  
  // get datahandle
  //DataHandle< CFreal > smoothness = socket_smoothness.getDataHandle();
  
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    //smoothness[(((*m_cellStates)[iSol]))->getLocalID()] = m_s;
  }
  
  if (m_s > m_Smax)
  {
      m_Smax = m_s;
  }
}

//////////////////////////////////////////////////////////////////////////////

void LLAVJacobFluxReconstructionMHD::computeSmoothness(const CFuint side)
{ 
  CFreal sNum = 0.0;
  
  CFreal sDenom = 0.0;
  
  // get datahandle
  //DataHandle< CFreal > monPhysVar = socket_monPhysVar.getDataHandle();
  
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    CFreal stateP = 0.0;
    CFreal diffStatesPPMinOne = 0.0;
    
    if (m_monitoredPhysVar < m_pData.size())
    {
      RealVector statePMinOne = *((*m_states[side])[iSol]->getData()) - m_statesPMinOne[iSol];
            
      m_varSet->computePhysicalData(*((*m_states[side])[iSol]),m_pData);
      m_varSet->computePhysicalData(m_statesPMinOne[iSol],m_pData2);
      
      stateP = m_pData[m_monitoredPhysVar];
      diffStatesPPMinOne = stateP - m_pData2[m_monitoredPhysVar];
    
      //monPhysVar[(((*m_states[side])[iSol]))->getLocalID()] = stateP;
    }
    else
    {
      stateP = (*((*m_states[side])[iSol]))[m_monitoredVar];
      diffStatesPPMinOne = stateP - m_statesPMinOne[iSol][m_monitoredVar];
    }
    
    sNum += diffStatesPPMinOne*diffStatesPPMinOne;
    sDenom += stateP*stateP;
  }
  if (sNum <= MathTools::MathConsts::CFrealEps() || sDenom <= MathTools::MathConsts::CFrealEps())
  {
    m_s = -100.0;
  }
  else
  {
    m_s = log10(sNum/sDenom);
  }
}

//////////////////////////////////////////////////////////////////////////////

void LLAVJacobFluxReconstructionMHD::setup()
{
  CFAUTOTRACE;
  
  // setup parent class
  DiffRHSJacobFluxReconstruction::setup();

  // get the update varset
  m_updateVarSet = getMethodData().getUpdateVar();
  
  // get cell-node connectivity
  m_cellNodesConn = MeshDataStack::getActive()->getConnectivity("cellNodes_InnerCells");
  
  // get the local spectral FD data
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
  cf_assert(frLocalData.size() > 0);
  // for now, there should be only one type of element
  cf_assert(frLocalData.size() == 1);
  
  const CFPolyOrder::Type order = frLocalData[0]->getPolyOrder();
  
  m_subcellRes = frLocalData[0]->getSubcellResolution();
  
  m_order = static_cast<CFuint>(order);
  
  // get the coefs for extrapolation of the node artificial viscosities to the flx pnts
  m_nodePolyValsAtFlxPnts = frLocalData[0]->getNodePolyValsAtPnt(*(frLocalData[0]->getFlxPntsLocalCoords()));
  
  // get the coefs for extrapolation of the node artificial viscosities to the sol pnts
  m_nodePolyValsAtSolPnts = frLocalData[0]->getNodePolyValsAtPnt(*(frLocalData[0]->getSolPntsLocalCoords()));
  
  // number of cell corner nodes
  /// @note in the future, hanging nodes should be taken into account here
  m_nbrCornerNodes = frLocalData[0]->getNbrCornerNodes();
  
  // get the number of nodes in the mesh
  const CFuint nbrNodes = MeshDataStack::getActive()->getNbNodes();
  
  // get the elementTypeData
  SafePtr< vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();
  
  // get the number of cells in the mesh
  const CFuint nbrCells = (*elemType)[0].getEndIdx();
  
  // get datahandle
  //DataHandle< CFreal > artVisc = socket_artVisc.getDataHandle();
  //DataHandle< CFreal > monPhysVar = socket_monPhysVar.getDataHandle();
  //DataHandle< CFreal > smoothness = socket_smoothness.getDataHandle();
  
  const CFuint nbStates = nbrCells*m_nbrSolPnts;

  // resize socket
  //artVisc.resize(nbStates);
  //monPhysVar.resize(nbStates);
  //smoothness.resize(nbStates);
  
  m_nodeEpsilons.resize(nbrNodes);
  m_nbNodeNeighbors.resize(nbrNodes);
  m_cellEpsilons.resize(nbrCells);
  m_solEpsilons.resize(2);
  m_solEpsilons[LEFT].resize(m_nbrSolPnts);
  m_solEpsilons[RIGHT].resize(m_nbrSolPnts);
  m_epsilonLR.resize(2);
  m_epsilonLR[LEFT].resize(m_nbrFaceFlxPnts);
  m_epsilonLR[RIGHT].resize(m_nbrFaceFlxPnts);
  m_unitNormalFlxPnts2.resize(m_nbrFaceFlxPnts);
  m_tempSolPntVec.resize(m_nbrSolPnts);
  m_tempSolPntVec2.resize(m_nbrSolPnts);
  
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    RealVector temp(m_nbrEqs);
    temp = 0.0;
    m_statesPMinOne.push_back(temp);
  }
  
  SafePtr<RealMatrix> vdm = frLocalData[0]->getVandermondeMatrix();
  
  SafePtr<RealMatrix> vdmInv = frLocalData[0]->getVandermondeMatrixInv();
  
  RealMatrix temp(m_nbrSolPnts,m_nbrSolPnts);
  temp = 0.0;

  if (elemShape == CFGeoShape::TRIAG){
    for (CFuint idx = 0; idx < (m_order)*(m_order+1)/2; ++idx)
    {
      temp(idx,idx) = 1.0;
    }
  }
  else if(elemShape == CFGeoShape::QUAD){
    for (CFuint idx = 0; idx < (m_order)*(m_order); ++idx)
    {
      temp(idx,idx) = 1.0;
    }
  }  
  else if (elemShape == CFGeoShape::TETRA){
    for (CFuint idx = 0; idx < (m_order)*(m_order+1)*(m_order+2)/6; ++idx)
    {
      temp(idx,idx) = 1.0;
    }
  }
  else if (elemShape == CFGeoShape::PRISM){
    for (CFuint idx = 0; idx < (m_order)*(m_order)*(m_order+1)/2; ++idx)
    {
      temp(idx,idx) = 1.0;
    }
  }  
  else{
    for (CFuint idx = 0; idx < (m_order)*(m_order)*(m_order); ++idx)
    {
      temp(idx,idx) = 1.0;
    }
  }
  
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_unitNormalFlxPnts2[iFlx].resize(m_dim);
  }
  
  m_transformationMatrix.resize(m_nbrSolPnts,m_nbrSolPnts);
  
  m_transformationMatrix = (*vdm)*temp*(*vdmInv);
  
  //m_s0 = -m_s0*log10(static_cast<CFreal>(m_order));
  
  m_Smax = m_s0 + m_kappa;
  
  m_SmaxGlobal = m_Smax;
  
  m_nbNodeNeighbors = 0.0;
  
  m_flagComputeNbNghb = true;
  
  cf_assert(m_monitoredVar < m_nbrEqs);
  
  // get damping coeff
  m_dampCoeff = getMethodData().getDiffDampCoefficient();

  // get 3D varset
  m_varSet = getMethodData().getUpdateVar().d_castTo< MHD3DProjectionVarSet >();

  if (m_varSet.isNull())
  {
    throw Common::ShouldNotBeHereException (FromHere(),"Update variable set is not MHD3DProjectionVarSet in LLAVJacobFRMHD!\n");
  }

  // resize the physical data for internal and ghost solution points
  m_varSet->getModel()->resizePhysicalData(m_pData);
  m_varSet->getModel()->resizePhysicalData(m_pData2);
  
}

//////////////////////////////////////////////////////////////////////////////

void LLAVJacobFluxReconstructionMHD::setFaceData(CFuint faceID)
{
  LLAVJacobFluxReconstruction::setFaceData(faceID);
  
  if (getMethodData().getUpdateVarStr() != "Puvt" && getMethodData().hasDiffTerm())
  {
    // get the gradients datahandle
    DataHandle< vector< RealVector > > gradientsAV = socket_gradientsAV.getDataHandle();

    for (CFuint iSide = 0; iSide < 2; ++iSide)
    {
      for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
      {
        const CFuint stateID = (*(m_states[iSide]))[iState]->getLocalID();
        m_cellGrads[iSide][iState] = &gradientsAV[stateID];
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void LLAVJacobFluxReconstructionMHD::computeInterfaceFlxCorrection()
{   
  // Loop over the flux points to calculate FI
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  { 
    const CFreal epsilon = 0.5*(m_epsilonLR[LEFT][iFlxPnt]+m_epsilonLR[RIGHT][iFlxPnt]);
    
    // compute the average sol and grad to use the BR2 scheme
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      *(m_avgGrad[iVar]) = (*(m_cellGradFlxPnt[LEFT][iFlxPnt][iVar]) + *(m_cellGradFlxPnt[RIGHT][iFlxPnt][iVar]))/2.0;
    }
    
    // damping factor
    const CFreal dampFactor = m_dampCoeff*m_faceInvCharLengths[iFlxPnt];

    // compute averaged (damped) gradients
    for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
    {
      // compute damping term
      const RealVector dGradVarXNormal = ((*m_cellStatesFlxPnt[LEFT][iFlxPnt])[iGrad] - (*m_cellStatesFlxPnt[RIGHT][iFlxPnt])[iGrad])*m_unitNormalFlxPnts[iFlxPnt];
      *m_avgGrad[iGrad] -= dampFactor*dGradVarXNormal;
    }
    
    m_flxPntRiemannFlux[iFlxPnt] = 0.0;
    
    for (CFuint iDim = 0; iDim < m_dim; ++iDim)
    {
      for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
      {
        m_flxPntRiemannFlux[iFlxPnt][iVar] += epsilon*((*(m_avgGrad[iVar]))[iDim])*m_unitNormalFlxPnts[iFlxPnt][iDim];
      }
    }
     
    // compute FI in the mapped coord frame
    m_cellFlx[LEFT][iFlxPnt] = (m_flxPntRiemannFlux[iFlxPnt])*m_faceJacobVecSizeFlxPnts[iFlxPnt][LEFT];
    m_cellFlx[RIGHT][iFlxPnt] = (m_flxPntRiemannFlux[iFlxPnt])*m_faceJacobVecSizeFlxPnts[iFlxPnt][RIGHT];
  }
}

//////////////////////////////////////////////////////////////////////////////

void LLAVJacobFluxReconstructionMHD::unsetup()
{
  CFAUTOTRACE;
  
  // unsetup parent class
  LLAVJacobFluxReconstruction::unsetup();
}


//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

