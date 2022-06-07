// Copyright (C) 2019 KU Leuven, von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/MethodCommandProvider.hh"

#include "Framework/CFSide.hh"
#include "Framework/MeshData.hh"
#include "Framework/BaseTerm.hh"

#include "MathTools/MathFunctions.hh"

#include "FluxReconstructionTurb/ConvDiffLLAVFluxReconstructionTurb.hh"
#include "FluxReconstructionNavierStokes/FluxReconstructionNavierStokes.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"
#include "NavierStokes/NavierStokesVarSet.hh"

#include "NavierStokes/EulerTerm.hh"
#include "NavierStokes/EulerVarSet.hh"

#include "NavierStokes/Euler2DVarSet.hh"

#include "KOmega/NavierStokesKLogOmegaVarSetTypes.hh"


//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::Physics::KOmega;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< ConvDiffLLAVFluxReconstructionTurb,
		       FluxReconstructionSolverData,
		       FluxReconstructionNavierStokesModule >
convDiffLLAVRHSTurbFluxReconstructionProvider("ConvDiffLLAVRHSTurb");
  
//////////////////////////////////////////////////////////////////////////////
  
ConvDiffLLAVFluxReconstructionTurb::ConvDiffLLAVFluxReconstructionTurb(const std::string& name) :
  ConvDiffLLAVFluxReconstructionNS(name),
  socket_wallDistance("wallDistance"),
  m_closestSolToFlxIdx(CFNULL),
  m_pData(),
  m_pData2()
{
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffLLAVFluxReconstructionTurb::configure ( Config::ConfigArgs& args )
{
  ConvDiffLLAVFluxReconstructionNS::configure(args);
} 

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > >
    ConvDiffLLAVFluxReconstructionTurb::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result = ConvDiffLLAVFluxReconstructionNS::needsSockets();

  result.push_back(&socket_wallDistance);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffLLAVFluxReconstructionTurb::computeInterfaceFlxCorrection()
{ 
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_tempStatesL[iFlx] = m_cellStatesFlxPnt[LEFT][iFlx]->getData();
    m_tempStatesR[iFlx] = m_cellStatesFlxPnt[RIGHT][iFlx]->getData();
  }
  
  m_diffusiveVarSetNS->setGradientVars(m_tempStatesL,m_tempGradTermL,m_nbrFaceFlxPnts);
  m_diffusiveVarSetNS->setGradientVars(m_tempStatesR,m_tempGradTermR,m_nbrFaceFlxPnts);
  
  // Get the wall distance
  DataHandle< CFreal > wallDist = socket_wallDistance.getDataHandle();
  
  //SafePtr< NavierStokes2DKLogOmega > navierStokesVarSet = m_diffusiveVarSet.d_castTo< NavierStokes2DKLogOmega >();
    
  // Loop over the flux points to calculate FI
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  { 
    // local flux point indices in the left and right cell
    const CFuint flxPntIdxL = (*m_faceFlxPntConnPerOrient)[m_orient][LEFT][iFlxPnt];
    const CFuint flxPntIdxR = (*m_faceFlxPntConnPerOrient)[m_orient][RIGHT][iFlxPnt];
    
    const CFuint closestSolIdxL = (*m_closestSolToFlxIdx)[flxPntIdxL];
    const CFuint closestSolIdxR = (*m_closestSolToFlxIdx)[flxPntIdxR];
    
    const CFuint stateIDL = (*(m_states[LEFT]))[closestSolIdxL]->getLocalID();
    const CFuint stateIDR = (*(m_states[RIGHT]))[closestSolIdxR]->getLocalID();
    
    const CFreal flxPntWallDist = 0.5*(wallDist[stateIDL]+wallDist[stateIDR]);
    
    // Set the wall distance before computing the turbulent viscosity
    if (m_dim == 2)
    {
      m_navierStokesVarSetTurb->setWallDistance(flxPntWallDist);
    }
    else
    {
      m_navierStokesVarSetTurb3D->setWallDistance(flxPntWallDist); 
    }
    
    // compute the average sol and grad to use the BR2 scheme
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      *(m_avgGrad[iVar]) = (*(m_cellGradFlxPnt[LEFT][iFlxPnt][iVar]) + *(m_cellGradFlxPnt[RIGHT][iFlxPnt][iVar]))/2.0;
      *(m_avgGradAV[iVar]) = (*(m_cellGradFlxPntAV[LEFT][iFlxPnt][iVar]) + *(m_cellGradFlxPntAV[RIGHT][iFlxPnt][iVar]))/2.0;
       
      m_avgSol[iVar] = ((*(m_cellStatesFlxPnt[LEFT][iFlxPnt]))[iVar] + (*(m_cellStatesFlxPnt[RIGHT][iFlxPnt]))[iVar])/2.0; 
    }

    // damping factor
    const CFreal dampFactor = m_dampCoeffDiff*m_faceInvCharLengths[iFlxPnt];
    
    m_tempSolVarState = static_cast<RealVector&>(*m_updateToSolutionVecTrans->transform(m_cellStatesFlxPnt[LEFT][iFlxPnt]));
    m_tempSolVarState2 = static_cast<RealVector&>(*m_updateToSolutionVecTrans->transform(m_cellStatesFlxPnt[RIGHT][iFlxPnt]));

    // compute averaged (damped) gradients
    for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
    {
      // compute damping term
      const RealVector dGradVarXNormal = (m_tempGradTermL(iGrad,iFlxPnt) - m_tempGradTermR(iGrad,iFlxPnt))*m_unitNormalFlxPnts[iFlxPnt];
      *m_avgGrad[iGrad] -= dampFactor*dGradVarXNormal;
      
      // compute damping term for LLAV
      const RealVector dGradVarXNormalAV = (m_tempSolVarState[iGrad] - m_tempSolVarState2[iGrad])*m_unitNormalFlxPnts[iFlxPnt];
      *m_avgGradAV[iGrad] -= dampFactor*dGradVarXNormalAV;
    }
    
    prepareFluxComputation();
     
    // compute the diff flux
    computeFlux(m_avgSol,m_avgGrad,m_unitNormalFlxPnts[iFlxPnt],0,m_flxPntRiemannFlux[iFlxPnt]);
    
    // compute the convective riemann flux
    m_flxPntRiemannFlux[iFlxPnt] -= m_riemannFluxComputer->computeFlux(*(m_cellStatesFlxPnt[LEFT][iFlxPnt]),
									*(m_cellStatesFlxPnt[RIGHT][iFlxPnt]),
									m_unitNormalFlxPnts[iFlxPnt]);
    
    // compute artificial part
    // get epsilon
    const CFreal epsilon = 0.5*(m_epsilonLR[LEFT][iFlxPnt]+m_epsilonLR[RIGHT][iFlxPnt]);
    
    for (CFuint iDim = 0; iDim < m_dim; ++iDim)
    {
      for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
      {
        m_flxPntRiemannFlux[iFlxPnt][iVar] += epsilon*((*(m_avgGradAV[iVar]))[iDim])*m_unitNormalFlxPnts[iFlxPnt][iDim];
        
//        if (m_cells[LEFT]->getID() == 0) printf("second flx: %d, var: %d, dim: %d, eps: %e, grad: %e, n: %e\n",iFlxPnt,iVar,iDim,epsilon,((*(m_avgGradAV[iVar]))[iDim]),m_unitNormalFlxPnts[iFlxPnt][iDim]*m_faceJacobVecSizeFlxPnts[iFlxPnt][LEFT]);
//      if (m_cells[RIGHT]->getID() == 0) printf("second flx: %d, var: %d, dim: %d, eps: %e, grad: %e, n: %e\n",iFlxPnt,iVar,iDim,epsilon,((*(m_avgGradAV[iVar]))[iDim]),m_unitNormalFlxPnts[iFlxPnt][iDim]*m_faceJacobVecSizeFlxPnts[iFlxPnt][RIGHT]);
      }
    }
     
    // compute FI in the mapped coord frame
    m_cellFlx[LEFT][iFlxPnt] = (m_flxPntRiemannFlux[iFlxPnt])*m_faceJacobVecSizeFlxPnts[iFlxPnt][LEFT];
    m_cellFlx[RIGHT][iFlxPnt] = (m_flxPntRiemannFlux[iFlxPnt])*m_faceJacobVecSizeFlxPnts[iFlxPnt][RIGHT];
  }
  
  // Reset wall distance in case this is the last face: afterwards diff BC is done which needs the wall dist to be zero!!
  if (m_dim == 2)
  {
    m_navierStokesVarSetTurb->setWallDistance(0.0);
  }
  else
  {
    m_navierStokesVarSetTurb3D->setWallDistance(0.0); 
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffLLAVFluxReconstructionTurb::computeUnpertCellDiffResiduals(const CFuint side)
{ 
  //SafePtr< NavierStokes2DKLogOmega > navierStokesVarSet = m_diffusiveVarSet.d_castTo< NavierStokes2DKLogOmega >();
  
  // Get the wall distance
  DataHandle< CFreal > wallDist = socket_wallDistance.getDataHandle();
  
  // get the sol pnt normals
  DataHandle< CFreal > solPntNormals = socket_solPntNormals.getDataHandle();
  
  // store the sol pnt normals
  for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
  {
    const CFuint solID = (*(m_states[side]))[iState]->getLocalID();
      
    for (CFuint iDim = 0; iDim < m_dim; ++iDim)
    {
      for (CFuint jDim = 0; jDim < m_dim; ++jDim)
      {
        m_cellFluxProjVects[iDim][iState][jDim] = solPntNormals[solID*m_dim*m_dim+iDim*m_dim+jDim];
      }
    }
  }
  
  m_cellNodes = m_cells[side]->getNodes();
  
  // get datahandle
  DataHandle< CFreal > artVisc = socket_artVisc.getDataHandle();
  
  // loop over flx pnts to extrapolate the states to the flux points
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {   
    // reset the states in the flx pnts
    m_solEpsilons[side][iSol] = 0.0;

    // loop over the sol pnts to compute the states and grads in the flx pnts
    for (CFuint iNode = 0; iNode < m_nbrCornerNodes; ++iNode)
    {
      // get node local index
      const CFuint nodeIdx = (*m_cellNodes)[iNode]->getLocalID();
      
      m_solEpsilons[side][iSol] += m_nodePolyValsAtSolPnts[iSol][iNode]*m_nodeEpsilons[nodeIdx]/m_nbNodeNeighbors[nodeIdx];
    }
        
    artVisc[(((*(m_states[side]))[iSol]))->getLocalID()] = m_solEpsilons[side][iSol];
    
    //if (m_cells[side]->getID() == 0) CFLog(INFO, "eps: " << m_solEpsilons[side][iSol] << "\n");
  }
  
  // reset the extrapolated fluxes
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrTotalFlxPnts; ++iFlxPnt)
  {
    m_extrapolatedFluxes[iFlxPnt] = 0.0;
  }

  // Loop over solution points to calculate the discontinuous flux.
  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
  { 
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      *(m_tempGrad[iVar]) = (*(m_cellGrads[side][iSolPnt]))[iVar];
      *(m_avgGradAV[iVar]) = (*(m_cellGradsAV[side][iSolPnt]))[iVar];
    }
    
    const CFuint stateID = (*(m_states[side]))[iSolPnt]->getLocalID();
    
    // Set the wall distance before computing the turbulent viscosity
    if (m_dim == 2)
    {
      m_navierStokesVarSetTurb->setWallDistance(wallDist[stateID]);
    }
    else
    {
      m_navierStokesVarSetTurb3D->setWallDistance(wallDist[stateID]); 
    }
    
    m_updateVarSet->computePhysicalData(*((*(m_states[side]))[iSolPnt]), m_pData); 

    m_avgSol = *((*(m_states[side]))[iSolPnt]->getData());

    prepareFluxComputation();

    // calculate the discontinuous flux projected on x, y, z-directions
    for (CFuint iDim = 0; iDim < m_dim; ++iDim)
    {
      // add diffusive part 
      computeFlux(m_avgSol,m_tempGrad,m_cellFluxProjVects[iDim][iSolPnt],0,m_contFlxWoLLAV[iSolPnt][iDim]);
      
      // add convective part
      m_contFlxWoLLAV[iSolPnt][iDim] -= m_updateVarSet->getFlux()(m_pData,m_cellFluxProjVects[iDim][iSolPnt]);
      
      m_contFlx[iSolPnt][iDim] = m_contFlxWoLLAV[iSolPnt][iDim];
      
      // add artificial part
      for (CFuint iDim2 = 0; iDim2 < m_dim; ++iDim2)
      {
        for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
        {
          m_contFlx[iSolPnt][iDim][iVar] += m_solEpsilons[side][iSolPnt]*((*m_avgGradAV[iVar])[iDim2])*m_cellFluxProjVects[iDim][iSolPnt][iDim2];
          
//       if (m_cells[side]->getID() == 0) CFLog(INFO,"first iSol: " << iSolPnt << ", iDir: " << iDim << ", iEq: " << iVar << 
//               ", gradAV: " << ((*m_avgGradAV[iVar])[iDim2]) << ", eps: " << m_solEpsilons[side][iSolPnt] << ", n: " << m_cellFluxProjVects[iDim][iSolPnt][iDim2] << "\n");  
        }
      }
    }

//    for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFlxDep; ++iFlxPnt)
//    {
//      const CFuint flxIdx = (*m_solFlxDep)[iSolPnt][iFlxPnt];
//      const CFuint dim = (*m_flxPntFlxDim)[flxIdx];
//
//      m_extrapolatedFluxes[flxIdx] += (*m_solPolyValsAtFlxPnts)[flxIdx][iSolPnt]*(m_contFlx[iSolPnt][dim]);
//    }
  }
  
  // add the contribution of the faces
  const CFuint nbrFaces = m_cells[side]->nbNeighborGeos();

  for (CFuint iFace = 0; iFace < nbrFaces; ++iFace)
  {
    if (!((*m_isFaceOnBoundaryCell)[iFace]) || m_LLAVBCZero)
    {
      for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
      {
        const CFuint currFlxIdx = (*m_faceFlxPntConn)[iFace][iFlxPnt];

        for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolDep; ++iSolPnt)
        {
          const CFuint solIdx = (*m_flxSolDep)[currFlxIdx][iSolPnt];
          const CFuint dim = (*m_flxPntFlxDim)[currFlxIdx];

           m_extrapolatedFluxes[currFlxIdx] += (*m_solPolyValsAtFlxPnts)[currFlxIdx][solIdx]*(m_contFlx[solIdx][dim]);
        }
      }
    }
    else
    {
      for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
      {
        const CFuint currFlxIdx = (*m_faceFlxPntConn)[iFace][iFlxPnt];

        for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolDep; ++iSolPnt)
        {
          const CFuint solIdx = (*m_flxSolDep)[currFlxIdx][iSolPnt];
          const CFuint dim = (*m_flxPntFlxDim)[currFlxIdx];

           m_extrapolatedFluxes[currFlxIdx] += (*m_solPolyValsAtFlxPnts)[currFlxIdx][solIdx]*(m_contFlxWoLLAV[solIdx][dim]);
        }
      } 
    }
  }

  // Loop over solution pnts to calculate the divergence of the discontinuous flux
  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
  {
    // reset the divergence of FC
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      m_unpertCellDiffRes[side][m_nbrEqs*iSolPnt+iEq] = 0.0;
    }

    // Loop over solution pnt to count factor of all sol pnt polys
    for (CFuint jSolPnt = 0; jSolPnt < m_nbrSolSolDep; ++jSolPnt)
    {
      const CFuint jSolIdx = (*m_solSolDep)[iSolPnt][jSolPnt];

      // Loop over deriv directions and sum them to compute divergence
      for (CFuint iDir = 0; iDir < m_dim; ++iDir)
      {
        const CFreal polyCoef = (*m_solPolyDerivAtSolPnts)[iSolPnt][iDir][jSolIdx]; 

        // Loop over conservative fluxes 
        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
        {
          // Store divFD in the vector that will be divFC
          m_unpertCellDiffRes[side][m_nbrEqs*iSolPnt+iEq] += polyCoef*(m_contFlx[jSolIdx][iDir][iEq]);
	}
      }
    }

    for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFlxDep; ++iFlxPnt)
    {
      const CFuint flxIdx = (*m_solFlxDep)[iSolPnt][iFlxPnt];

      // get the divergence of the correction function
      const CFreal divh = m_corrFctDiv[iSolPnt][flxIdx];
  
      // Fill in the corrections
      for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
      {
        m_unpertCellDiffRes[side][m_nbrEqs*iSolPnt+iVar] += -m_extrapolatedFluxes[flxIdx][iVar] * divh; 
      }
    }
  }
  
  // Reset wall distance in case this is the last face: afterwards diff BC is done which needs the wall dist to be zero!!
  if (m_dim == 2)
  {
    m_navierStokesVarSetTurb->setWallDistance(0.0);
  }
  else
  {
    m_navierStokesVarSetTurb3D->setWallDistance(0.0); 
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffLLAVFluxReconstructionTurb::setup()
{
  CFAUTOTRACE;
  
  // setup parent class
  ConvDiffLLAVFluxReconstructionNS::setup();
  
  const bool RhoivtTv = getMethodData().getUpdateVarStr() == "RhoivtTv";
  
  if(!RhoivtTv)
  {
    // get Euler varset
    m_eulerVarSet = getMethodData().getUpdateVar().d_castTo<EulerVarSet>();
    m_eulerVarSet->getModel()->resizePhysicalData(m_pData);
    m_eulerVarSet->getModel()->resizePhysicalData(m_pData2);
  } 
  else
  {
    m_eulerVarSet2 = getMethodData().getUpdateVar().d_castTo< MultiScalarVarSet< Euler2DVarSet > >();  
    
    m_msEulerTerm = PhysicalModelStack::getActive()-> getImplementor()->getConvectiveTerm().d_castTo< MultiScalarTerm< EulerTerm > >();
    if (m_msEulerTerm.isNull())
    {
      throw Common::ShouldNotBeHereException (FromHere(),"Update variable set is not MultiScalar EulerTerm in BCNoSlipWallrvt!");
    }
  
    m_nbrSpecies = m_msEulerTerm->getNbScalarVars(0);
    
    m_eulerVarSet2->getModel()->resizePhysicalData(m_pData);
    m_eulerVarSet2->getModel()->resizePhysicalData(m_pData2);
  } 
  
  // get the local FR data
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
  cf_assert(frLocalData.size() > 0);
  // for now, there should be only one type of element
  cf_assert(frLocalData.size() == 1);
  
  // get closest sol indices
  m_closestSolToFlxIdx = frLocalData[0]->getClosestSolToFlxIdx();
  
  if (m_dim == 2)
  {
    m_navierStokesVarSetTurb = m_diffusiveVarSet.d_castTo< NavierStokes2DKLogOmega >();
  }
  else
  {  
    m_navierStokesVarSetTurb3D = m_diffusiveVarSet.d_castTo< NavierStokes3DKLogOmega >();
  }
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD
