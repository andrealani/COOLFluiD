// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/MethodCommandProvider.hh"

#include "Framework/CFSide.hh"
#include "Framework/MeshData.hh"
#include "Framework/BaseTerm.hh"

#include "MathTools/MathFunctions.hh"

#include "FluxReconstructionNavierStokes/ConvDiffLLAVFluxReconstructionNS.hh"
#include "FluxReconstructionNavierStokes/FluxReconstructionNavierStokes.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"
#include "NavierStokes/NavierStokesVarSet.hh"

#include "NavierStokes/EulerTerm.hh"
#include "NavierStokes/EulerVarSet.hh"

#include "NavierStokes/Euler2DVarSet.hh"


//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< ConvDiffLLAVFluxReconstructionNS,
		       FluxReconstructionSolverData,
		       FluxReconstructionNavierStokesModule >
convDiffLLAVRHSNSFluxReconstructionProvider("ConvDiffLLAVRHSNS");
  
//////////////////////////////////////////////////////////////////////////////
  
ConvDiffLLAVFluxReconstructionNS::ConvDiffLLAVFluxReconstructionNS(const std::string& name) :
  ConvDiffLLAVFluxReconstruction(name),
  m_tempGradTerm(),
  m_tempGradTermL(),
  m_tempGradTermR(),
  m_diffusiveVarSetNS(CFNULL),
  m_tempStatesL(),
  m_tempStatesR(),
  m_tempStatesL2(),
  m_tempStatesR2(),
  m_tempStatesCell(),
  m_dampCoeffDiff(),
  m_eulerVarSet(CFNULL),
  m_msEulerTerm(CFNULL),
  m_nbrSpecies(),
  m_pData(),
  m_pData2(),
  m_eulerVarSet2(CFNULL),
  m_tempSolVarState(),
  m_tempSolVarState2()
{
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffLLAVFluxReconstructionNS::configure ( Config::ConfigArgs& args )
{
  ConvDiffLLAVFluxReconstruction::configure(args);
} 
  
//////////////////////////////////////////////////////////////////////////////

void ConvDiffLLAVFluxReconstructionNS::computeWaveSpeedUpdates(vector< CFreal >& waveSpeedUpd)
{
  // compute the wave speed updates for the neighbouring cells
  cf_assert(waveSpeedUpd.size() == 2);
  
  // here convective and artificial parts are added!
  ConvDiffLLAVFluxReconstruction::computeWaveSpeedUpdates(waveSpeedUpd);
          
  // now add diffusive part
  CFreal visc = 1.0;

  const CFreal dynVisc = m_diffusiveVarSetNS->getCurrDynViscosity();
  
  const CFreal factorPr = min(m_diffusiveVarSetNS->getModel().getPrandtl(),1.0);
  cf_assert(factorPr>0.0);
  
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    for (CFuint iFlx = 0; iFlx < m_cellFlx[iSide].size(); ++iFlx)
    {
      const CFreal jacobXJacobXIntCoef = m_faceJacobVecAbsSizeFlxPnts[iFlx]*
                                 m_faceJacobVecAbsSizeFlxPnts[iFlx]*
                                   (*m_faceIntegrationCoefs)[iFlx]*
                                   m_cflConvDiffRatio;
      const CFreal rho = m_diffusiveVarSetNS->getDensity(*(m_cellStatesFlxPnt[iSide][iFlx]));
      visc = dynVisc/rho/factorPr;
      
      // transform update states to physical data to calculate eigenvalues
      waveSpeedUpd[iSide] += visc*jacobXJacobXIntCoef/m_cellVolume[iSide];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffLLAVFluxReconstructionNS::computeInterfaceFlxCorrection()
{
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_tempStatesL[iFlx] = m_cellStatesFlxPnt[LEFT][iFlx]->getData();
    m_tempStatesR[iFlx] = m_cellStatesFlxPnt[RIGHT][iFlx]->getData();
  }
  
  m_diffusiveVarSetNS->setGradientVars(m_tempStatesL,m_tempGradTermL,m_nbrFaceFlxPnts);
  m_diffusiveVarSetNS->setGradientVars(m_tempStatesR,m_tempGradTermR,m_nbrFaceFlxPnts);
    
  // Loop over the flux points to calculate FI
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  { 
    // compute the average sol and grad to use the BR2 scheme
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      *(m_avgGrad[iVar]) = (*(m_cellGradFlxPnt[LEFT][iFlxPnt][iVar]) + *(m_cellGradFlxPnt[RIGHT][iFlxPnt][iVar]))/2.0;
      *(m_avgGradAV[iVar]) = (*(m_cellGradFlxPntAV[LEFT][iFlxPnt][iVar]) + *(m_cellGradFlxPntAV[RIGHT][iFlxPnt][iVar]))/2.0;
       
      m_avgSol[iVar] = ((*(m_cellStatesFlxPnt[LEFT][iFlxPnt]))[iVar] + (*(m_cellStatesFlxPnt[RIGHT][iFlxPnt]))[iVar])/2.0; 
    }

    // damping factor
    const CFreal dampFactor = m_dampCoeffDiff*m_faceInvCharLengths[iFlxPnt];

    // compute averaged (damped) gradients
    for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
    {
      // compute damping term
      const RealVector dGradVarXNormal = (m_tempGradTermL(iGrad,iFlxPnt) - m_tempGradTermR(iGrad,iFlxPnt))*m_unitNormalFlxPnts[iFlxPnt];
      *m_avgGrad[iGrad] -= dampFactor*dGradVarXNormal;
      
      // compute damping term for LLAV
      m_tempSolVarState = static_cast<RealVector&>(*m_updateToSolutionVecTrans->transform(m_cellStatesFlxPnt[LEFT][iFlxPnt]));
      m_tempSolVarState2 = static_cast<RealVector&>(*m_updateToSolutionVecTrans->transform(m_cellStatesFlxPnt[RIGHT][iFlxPnt]));
    
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
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffLLAVFluxReconstructionNS::computeCellGradTerm(RealMatrix& gradTerm)
{   
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    m_tempStatesCell[iSol] = (*m_cellStates)[iSol]->getData();
  }
  
  m_diffusiveVarSetNS->setGradientVars(m_tempStatesCell,gradTerm,m_nbrSolPnts);
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffLLAVFluxReconstructionNS::computeFaceGradTerms(RealMatrix& gradTermL, RealMatrix& gradTermR)
{
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_tempStatesL2[iFlx] = m_cellStatesFlxPnt[LEFT][iFlx]->getData();
    m_tempStatesR2[iFlx] = m_cellStatesFlxPnt[RIGHT][iFlx]->getData();
  }
  
  m_diffusiveVarSetNS->setGradientVars(m_tempStatesL2,gradTermL,m_nbrFaceFlxPnts);
  m_diffusiveVarSetNS->setGradientVars(m_tempStatesR2,gradTermR,m_nbrFaceFlxPnts);
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffLLAVFluxReconstructionNS::prepareFluxComputation()
{
  const bool isPerturb = this->getMethodData().isPerturb();
  const CFuint iPerturbVar = this->getMethodData().iPerturbVar();

  m_diffusiveVarSetNS->setComposition(m_avgSol, isPerturb, iPerturbVar);
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffLLAVFluxReconstructionNS::computeGradients()
{
  // get the sol pnt normals
  DataHandle< CFreal > solPntNormals = socket_solPntNormals.getDataHandle();

  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    m_tempStatesCell[iSol] = (*m_cellStates)[iSol]->getData();
    
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      //set the grad updates to 0 
      m_gradUpdates[0][iSol][iEq] = 0.0;
    }
  }

  m_diffusiveVarSetNS->setGradientVars(m_tempStatesCell,m_tempGradTerm,m_nbrSolPnts);
  
  // Loop over solution pnts to calculate the grad updates
  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
  {
    const CFuint solID = (*m_cellStates)[iSolPnt]->getLocalID();

    // Loop over  variables
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      // Loop over gradient directions
      for (CFuint iDir = 0; iDir < m_dim; ++iDir)
      {
        for (CFuint jDir = 0; jDir < m_dim; ++jDir)
        {
	  // project the state on a normal and reuse a RealVector variable of the class to store
	  m_projectedCorrL[jDir] = m_tempGradTerm(iEq,iSolPnt) * solPntNormals[solID*(m_dim+m_ndimplus)*m_dim+(iDir+m_ndimplus)*m_dim+jDir];
        }

        // Loop over solution pnts to count factor of all sol pnt polys
        for (CFuint jSolPnt = 0; jSolPnt < m_nbrSolSolDep; ++jSolPnt)
        { 
          const CFuint jSolIdx = (*m_solSolDep)[iSolPnt][jSolPnt];

          // compute the grad updates
          m_gradUpdates[0][jSolIdx][iEq] += (*m_solPolyDerivAtSolPnts)[jSolIdx][iDir][iSolPnt]*m_projectedCorrL;
	}
      }
    }
  }

  // get the gradients
  DataHandle< vector< RealVector > > gradients = socket_gradients.getDataHandle();
  
  // get the volumes
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();

  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    // get state ID
    const CFuint solID = (*m_cellStates)[iSol]->getLocalID();

    // inverse Jacobian determinant
    const CFreal invJacobDet = 1.0/volumes[solID];

    // update gradients
    for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
    {
      gradients[solID][iGrad] += m_gradUpdates[0][iSol][iGrad];
      gradients[solID][iGrad] *= invJacobDet;
    }
  }

  // compute the gradients for the artificial viscosity
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      {
        //set the grad updates to 0 
        m_gradUpdates[0][iSol][iEq] = 0.0;
      }
    }

    // Loop over solution pnts to calculate the grad updates
    for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
    {
      const CFuint solID = (*m_cellStates)[iSolPnt]->getLocalID();
      
      m_tempSolVarState = static_cast<RealVector&>(*m_updateToSolutionVecTrans->transform((*m_cellStates)[iSolPnt]));
              
      // Loop over  variables
      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      {
        // Loop over gradient directions
        for (CFuint iDir = 0; iDir < m_dim; ++iDir)
        {
          for (CFuint jDir = 0; jDir < m_dim; ++jDir)
          {
	    // project the state on a normal and reuse a RealVector variable of the class to store
	    m_projectedCorrL[jDir] = m_tempSolVarState[iEq] * solPntNormals[solID*(m_dim+m_ndimplus)*m_dim+(iDir+m_ndimplus)*m_dim+jDir];
          }

          // Loop over solution pnts to count factor of all sol pnt polys
          for (CFuint jSolPnt = 0; jSolPnt < m_nbrSolSolDep; ++jSolPnt)
          {
            const CFuint jSolIdx = (*m_solSolDep)[iSolPnt][jSolPnt];
	  
            // compute the grad updates
            m_gradUpdates[0][jSolIdx][iEq] += (*m_solPolyDerivAtSolPnts)[jSolIdx][iDir][iSolPnt]*m_projectedCorrL;
	  }
        }
      }
    }
    
    // get the gradientsAV
    DataHandle< vector< RealVector > > gradientsAV = socket_gradientsAV.getDataHandle();

    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      // get state ID
      const CFuint solID = (*m_cellStates)[iSol]->getLocalID();

      // inverse Jacobian determinant
      const CFreal invJacobDet = 1.0/volumes[solID];

      // update gradientsAV
      for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
      {
        gradientsAV[solID][iGrad] += m_gradUpdates[0][iSol][iGrad];
        gradientsAV[solID][iGrad] *= invJacobDet;
      }
    }
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffLLAVFluxReconstructionNS::computeGradientFaceCorrections()
{
  // get the face flux point normals
  DataHandle< CFreal > flxPntNormals = socket_flxPntNormals.getDataHandle();
  
  // Loop over flux points to set the normal vectors
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  { 
    for (CFuint iDim = 0; iDim < m_dim; ++iDim)
    {
      m_faceJacobVecs[iFlxPnt][iDim] = flxPntNormals[m_face->getID()*m_nbrFaceFlxPnts*m_dim+iFlxPnt*m_dim+iDim];
    }
  }
  
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_tempStatesL2[iFlx] = m_cellStatesFlxPnt[LEFT][iFlx]->getData();
    m_tempStatesR2[iFlx] = m_cellStatesFlxPnt[RIGHT][iFlx]->getData();
  }
  
  m_diffusiveVarSetNS->setGradientVars(m_tempStatesL2,m_tempGradTermL,m_nbrFaceFlxPnts);
  m_diffusiveVarSetNS->setGradientVars(m_tempStatesR2,m_tempGradTermR,m_nbrFaceFlxPnts);
  
  // Loop over solution pnts to reset the grad updates
  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
  {
    // Loop over  variables
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      //set the grad updates to 0 
      m_gradUpdates[LEFT][iSolPnt][iEq] = 0.0;
      m_gradUpdates[RIGHT][iSolPnt][iEq] = 0.0;
    }
  }
  
  // compute the face corrections to the gradients
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    const CFuint flxIdxL = (*m_faceFlxPntConnPerOrient)[m_orient][LEFT][iFlx];
    const CFuint flxIdxR = (*m_faceFlxPntConnPerOrient)[m_orient][RIGHT][iFlx];

    // Loop over  variables
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      const CFreal avgSol = (m_tempGradTermL(iEq,iFlx)+m_tempGradTermR(iEq,iFlx))/2.0;
      m_projectedCorrL = (avgSol-m_tempGradTermL(iEq,iFlx))*m_faceJacobVecSizeFlxPnts[iFlx][LEFT]*m_unitNormalFlxPnts[iFlx];
      m_projectedCorrR = (avgSol-m_tempGradTermR(iEq,iFlx))*m_faceJacobVecSizeFlxPnts[iFlx][RIGHT]*m_unitNormalFlxPnts[iFlx];

      // Loop over solution pnts to calculate the grad updates
      for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolDep; ++iSolPnt)
      {
        const CFuint iSolIdxL = (*m_flxSolDep)[flxIdxL][iSolPnt];
        const CFuint iSolIdxR = (*m_flxSolDep)[flxIdxR][iSolPnt];

	/// @todo Check if this is also OK for triangles!!
	m_gradUpdates[LEFT][iSolIdxL][iEq] += m_projectedCorrL*m_corrFctDiv[iSolIdxL][flxIdxL];
	m_gradUpdates[RIGHT][iSolIdxR][iEq] += m_projectedCorrR*m_corrFctDiv[iSolIdxR][flxIdxR];
      }
    }
  }
  
  // get the gradients
  DataHandle< vector< RealVector > > gradients = socket_gradients.getDataHandle();

  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      // get state ID
      const CFuint solID = (*m_states[iSide])[iSol]->getLocalID();

      // update gradients
      for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
      {
        gradients[solID][iGrad] += m_gradUpdates[iSide][iSol][iGrad];
      }
    }
  }
  
  // compute the gradients for the artificial viscosity
    // Loop over solution pnts to reset the grad updates
    for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
    {
      // Loop over  variables
      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      {
        //set the grad updates to 0 
        m_gradUpdates[LEFT][iSolPnt][iEq] = 0.0;
        m_gradUpdates[RIGHT][iSolPnt][iEq] = 0.0;
      }
    }

    // Loop over  variables
    for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
    {
      const CFuint flxIdxL = (*m_faceFlxPntConnPerOrient)[m_orient][LEFT][iFlx];
      const CFuint flxIdxR = (*m_faceFlxPntConnPerOrient)[m_orient][RIGHT][iFlx];
      
      
      m_tempSolVarState = static_cast<RealVector&>(*m_updateToSolutionVecTrans->transform(m_cellStatesFlxPnt[LEFT][iFlx]));
      m_tempSolVarState2 = static_cast<RealVector&>(*m_updateToSolutionVecTrans->transform(m_cellStatesFlxPnt[RIGHT][iFlx]));

      // compute the face corrections to the gradients
      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      {
        const CFreal avgSol = (m_tempSolVarState[iEq]+m_tempSolVarState2[iEq])/2.0;
	m_projectedCorrL = (avgSol-m_tempSolVarState[iEq])*m_faceJacobVecSizeFlxPnts[iFlx][LEFT]*m_unitNormalFlxPnts[iFlx];
	m_projectedCorrR = (avgSol-m_tempSolVarState2[iEq])*m_faceJacobVecSizeFlxPnts[iFlx][RIGHT]*m_unitNormalFlxPnts[iFlx];

        // Loop over solution pnts to calculate the grad updates
        for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolDep; ++iSolPnt)
        {
          const CFuint iSolIdxL = (*m_flxSolDep)[flxIdxL][iSolPnt];
          const CFuint iSolIdxR = (*m_flxSolDep)[flxIdxR][iSolPnt];

	  /// @todo Check if this is also OK for triangles!!
	  m_gradUpdates[LEFT][iSolIdxL][iEq] += m_projectedCorrL*m_corrFctDiv[iSolIdxL][flxIdxL];
	  m_gradUpdates[RIGHT][iSolIdxR][iEq] += m_projectedCorrR*m_corrFctDiv[iSolIdxR][flxIdxR];
        }
      }
    }
  
    // get the gradientsAV
    DataHandle< vector< RealVector > > gradientsAV = socket_gradientsAV.getDataHandle();

    for (CFuint iSide = 0; iSide < 2; ++iSide)
    {
      for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
      {
        // get state ID
        const CFuint solID = (*m_states[iSide])[iSol]->getLocalID();

        // update gradientsAV
        for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
        {
          gradientsAV[solID][iGrad] += m_gradUpdates[iSide][iSol][iGrad];
        }
      }
    }
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffLLAVFluxReconstructionNS::computeSmoothness()
{ 
  CFreal sNum = 0.0;
  
  CFreal sDenom = 0.0;
  
  // get datahandle
  DataHandle< CFreal > monPhysVar = socket_monPhysVar.getDataHandle();
  
  const bool RhoivtTv = getMethodData().getUpdateVarStr() == "RhoivtTv";
  
  if (RhoivtTv && m_monitoredPhysVar < m_pData.size())
  {
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      CFreal stateP = 0.0;
      CFreal diffStatesPPMinOne = 0.0;

      m_eulerVarSet2->computePhysicalData(*((*m_cellStates)[iSol]),m_pData);
      m_eulerVarSet2->computePhysicalData(m_statesPMinOne[iSol],m_pData2);

      stateP = m_pData[m_monitoredPhysVar];
      diffStatesPPMinOne = stateP - m_pData2[m_monitoredPhysVar];
    
      monPhysVar[(((*m_cellStates)[iSol]))->getLocalID()] = stateP;

      sNum += diffStatesPPMinOne*diffStatesPPMinOne;
      sDenom += stateP*stateP;
    }
  }
  else if (!RhoivtTv && m_monitoredPhysVar < m_pData.size())
  {
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      CFreal stateP = 0.0;
      CFreal diffStatesPPMinOne = 0.0;

      m_eulerVarSet->computePhysicalData(*((*m_cellStates)[iSol]),m_pData);
      m_eulerVarSet->computePhysicalData(m_statesPMinOne[iSol],m_pData2);

      stateP = m_pData[m_monitoredPhysVar];
      diffStatesPPMinOne = stateP - m_pData2[m_monitoredPhysVar];
    
      monPhysVar[(((*m_cellStates)[iSol]))->getLocalID()] = stateP;

      sNum += diffStatesPPMinOne*diffStatesPPMinOne;
      sDenom += stateP*stateP;
    }
  }
  else
  {
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      CFreal stateP = 0.0;
      CFreal diffStatesPPMinOne = 0.0;

      stateP = (*((*m_cellStates)[iSol]))[m_monitoredVar];
      diffStatesPPMinOne = stateP - m_statesPMinOne[iSol][m_monitoredVar];

      sNum += diffStatesPPMinOne*diffStatesPPMinOne;
      sDenom += stateP*stateP;
    }
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
  DataHandle< CFreal > smoothness = socket_smoothness.getDataHandle();
  
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    smoothness[(((*m_cellStates)[iSol]))->getLocalID()] = m_s;
  }
  
  if (m_s > m_Smax)
  {
    m_Smax = m_s;
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffLLAVFluxReconstructionNS::computeSmoothness(const CFuint side)
{ 
  CFreal sNum = 0.0;
  
  CFreal sDenom = 0.0;
  
  // get datahandle
  DataHandle< CFreal > monPhysVar = socket_monPhysVar.getDataHandle();
  
  const bool RhoivtTv = getMethodData().getUpdateVarStr() == "RhoivtTv";
  
  if (RhoivtTv && m_monitoredPhysVar < m_pData.size())
  {
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      CFreal stateP = 0.0;
      CFreal diffStatesPPMinOne = 0.0;

      m_eulerVarSet2->computePhysicalData(*((*m_states[side])[iSol]),m_pData);
      m_eulerVarSet2->computePhysicalData(m_statesPMinOne[iSol],m_pData2);

      stateP = m_pData[m_monitoredPhysVar];
      diffStatesPPMinOne = stateP - m_pData2[m_monitoredPhysVar];
    
      monPhysVar[(((*m_cellStates)[iSol]))->getLocalID()] = stateP;

      sNum += diffStatesPPMinOne*diffStatesPPMinOne;
      sDenom += stateP*stateP;
    }
  }
  else if (!RhoivtTv && m_monitoredPhysVar < m_pData.size())
  {
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      CFreal stateP = 0.0;
      CFreal diffStatesPPMinOne = 0.0;

      m_eulerVarSet->computePhysicalData(*((*m_states[side])[iSol]),m_pData);
      m_eulerVarSet->computePhysicalData(m_statesPMinOne[iSol],m_pData2);

      stateP = m_pData[m_monitoredPhysVar];
      diffStatesPPMinOne = stateP - m_pData2[m_monitoredPhysVar];
    
      monPhysVar[(((*m_cellStates)[iSol]))->getLocalID()] = stateP;

      sNum += diffStatesPPMinOne*diffStatesPPMinOne;
      sDenom += stateP*stateP;
    }
  }
  else
  {
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      CFreal stateP = 0.0;
      CFreal diffStatesPPMinOne = 0.0;

      stateP = (*((*m_cellStates)[iSol]))[m_monitoredVar];
      diffStatesPPMinOne = stateP - m_statesPMinOne[iSol][m_monitoredVar];

      sNum += diffStatesPPMinOne*diffStatesPPMinOne;
      sDenom += stateP*stateP;
    }
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

void ConvDiffLLAVFluxReconstructionNS::computeEpsilon0()
{ 
  // compute a cell average characteristic flow speed. Note that a straight average is used, not a weighted one, maybe change this
  CFreal wavespeed = 0.0;
  
  const bool RhoivtTv = getMethodData().getUpdateVarStr() == "RhoivtTv";
  
  if (RhoivtTv)
  {
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      m_eulerVarSet2->computePhysicalData(*((*m_cellStates)[iSol]),m_pData);
      
      wavespeed += m_pData[EulerTerm::V] + m_pData[EulerTerm::A];
      cf_assert(m_pData[EulerTerm::V] > 0.0 && m_pData[EulerTerm::A] > 0.0);
    }
  }
  else
  {
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      m_eulerVarSet->computePhysicalData(*((*m_cellStates)[iSol]),m_pData);

      wavespeed += m_pData[EulerTerm::V] + m_pData[EulerTerm::A];
      cf_assert(m_pData[EulerTerm::V] > 0.0 && m_pData[EulerTerm::A] > 0.0);
    }
  }
  
  cf_assert(wavespeed > 0.0);
  
  wavespeed /= m_nbrSolPnts;
    
  const CFreal peclet = computePeclet();
  
  const CFreal oneOverDim = 1./m_dim;
  
  // get the cell volumes
  DataHandle< CFreal > cellVolumes = socket_cellVolumes.getDataHandle();
  
  const CFreal h = pow(cellVolumes[m_cell->getID()],oneOverDim);
  
  m_epsilon0 = max(h*wavespeed*(2.0/peclet - m_subcellRes/peclet),0.0);
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffLLAVFluxReconstructionNS::computeEpsilon0(const CFuint side)
{ 
  // compute a cell average characteristic flow speed. Note that a straight average is used, not a weighted one, maybe change this
  CFreal wavespeed = 0.0;
  
  const bool RhoivtTv = getMethodData().getUpdateVarStr() == "RhoivtTv";
  
  if (RhoivtTv)
  {
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      m_eulerVarSet2->computePhysicalData(*((*m_states[side])[iSol]),m_pData);
      
      wavespeed += m_pData[EulerTerm::V] + m_pData[EulerTerm::A];
      cf_assert(m_pData[EulerTerm::V] > 0.0 && m_pData[EulerTerm::A] > 0.0);
    }
  }
  else
  {
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      m_eulerVarSet->computePhysicalData(*((*m_states[side])[iSol]),m_pData);

      wavespeed += m_pData[EulerTerm::V] + m_pData[EulerTerm::A];
      cf_assert(m_pData[EulerTerm::V] > 0.0 && m_pData[EulerTerm::A] > 0.0);
    }
  }
  
  cf_assert(wavespeed > 0.0);
  
  wavespeed /= m_nbrSolPnts;
    
  const CFreal peclet = computePeclet();
  
  const CFreal oneOverDim = 1./m_dim;
  
  // get the cell volumes
  DataHandle< CFreal > cellVolumes = socket_cellVolumes.getDataHandle();
  
  const CFreal h = pow(cellVolumes[m_cells[side]->getID()],oneOverDim);
  
  m_epsilon0 = max(h*wavespeed*(2.0/peclet - m_subcellRes/peclet),0.0);
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffLLAVFluxReconstructionNS::setup()
{
  CFAUTOTRACE;
  
  // setup parent class
  ConvDiffLLAVFluxReconstruction::setup();
  
  // get damping coeff
  m_dampCoeffDiff = getMethodData().getDiffDampCoefficient();
  
  // get the diffusive varset
  m_diffusiveVarSetNS = (getMethodData().getDiffusiveVar()).d_castTo< NavierStokesVarSet >();
  
  m_tempSolVarState.resize(m_nbrEqs);
  m_tempSolVarState2.resize(m_nbrEqs);
  
  m_tempGradTerm.resize(m_nbrEqs,m_nbrSolPnts);
  
  m_tempGradTermL.resize(m_nbrEqs,m_nbrFaceFlxPnts);
  m_tempGradTermR.resize(m_nbrEqs,m_nbrFaceFlxPnts);
  
  m_tempStatesL.resize(m_nbrFaceFlxPnts);
  m_tempStatesR.resize(m_nbrFaceFlxPnts);
  m_tempStatesL2.resize(m_nbrFaceFlxPnts);
  m_tempStatesR2.resize(m_nbrFaceFlxPnts);
  
  m_tempStatesCell.resize(m_nbrSolPnts);
  
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
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

