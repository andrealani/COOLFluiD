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

#include "FluxReconstructionNavierStokes/ConvDiffJacobFluxReconstructionNS.hh"
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

MethodCommandProvider< ConvDiffJacobFluxReconstructionNS,
		       FluxReconstructionSolverData,
		       FluxReconstructionNavierStokesModule >
convDiffRHSJacobNSFluxReconstructionProvider("ConvDiffRHSJacobNS");
  
//////////////////////////////////////////////////////////////////////////////
  
ConvDiffJacobFluxReconstructionNS::ConvDiffJacobFluxReconstructionNS(const std::string& name) :
  ConvDiffJacobFluxReconstruction(name),
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
  m_tempGradTermJacob(),
  m_tempStatesJacob(),
  m_tempGradTermJacob2(),
  m_tempStatesJacob2(),
  m_unpertGradVars(),
  m_pertGradVars()
{
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffJacobFluxReconstructionNS::configure ( Config::ConfigArgs& args )
{
  ConvDiffJacobFluxReconstruction::configure(args);
} 
  
//////////////////////////////////////////////////////////////////////////////

void ConvDiffJacobFluxReconstructionNS::computeWaveSpeedUpdates(vector< CFreal >& waveSpeedUpd)
{
  // compute the wave speed updates for the neighbouring cells
  cf_assert(waveSpeedUpd.size() == 2);
  
  // here convective and artificial parts are added!
  ConvDiffJacobFluxReconstruction::computeWaveSpeedUpdates(waveSpeedUpd);
          
  // now add diffusive part
  CFreal visc = 1.0;

  const CFreal dynVisc = m_diffusiveVarSetNS->getCurrDynViscosity();
  
  const CFreal factorPr = min(m_diffusiveVarSetNS->getModel().getPrandtl(),1.0);
  cf_assert(factorPr>0.0);
    
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    //for (CFuint iFlx = 0; iFlx < m_cellFlx[iSide].size(); ++iFlx)
    for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
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

void ConvDiffJacobFluxReconstructionNS::computeInterfaceFlxCorrection()
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
    }
    
    prepareFluxComputation();
     
    // compute the diff flux
    computeFlux(m_avgSol,m_avgGrad,m_unitNormalFlxPnts[iFlxPnt],0,m_flxPntRiemannFluxDiff[iFlxPnt]);
    
    m_flxPntRiemannFlux[iFlxPnt] = m_flxPntRiemannFluxDiff[iFlxPnt];
    
    // compute the convective riemann flux
    m_flxPntRiemannFlux[iFlxPnt] -= m_riemannFluxComputer->computeFlux(*(m_cellStatesFlxPnt[LEFT][iFlxPnt]),
									    *(m_cellStatesFlxPnt[RIGHT][iFlxPnt]),
									    m_unitNormalFlxPnts[iFlxPnt]);
     
    // compute FI in the mapped coord frame
    m_cellFlx[LEFT][iFlxPnt] = (m_flxPntRiemannFlux[iFlxPnt])*m_faceJacobVecSizeFlxPnts[iFlxPnt][LEFT];
    m_cellFlx[RIGHT][iFlxPnt] = (m_flxPntRiemannFlux[iFlxPnt])*m_faceJacobVecSizeFlxPnts[iFlxPnt][RIGHT];
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffJacobFluxReconstructionNS::computeRiemannFluxJacobianNum(const CFreal resFactor)
{
  CFLog(VERBOSE, "NS computeRiemannFluxJacobianNum\n");
  
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_tempStatesL[iFlx] = m_cellStatesFlxPnt[LEFT][iFlx]->getData();
    m_tempStatesR[iFlx] = m_cellStatesFlxPnt[RIGHT][iFlx]->getData();
  }
      
  // loop over the face flux points to compute the Riemann flux jacobian
  for (m_pertSide = 0; m_pertSide < 2; ++m_pertSide)
  {
    // variable for the other side
    const CFuint iOtherSide = m_pertSide == LEFT ? RIGHT : LEFT;
    
    for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
    {
      // dereference state
      State& pertState = *(m_cellStatesFlxPnt[m_pertSide][iFlxPnt]);
      
      // damping factor
      const CFreal dampFactor = m_dampCoeffDiff*m_faceInvCharLengths[iFlxPnt];
                
      // loop over the variables in the state
      for (m_pertVar = 0; m_pertVar < m_nbrEqs; ++m_pertVar)
      {
        // perturb physical variable in state
        m_numJacob->perturb(m_pertVar,pertState[m_pertVar]);
        
        if (m_pertSide == LEFT)
        {
          *(m_tempStatesL[iFlxPnt]) = pertState;
        }
        else
        {
          *(m_tempStatesR[iFlxPnt]) = pertState;
        }
  
        m_diffusiveVarSetNS->setGradientVars(m_tempStatesL,m_tempGradTermL,m_nbrFaceFlxPnts);
        m_diffusiveVarSetNS->setGradientVars(m_tempStatesR,m_tempGradTermR,m_nbrFaceFlxPnts);
        
        // compute the average grad to use the BR2 scheme
        for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
        {
          *(m_avgGrad[iVar]) = (*(m_cellGradFlxPnt[LEFT][iFlxPnt][iVar]) + *(m_cellGradFlxPnt[RIGHT][iFlxPnt][iVar]))/2.0;          
        }

        // compute averaged (damped) gradients
        for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
        {
          // compute damping term
          const RealVector dGradVarXNormal = (m_tempGradTermL(iGrad,iFlxPnt) - m_tempGradTermR(iGrad,iFlxPnt))*m_unitNormalFlxPnts[iFlxPnt];
          *m_avgGrad[iGrad] -= dampFactor*dGradVarXNormal;
        }
        
        // compute the average sol
        for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
        {        
          m_avgSol[iVar] = (pertState[iVar] + (*(m_cellStatesFlxPnt[iOtherSide][iFlxPnt]))[iVar])/2.0; 
        }
    
        prepareFluxComputation();
     
        // compute diffusive flux
        computeFlux(m_avgSol,m_avgGrad,m_unitNormalFlxPnts[iFlxPnt],0,m_flxPntRiemannFluxPert[iFlxPnt]);
        
        if (m_pertSide == LEFT)
        {
          // compute the convective riemann flux
          m_flxPntRiemannFluxPert[iFlxPnt] -= m_riemannFluxComputer->computeFlux(pertState,
									              *(m_cellStatesFlxPnt[RIGHT][iFlxPnt]),
									              m_unitNormalFlxPnts[iFlxPnt]);
        }
        else
        {
          // compute the convective riemann flux
          m_flxPntRiemannFluxPert[iFlxPnt] -= m_riemannFluxComputer->computeFlux(*(m_cellStatesFlxPnt[LEFT][iFlxPnt]),
                                                                                      pertState,
									              m_unitNormalFlxPnts[iFlxPnt]);
        }
       
        // compute the flux current jacobian term
        // compute the finite difference derivative of the face term (note an implicit minus sign is added here, by the ordering of arguments)
        m_numJacob->computeDerivative(m_flxPntRiemannFluxPert[iFlxPnt],m_flxPntRiemannFlux[iFlxPnt],m_riemannFluxJacobian[m_pertSide][iFlxPnt][m_pertVar]);

        // multiply residual update derivatives with residual factor so it is taken into the final jacobian
        m_riemannFluxJacobian[m_pertSide][iFlxPnt][m_pertVar] *= resFactor;

        // restore physical variable in state
        m_numJacob->restore(pertState[m_pertVar]);
        
        // restore tempStates
        if (m_pertSide == LEFT)
        {
          m_tempStatesL[iFlxPnt] = m_cellStatesFlxPnt[m_pertSide][iFlxPnt]->getData();
        }
        else
        {
          m_tempStatesR[iFlxPnt] = m_cellStatesFlxPnt[m_pertSide][iFlxPnt]->getData();
        }
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffJacobFluxReconstructionNS::computeRiemannFluxToGradJacobianNum(const CFreal resFactor)
{  
  CFLog(VERBOSE, "NS computeRiemannFluxToGradJacobianNum\n");

  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_tempStatesL[iFlx] = m_cellStatesFlxPnt[LEFT][iFlx]->getData();
    m_tempStatesR[iFlx] = m_cellStatesFlxPnt[RIGHT][iFlx]->getData();
  }
  
  m_diffusiveVarSetNS->setGradientVars(m_tempStatesL,m_tempGradTermL,m_nbrFaceFlxPnts);
  m_diffusiveVarSetNS->setGradientVars(m_tempStatesR,m_tempGradTermR,m_nbrFaceFlxPnts); 

  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {   
    // damping factor
    const CFreal dampFactor = m_dampCoeffDiff*m_faceInvCharLengths[iFlxPnt]; 
    
    // compute the average grad to use the BR2 scheme
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      *(m_avgGrad[iVar]) = (*(m_cellGradFlxPnt[LEFT][iFlxPnt][iVar]) + *(m_cellGradFlxPnt[RIGHT][iFlxPnt][iVar]))/2.0;            
    }

    // compute the average sol
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {        
      m_avgSol[iVar] = ((*(m_cellStatesFlxPnt[LEFT][iFlxPnt]))[iVar] + (*(m_cellStatesFlxPnt[RIGHT][iFlxPnt]))[iVar])/2.0; 
    }
  
    // loop over the variables in the state
    for (m_pertVar = 0; m_pertVar < m_nbrEqs; ++m_pertVar)
    {
      for (CFuint pertDir = 0; pertDir < m_dim; ++pertDir)
      {
        // perturb physical variable in state
        m_numJacob->perturb(m_pertVar,(*(m_avgGrad[m_pertVar]))[pertDir]);

        // compute averaged (damped) gradients
        for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
        {
          // compute damping term
          const RealVector dGradVarXNormal = (m_tempGradTermL(iGrad,iFlxPnt) - m_tempGradTermR(iGrad,iFlxPnt))*m_unitNormalFlxPnts[iFlxPnt];
          *m_tempGrad[iGrad] = *(m_avgGrad[iGrad]) - dampFactor*dGradVarXNormal;
        }

        prepareFluxComputation();

        // compute diffusive flux
        computeFlux(m_avgSol,m_tempGrad,m_unitNormalFlxPnts[iFlxPnt],0,m_flxPntRiemannFluxPert[iFlxPnt]);

        // compute the flux current jacobian term
        // compute the finite difference derivative of the face term (note an implicit minus sign is added here, by the ordering of arguments)
        m_numJacob->computeDerivative(m_flxPntRiemannFluxPert[iFlxPnt],m_flxPntRiemannFluxDiff[iFlxPnt],m_riemannFluxGradJacobian[iFlxPnt][m_pertVar][pertDir]);

        // multiply residual update derivatives with residual factor so it is taken into the final jacobian
        m_riemannFluxGradJacobian[iFlxPnt][m_pertVar][pertDir] *= resFactor;

        // restore physical variable in state
        m_numJacob->restore((*(m_avgGrad[m_pertVar]))[pertDir]);
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffJacobFluxReconstructionNS::computeGradVarsToStateJacobianNum()
{
  CFLog(VERBOSE, "NS computeGradVarsToStateJacobianNum\n");

  for (m_pertSide = 0; m_pertSide < 2; ++m_pertSide)
  { 
    for (m_pertSol = 0; m_pertSol < m_nbrSolPnts; ++m_pertSol)
    {
      m_tempStatesJacob[0] = (*m_states[m_pertSide])[m_pertSol];
    
      m_diffusiveVarSetNS->setGradientVars(m_tempStatesJacob,m_tempGradTermJacob,1);
      
      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      {
        m_unpertGradVars[iEq] = m_tempGradTermJacob(iEq,0);
      }

      // dereference state
      State& pertState = *(*m_states[m_pertSide])[m_pertSol];

      for (m_pertVar = 0; m_pertVar < m_nbrEqs; ++m_pertVar)
      {
        // perturb physical variable in state
        m_numJacob->perturb(m_pertVar,pertState[m_pertVar]);
        
        m_tempStatesJacob2[0] = &pertState;
  
        m_diffusiveVarSetNS->setGradientVars(m_tempStatesJacob2,m_tempGradTermJacob2,1);
        
        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
        {
          m_pertGradVars[iEq] = m_tempGradTermJacob2(iEq,0);
        }
        
        m_numJacob->computeDerivative(m_unpertGradVars,m_pertGradVars,m_gradVarsToStateJacobian[m_pertSide][m_pertSol][m_pertVar]);

        // restore physical variable in state
        m_numJacob->restore(pertState[m_pertVar]);
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffJacobFluxReconstructionNS::prepareFluxComputation()
{
  const bool isPerturb = this->getMethodData().isPerturb();
  const CFuint iPerturbVar = this->getMethodData().iPerturbVar();

  m_diffusiveVarSetNS->setComposition(m_avgSol, isPerturb, iPerturbVar);
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffJacobFluxReconstructionNS::computeGradients()
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
	  m_projectedCorrL[jDir] = m_tempGradTerm(iEq,iSolPnt) * solPntNormals[solID*(m_dim+m_ndimplus)*m_dim+iDir*m_dim+jDir]; //Modified
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
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffJacobFluxReconstructionNS::computeGradientFaceCorrections()
{
  // get the face flux point normals
  DataHandle< CFreal > flxPntNormals = socket_flxPntNormals.getDataHandle();
  
  // Loop over flux points to set the normal vectors
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  { 
    for (CFuint iDim = 0; iDim < m_dim; ++iDim)
    {
      m_faceJacobVecs[iFlxPnt][iDim] = flxPntNormals[m_face->getID()*m_nbFaceFlxPntsMax*m_dim+iFlxPnt*m_dim+iDim];
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
    
    m_nbrSolDep = ((*m_flxSolDep)[flxIdxL]).size();
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
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffJacobFluxReconstructionNS::setup()
{
  CFAUTOTRACE;
  
  // setup parent class
  ConvDiffJacobFluxReconstruction::setup();
  
  // get damping coeff
  m_dampCoeffDiff = getMethodData().getDiffDampCoefficient();
  
  // get the diffusive varset
  m_diffusiveVarSetNS = (getMethodData().getDiffusiveVar()).d_castTo< NavierStokesVarSet >();
  
  m_updateToSolutionVecTrans = getMethodData().getUpdateToSolutionVecTrans();
  
  m_updateToSolutionVecTrans->setup(2);
  
  m_tempGradTerm.resize(m_nbrEqs,m_nbrSolPnts);
  m_tempGradTermJacob.resize(m_nbrEqs,1);
  m_tempGradTermJacob2.resize(m_nbrEqs,1);

  
  m_tempGradTermL.resize(m_nbrEqs,m_nbrFaceFlxPnts);
  m_tempGradTermR.resize(m_nbrEqs,m_nbrFaceFlxPnts);
  
  m_tempStatesL.resize(m_nbrFaceFlxPnts);
  m_tempStatesR.resize(m_nbrFaceFlxPnts);
  m_tempStatesL2.resize(m_nbrFaceFlxPnts);
  m_tempStatesR2.resize(m_nbrFaceFlxPnts);
  m_tempStatesJacob.resize(1);
  m_tempStatesJacob2.resize(1);
  m_unpertGradVars.resize(m_nbrEqs);
  m_pertGradVars.resize(m_nbrEqs);
  
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
