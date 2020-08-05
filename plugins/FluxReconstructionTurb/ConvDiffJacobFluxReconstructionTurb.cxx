// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/MethodCommandProvider.hh"

#include "Framework/CFSide.hh"
#include "Framework/MeshData.hh"
#include "Framework/BaseTerm.hh"

#include "MathTools/MathFunctions.hh"

#include "FluxReconstructionTurb/ConvDiffJacobFluxReconstructionTurb.hh"
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

MethodCommandProvider< ConvDiffJacobFluxReconstructionTurb,
		       FluxReconstructionSolverData,
		       FluxReconstructionNavierStokesModule >
convDiffRHSJacobTurbFluxReconstructionProvider("ConvDiffRHSJacobTurb");
  
//////////////////////////////////////////////////////////////////////////////
  
ConvDiffJacobFluxReconstructionTurb::ConvDiffJacobFluxReconstructionTurb(const std::string& name) :
  ConvDiffJacobFluxReconstructionNS(name),
  socket_wallDistance("wallDistance"),
  m_closestSolToFlxIdx(CFNULL),
  m_pData(),
  m_pData2()
{
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffJacobFluxReconstructionTurb::configure ( Config::ConfigArgs& args )
{
  ConvDiffJacobFluxReconstructionNS::configure(args);
} 

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > >
    ConvDiffJacobFluxReconstructionTurb::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result = ConvDiffJacobFluxReconstructionNS::needsSockets();

  result.push_back(&socket_wallDistance);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffJacobFluxReconstructionTurb::computeInterfaceFlxCorrection()
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
  
  SafePtr< NavierStokes2DKLogOmega > navierStokesVarSet = m_diffusiveVarSet.d_castTo< NavierStokes2DKLogOmega >();
    
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
    navierStokesVarSet->setWallDistance(flxPntWallDist);
    
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

void ConvDiffJacobFluxReconstructionTurb::computeUnpertCellDiffResiduals(const CFuint side)
{
  // get the sol pnt normals
  DataHandle< CFreal > solPntNormals = socket_solPntNormals.getDataHandle();
  
  SafePtr< NavierStokes2DKLogOmega > navierStokesVarSet = m_diffusiveVarSet.d_castTo< NavierStokes2DKLogOmega >();
  
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
  
  // reset the extrapolated fluxes
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrTotalFlxPnts; ++iFlxPnt)
  {
    m_extrapolatedFluxes[iFlxPnt] = 0.0;
  }
  
  // Get the wall distance
  DataHandle< CFreal > wallDist = socket_wallDistance.getDataHandle();

  // Loop over solution points to calculate the discontinuous flux.
  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
  { 
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      *(m_tempGrad[iVar]) = (*(m_cellGrads[side][iSolPnt]))[iVar];
    }
    
    const CFuint stateID = (*(m_states[side]))[iSolPnt]->getLocalID();
    
    // Set the wall distance before computing the turbulent viscosity
    navierStokesVarSet->setWallDistance(wallDist[stateID]);
    
    m_updateVarSet->computePhysicalData(*((*(m_states[side]))[iSolPnt]), m_pData); 

    m_avgSol = *((*(m_states[side]))[iSolPnt]->getData());

    prepareFluxComputation();

    // calculate the discontinuous flux projected on x, y, z-directions
    for (CFuint iDim = 0; iDim < m_dim; ++iDim)
    {
      // add diffusive part 
      computeFlux(m_avgSol,m_tempGrad,m_cellFluxProjVects[iDim][iSolPnt],0,m_contFlx[iSolPnt][iDim]);
      
      // add convective part
      m_contFlx[iSolPnt][iDim] -= m_updateVarSet->getFlux()(m_pData,m_cellFluxProjVects[iDim][iSolPnt]);
    }

    for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFlxDep; ++iFlxPnt)
    {
      const CFuint flxIdx = (*m_solFlxDep)[iSolPnt][iFlxPnt];
      const CFuint dim = (*m_flxPntFlxDim)[flxIdx];

      m_extrapolatedFluxes[flxIdx] += (*m_solPolyValsAtFlxPnts)[flxIdx][iSolPnt]*(m_contFlx[iSolPnt][dim]);
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
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffJacobFluxReconstructionTurb::initJacobianComputation()
{
  CFLog(VERBOSE, "initJacobianComputation\n");
    
  // get the sol pnt normals
  DataHandle< CFreal > solPntNormals = socket_solPntNormals.getDataHandle();
  
  // Get the wall distance
  DataHandle< CFreal > wallDist = socket_wallDistance.getDataHandle();
  
  SafePtr< NavierStokes2DKLogOmega > navierStokesVarSet = m_diffusiveVarSet.d_castTo< NavierStokes2DKLogOmega >();
  
  // set block row and column indices, proj vectors and make a backup of discontinuous fluxes
  for (m_pertSide = 0; m_pertSide < 2; ++m_pertSide)
  {
    // store the sol pnt normals
    for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
    {
      const CFuint solID = (*(m_states[m_pertSide]))[iState]->getLocalID();
      
      for (CFuint iDim = 0; iDim < m_dim; ++iDim)
      {
        for (CFuint jDim = 0; jDim < m_dim; ++jDim)
        {
          m_neighbCellFluxProjVects[m_pertSide][iDim][iState][jDim] = solPntNormals[solID*m_dim*m_dim+iDim*m_dim+jDim];
        }
      }
    }
    
    // Loop over solution points to calculate the discontinuous flux.
    for (m_pertSol = 0; m_pertSol < m_nbrSolPnts; ++m_pertSol)
    { 
      for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
      {
        *(m_tempGrad[iVar]) = (*(m_cellGrads[m_pertSide][m_pertSol]))[iVar];
      }
      
      const CFuint stateID = (*(m_states[m_pertSide]))[m_pertSol]->getLocalID();
    
      // Set the wall distance before computing the turbulent viscosity
      navierStokesVarSet->setWallDistance(wallDist[stateID]);

      m_avgSol = *((*(m_states[m_pertSide]))[m_pertSol]->getData());
      
      m_updateVarSet->computePhysicalData(*((*(m_states[m_pertSide]))[m_pertSol]), m_pData); 

      prepareFluxComputation();

      // calculate the discontinuous flux projected on x, y, z-directions
      for (CFuint iDim = 0; iDim < m_dim; ++iDim)
      {
        // diffusive part 
        computeFlux(m_avgSol,m_tempGrad,m_neighbCellFluxProjVects[m_pertSide][iDim][m_pertSol],0,m_contFlxBackupDiff[m_pertSide][m_pertSol][iDim]);
        
        m_contFlxBackup[m_pertSide][m_pertSol][iDim] = m_contFlxBackupDiff[m_pertSide][m_pertSol][iDim];
        
        // convective part
        m_contFlxBackup[m_pertSide][m_pertSol][iDim] -= m_updateVarSet->getFlux()(m_pData,m_neighbCellFluxProjVects[m_pertSide][iDim][m_pertSol]);
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffJacobFluxReconstructionTurb::computeCellFluxJacobianNum(const CFreal resFactor)
{
  CFLog(VERBOSE, "computeCellFluxJacobianNum\n");
  
  // Get the wall distance
  DataHandle< CFreal > wallDist = socket_wallDistance.getDataHandle();
  
  SafePtr< NavierStokes2DKLogOmega > navierStokesVarSet = m_diffusiveVarSet.d_castTo< NavierStokes2DKLogOmega >();
    
  // loop over left and right cell to compute flux jacobian
  for (m_pertSide = 0; m_pertSide < 2; ++m_pertSide)
  {
    // variable for the other side
    const CFuint iOtherSide = m_pertSide == LEFT ? RIGHT : LEFT;
    
    // cell ID of the cell at the non-perturbed side
    const CFuint otherCellID = m_cells[iOtherSide]->getID();

    // loop over the states to perturb the states
    for (m_pertSol = 0; m_pertSol < m_nbrSolPnts; ++m_pertSol)
    {
      // dereference state
      State& pertState = *(*m_states[m_pertSide])[m_pertSol];
      
      for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
      {
        *(m_tempGrad[iVar]) = (*(m_cellGrads[m_pertSide][m_pertSol]))[iVar];
      }
      
      const CFuint stateID = (*(m_states[m_pertSide]))[m_pertSol]->getLocalID();
    
      // Set the wall distance before computing the turbulent viscosity
      navierStokesVarSet->setWallDistance(wallDist[stateID]);

      // loop over the variables in the state
      for (m_pertVar = 0; m_pertVar < m_nbrEqs; ++m_pertVar)
      {
        // perturb physical variable in state
        m_numJacob->perturb(m_pertVar,pertState[m_pertVar]);
        
        // compute perturbed physical data
        m_updateVarSet->computePhysicalData(*((*(m_states[m_pertSide]))[m_pertSol]), m_pData); 

        m_avgSol = *((*(m_states[m_pertSide]))[m_pertSol]->getData());
        
        // compute perturbed fluxes
        prepareFluxComputation();

        // calculate the discontinuous flux projected on x, y, z-directions
        for (CFuint iDim = 0; iDim < m_dim; ++iDim)
        {       
          // diffusive part
          computeFlux(m_avgSol,m_tempGrad,m_neighbCellFluxProjVects[m_pertSide][iDim][m_pertSol],0,m_contFlxNeighb[m_pertSide][m_pertSol][iDim]);
                    
          // convective part
          m_contFlxNeighb[m_pertSide][m_pertSol][iDim] -= m_updateVarSet->getFlux()(m_pData,m_neighbCellFluxProjVects[m_pertSide][iDim][m_pertSol]);

          // compute the flux current jacobian term
          // compute the finite difference derivative of the face term (note an implicit minus sign is added here, by the ordering of arguments)
          m_numJacob->computeDerivative(m_contFlxNeighb[m_pertSide][m_pertSol][iDim],m_contFlxBackup[m_pertSide][m_pertSol][iDim],m_tempFlux);

          // multiply residual update derivatives with residual factor so it is taken into the final jacobian
          m_tempFlux *= resFactor;
        
          // store the flux jacobian
          m_fluxJacobian[m_pertSide][m_pertSol][m_pertVar][iDim] = m_tempFlux;
        }

        // restore physical variable in state
        m_numJacob->restore(pertState[m_pertVar]);
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffJacobFluxReconstructionTurb::computeRiemannFluxJacobianNum(const CFreal resFactor)
{
  CFLog(VERBOSE, "Turb computeRiemannFluxJacobianNum\n");
  
  // Get the wall distance
  DataHandle< CFreal > wallDist = socket_wallDistance.getDataHandle();
  
  SafePtr< NavierStokes2DKLogOmega > navierStokesVarSet = m_diffusiveVarSet.d_castTo< NavierStokes2DKLogOmega >();
  
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
      // local flux point indices in the left and right cell
      const CFuint flxPntIdxL = (*m_faceFlxPntConnPerOrient)[m_orient][LEFT][iFlxPnt];
      const CFuint flxPntIdxR = (*m_faceFlxPntConnPerOrient)[m_orient][RIGHT][iFlxPnt];
    
      const CFuint closestSolIdxL = (*m_closestSolToFlxIdx)[flxPntIdxL];
      const CFuint closestSolIdxR = (*m_closestSolToFlxIdx)[flxPntIdxR];
    
      const CFuint stateIDL = (*(m_states[LEFT]))[closestSolIdxL]->getLocalID();
      const CFuint stateIDR = (*(m_states[RIGHT]))[closestSolIdxR]->getLocalID();
    
      const CFreal flxPntWallDist = 0.5*(wallDist[stateIDL]+wallDist[stateIDR]);
    
      // Set the wall distance before computing the turbulent viscosity
      navierStokesVarSet->setWallDistance(flxPntWallDist);
    
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

void ConvDiffJacobFluxReconstructionTurb::computeRiemannFluxToGradJacobianNum(const CFreal resFactor)
{  
  CFLog(VERBOSE, "NS computeRiemannFluxToGradJacobianNum\n");
  
  // Get the wall distance
  DataHandle< CFreal > wallDist = socket_wallDistance.getDataHandle();
  
  SafePtr< NavierStokes2DKLogOmega > navierStokesVarSet = m_diffusiveVarSet.d_castTo< NavierStokes2DKLogOmega >();
  
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_tempStatesL[iFlx] = m_cellStatesFlxPnt[LEFT][iFlx]->getData();
    m_tempStatesR[iFlx] = m_cellStatesFlxPnt[RIGHT][iFlx]->getData();
  }
  
  m_diffusiveVarSetNS->setGradientVars(m_tempStatesL,m_tempGradTermL,m_nbrFaceFlxPnts);
  m_diffusiveVarSetNS->setGradientVars(m_tempStatesR,m_tempGradTermR,m_nbrFaceFlxPnts); 

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
    navierStokesVarSet->setWallDistance(flxPntWallDist);
      
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
  
  // Reset wall distance in case this is the last face: afterwards diff BC is done which needs the wall dist to be zero!!
  navierStokesVarSet->setWallDistance(0.0);
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffJacobFluxReconstructionTurb::computeFluxToGradJacobianNum(const CFreal resFactor)
{
  CFLog(VERBOSE, "computeFluxToGradJacobianNum\n");
  
  // Get the wall distance
  DataHandle< CFreal > wallDist = socket_wallDistance.getDataHandle();
  
  SafePtr< NavierStokes2DKLogOmega > navierStokesVarSet = m_diffusiveVarSet.d_castTo< NavierStokes2DKLogOmega >();
    
  // loop over left and right cell to compute flux jacobian
  for (m_pertSide = 0; m_pertSide < 2; ++m_pertSide)
  {
    // variable for the other side
    const CFuint iOtherSide = m_pertSide == LEFT ? RIGHT : LEFT;
    
    // cell ID of the cell at the non-perturbed side
    const CFuint otherCellID = m_cells[iOtherSide]->getID();

    // loop over the states to perturb the states
    for (m_pertSol = 0; m_pertSol < m_nbrSolPnts; ++m_pertSol)
    {
      const CFuint stateID = (*(m_states[m_pertSide]))[m_pertSol]->getLocalID();
    
      // Set the wall distance before computing the turbulent viscosity
      navierStokesVarSet->setWallDistance(wallDist[stateID]);
      
      // dereference state
      m_avgSol = *((*(m_states[m_pertSide]))[m_pertSol]->getData());
      
      for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
      {
        *(m_tempGrad[iVar]) = (*(m_cellGrads[m_pertSide][m_pertSol]))[iVar];
      }

      // loop over the variables in the state
      for (m_pertVar = 0; m_pertVar < m_nbrEqs; ++m_pertVar)
      {
        for (CFuint pertDim = 0; pertDim < m_dim; ++pertDim)
        {
          // perturb physical variable in state
          m_numJacob->perturb(m_pertVar,(*(m_tempGrad[m_pertVar]))[pertDim]);

          // compute perturbed fluxes
          prepareFluxComputation();

          // calculate the discontinuous flux projected on x, y, z-directions
          for (CFuint iDim = 0; iDim < m_dim; ++iDim)
          {       
            // diffusive part
            computeFlux(m_avgSol,m_tempGrad,m_neighbCellFluxProjVects[m_pertSide][iDim][m_pertSol],0,m_contFlxNeighb[m_pertSide][m_pertSol][iDim]);

            // compute the flux current jacobian term
            // compute the finite difference derivative of the face term (note an implicit minus sign is added here, by the ordering of arguments)
            m_numJacob->computeDerivative(m_contFlxNeighb[m_pertSide][m_pertSol][iDim],m_contFlxBackupDiff[m_pertSide][m_pertSol][iDim],m_tempFlux);

            // multiply residual update derivatives with residual factor so it is taken into the final jacobian
            m_tempFlux *= resFactor;

            // store the flux jacobian
            m_gradientFluxJacobian[m_pertSide][m_pertSol][m_pertVar][pertDim][iDim] = m_tempFlux;
          }

          // restore physical variable in state
          m_numJacob->restore((*(m_tempGrad[m_pertVar]))[pertDim]);
        }
      }
    }
  }  
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffJacobFluxReconstructionTurb::setup()
{
  CFAUTOTRACE;
  
  // setup parent class
  ConvDiffJacobFluxReconstructionNS::setup();
  
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
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD
