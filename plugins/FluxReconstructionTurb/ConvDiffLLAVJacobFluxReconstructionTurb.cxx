// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/MethodCommandProvider.hh"

#include "Framework/CFSide.hh"
#include "Framework/MeshData.hh"
#include "Framework/BaseTerm.hh"
#include "Framework/BlockAccumulator.hh"
#include "Framework/LSSMatrix.hh"

#include "MathTools/MathFunctions.hh"

#include "FluxReconstructionTurb/ConvDiffLLAVJacobFluxReconstructionTurb.hh"
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

MethodCommandProvider< ConvDiffLLAVJacobFluxReconstructionTurb,
		       FluxReconstructionSolverData,
		       FluxReconstructionNavierStokesModule >
convDiffLLAVRHSJacobTurbFluxReconstructionProvider("ConvDiffLLAVRHSJacobTurb");
  
//////////////////////////////////////////////////////////////////////////////
  
ConvDiffLLAVJacobFluxReconstructionTurb::ConvDiffLLAVJacobFluxReconstructionTurb(const std::string& name) :
  ConvDiffLLAVJacobFluxReconstructionNS(name),
  socket_wallDistance("wallDistance"),
  m_closestSolToFlxIdx(CFNULL),
  m_pData(),
  m_pData2()
{
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffLLAVJacobFluxReconstructionTurb::configure ( Config::ConfigArgs& args )
{
  ConvDiffLLAVJacobFluxReconstructionNS::configure(args);
} 

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > >
    ConvDiffLLAVJacobFluxReconstructionTurb::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result = ConvDiffLLAVJacobFluxReconstructionNS::needsSockets();

  result.push_back(&socket_wallDistance);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffLLAVJacobFluxReconstructionTurb::computeInterfaceFlxCorrection()
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
    computeFlux(m_avgSol,m_avgGrad,m_unitNormalFlxPnts[iFlxPnt],0,m_flxPntRiemannFluxDiff[iFlxPnt]);
    
    m_flxPntRiemannFluxDiffConv[iFlxPnt] = m_flxPntRiemannFluxDiff[iFlxPnt];
    
    // compute the convective riemann flux
    m_flxPntRiemannFluxDiffConv[iFlxPnt] -= m_riemannFluxComputer->computeFlux(*(m_cellStatesFlxPnt[LEFT][iFlxPnt]),
									    *(m_cellStatesFlxPnt[RIGHT][iFlxPnt]),
									    m_unitNormalFlxPnts[iFlxPnt]);
    
    m_flxPntRiemannFlux[iFlxPnt] = m_flxPntRiemannFluxDiffConv[iFlxPnt];
     
    // compute artificial part
    // get epsilon
    const CFreal epsilon = 0.5*(m_epsilonLR[LEFT][iFlxPnt]+m_epsilonLR[RIGHT][iFlxPnt]);
    
    for (CFuint iDim = 0; iDim < m_dim; ++iDim)
    {
      for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
      {
        m_flxPntRiemannFlux[iFlxPnt][iVar] += epsilon*((*(m_avgGradAV[iVar]))[iDim])*m_unitNormalFlxPnts[iFlxPnt][iDim];
      }
    }
     
    // compute FI in the mapped coord frame
    m_cellFlx[LEFT][iFlxPnt] = (m_flxPntRiemannFlux[iFlxPnt])*m_faceJacobVecSizeFlxPnts[iFlxPnt][LEFT];
    m_cellFlx[RIGHT][iFlxPnt] = (m_flxPntRiemannFlux[iFlxPnt])*m_faceJacobVecSizeFlxPnts[iFlxPnt][RIGHT];
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffLLAVJacobFluxReconstructionTurb::computeUnpertCellDiffResiduals(const CFuint side)
{ 
  // get datahandle
  DataHandle< CFreal > artVisc = socket_artVisc.getDataHandle();
  
  //SafePtr< NavierStokes2DKLogOmega > navierStokesVarSet = m_diffusiveVarSet.d_castTo< NavierStokes2DKLogOmega >();
  
  // Get the wall distance
  DataHandle< CFreal > wallDist = socket_wallDistance.getDataHandle();
  
  // loop over flx pnts to extrapolate the states to the flux points
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {     
    artVisc[(((*(m_states[side]))[iSol]))->getLocalID()] = m_solEpsilons[side][iSol];
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
      computeFlux(m_avgSol,m_tempGrad,m_neighbCellFluxProjVects[side][iDim][iSolPnt],0,m_contFlxWoLLAV[iSolPnt][iDim]);
      
      // add convective part
      m_contFlxWoLLAV[iSolPnt][iDim] -= m_updateVarSet->getFlux()(m_pData,m_neighbCellFluxProjVects[side][iDim][iSolPnt]);
      
      m_contFlx[iSolPnt][iDim] = m_contFlxWoLLAV[iSolPnt][iDim];
      
      // add artificial part
      for (CFuint iDim2 = 0; iDim2 < m_dim; ++iDim2)
      {
        for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
        {
          m_contFlx[iSolPnt][iDim][iVar] += m_solEpsilons[side][iSolPnt]*((*m_avgGradAV[iVar])[iDim2])*m_neighbCellFluxProjVects[side][iDim][iSolPnt][iDim2];
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
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffLLAVJacobFluxReconstructionTurb::initJacobianComputation()
{
  CFLog(VERBOSE, "initJacobianComputation\n");
  
  // Get the wall distance
  DataHandle< CFreal > wallDist = socket_wallDistance.getDataHandle();
  
  //SafePtr< NavierStokes2DKLogOmega > navierStokesVarSet = m_diffusiveVarSet.d_castTo< NavierStokes2DKLogOmega >();
  
  // set block row and column indices, proj vectors and make a backup of discontinuous fluxes
  for (m_pertSide = 0; m_pertSide < 2; ++m_pertSide)
  {    
    // Loop over solution points to calculate the discontinuous flux.
    for (m_pertSol = 0; m_pertSol < m_nbrSolPnts; ++m_pertSol)
    { 
      for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
      {
        *(m_tempGrad[iVar]) = (*(m_cellGrads[m_pertSide][m_pertSol]))[iVar];
      }
      
      const CFuint stateID = (*(m_states[m_pertSide]))[m_pertSol]->getLocalID();
    
      // Set the wall distance before computing the turbulent viscosity
      if (m_dim == 2)
      {
        m_navierStokesVarSetTurb->setWallDistance(wallDist[stateID]);
      }
      else
      {
        m_navierStokesVarSetTurb3D->setWallDistance(wallDist[stateID]); 
      }

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

void ConvDiffLLAVJacobFluxReconstructionTurb::computeCellFluxJacobianNum(const CFreal resFactor)
{
  CFLog(VERBOSE, "computeCellFluxJacobianNum\n");
  
  // Get the wall distance
  DataHandle< CFreal > wallDist = socket_wallDistance.getDataHandle();
  
  //SafePtr< NavierStokes2DKLogOmega > navierStokesVarSet = m_diffusiveVarSet.d_castTo< NavierStokes2DKLogOmega >();
  
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
      if (m_dim == 2)
      {
        m_navierStokesVarSetTurb->setWallDistance(wallDist[stateID]);
      }
      else
      {
        m_navierStokesVarSetTurb3D->setWallDistance(wallDist[stateID]); 
      }

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

void ConvDiffLLAVJacobFluxReconstructionTurb::computeRiemannFluxJacobianNum(const CFreal resFactor)
{
  CFLog(VERBOSE, "Turb computeRiemannFluxJacobianNum\n");
  
  // Get the wall distance
  DataHandle< CFreal > wallDist = socket_wallDistance.getDataHandle();
  
  //SafePtr< NavierStokes2DKLogOmega > navierStokesVarSet = m_diffusiveVarSet.d_castTo< NavierStokes2DKLogOmega >();
  
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
      if (m_dim == 2)
      {
        m_navierStokesVarSetTurb->setWallDistance(flxPntWallDist);
      }
      else
      {
        m_navierStokesVarSetTurb3D->setWallDistance(flxPntWallDist); 
      }
      
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
        m_numJacob->computeDerivative(m_flxPntRiemannFluxPert[iFlxPnt],m_flxPntRiemannFluxDiffConv[iFlxPnt],m_riemannFluxJacobian[m_pertSide][iFlxPnt][m_pertVar]);

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

void ConvDiffLLAVJacobFluxReconstructionTurb::computeRiemannFluxToGradJacobianNum(const CFreal resFactor)
{  
  CFLog(VERBOSE, "NS computeRiemannFluxToGradJacobianNum\n");
  
  // Get the wall distance
  DataHandle< CFreal > wallDist = socket_wallDistance.getDataHandle();
  
  //SafePtr< NavierStokes2DKLogOmega > navierStokesVarSet = m_diffusiveVarSet.d_castTo< NavierStokes2DKLogOmega >();
  
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
    if (m_dim == 2)
    {
      m_navierStokesVarSetTurb->setWallDistance(flxPntWallDist);
    }
    else
    {
      m_navierStokesVarSetTurb3D->setWallDistance(flxPntWallDist); 
    }
    
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

void ConvDiffLLAVJacobFluxReconstructionTurb::computeFluxToGradJacobianNum(const CFreal resFactor)
{
  CFLog(VERBOSE, "computeFluxToGradJacobianNum\n");
  
  // Get the wall distance
  DataHandle< CFreal > wallDist = socket_wallDistance.getDataHandle();
  
  //SafePtr< NavierStokes2DKLogOmega > navierStokesVarSet = m_diffusiveVarSet.d_castTo< NavierStokes2DKLogOmega >();
  
  // loop over left and right cell to compute flux jacobian
  for (m_pertSide = 0; m_pertSide < 2; ++m_pertSide)
  {
    // variable for the other side
    const CFuint iOtherSide = m_pertSide == LEFT ? RIGHT : LEFT;

    // loop over the states to perturb the states
    for (m_pertSol = 0; m_pertSol < m_nbrSolPnts; ++m_pertSol)
    {
      const CFuint stateID = (*(m_states[m_pertSide]))[m_pertSol]->getLocalID();
    
      // Set the wall distance before computing the turbulent viscosity
      if (m_dim == 2)
      {
        m_navierStokesVarSetTurb->setWallDistance(wallDist[stateID]);
      }
      else
      {
        m_navierStokesVarSetTurb3D->setWallDistance(wallDist[stateID]); 
      }
      
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

void ConvDiffLLAVJacobFluxReconstructionTurb::computeBothJacobsDiffFaceTerm()
{
  CFLog(VERBOSE, "computeBothJacobsDiffFaceTerm\n");
    
  // get residual factor
  const CFreal resFactor = getMethodData().getResFactor();

  // dereference accumulator
  BlockAccumulator& acc = *m_acc;
  
  CFuint solIdx = 0;
  for (m_pertSide = 0; m_pertSide < 2; ++m_pertSide)
  {
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol, ++solIdx)
    {
      acc.setRowColIndex(solIdx,(*m_states[m_pertSide])[iSol]->getLocalID());
    }
  }

  //// compute the needed flux jacobians
  
  initJacobianComputation();
  
  computeCellFluxJacobianNum(resFactor);
  
  computeRiemannFluxJacobianNum(resFactor);
  
  computeFluxToGradJacobianNum(resFactor);
  
  if (m_addRiemannToGradJacob || m_addRiemannToGradCrossCellJacob) computeRiemannFluxToGradJacobianNum(resFactor);
  
  computeGradToStateJacobianAna();
  
  computeGradVarsToStateJacobianNum();
  
  computeEpsToStateJacobianAna();
  
  computeLLAVCellFluxJacobianAna(resFactor);
  
  computeLLAVRiemannFluxJacobianAna(resFactor);
  
  //// add the total jacobians to the system jacobian
  
  // loop over left and right cell to add the discontinuous (cell) part to the jacobian
  for (m_pertSide = 0; m_pertSide < 2; ++m_pertSide)
  {
    // make sure this is only done once per cell
    if (!m_cellFlags[m_cells[m_pertSide]->getID()]) 
    {
      // term depending on iSide
      const CFuint pertSideTerm = m_pertSide*m_nbrSolPnts;

      // loop over the states to which to derive (l)
      for (m_pertSol = 0; m_pertSol < m_nbrSolPnts; ++m_pertSol)
      {
        // loop over the variables in the state (k)
        for (m_pertVar = 0; m_pertVar < m_nbrEqs; ++m_pertVar)
        {
          const CFuint nbDepGradVar = m_nbrVarToGradVarDep[m_pertVar];
            
          // add the discontinuous part of the jacobian related to the sol pnt (i)
          for (CFuint jSolPnt = 0; jSolPnt < m_nbrSolSolDep; ++jSolPnt)
          {
            const CFuint jSolIdx = (*m_solSolDep)[m_pertSol][jSolPnt];
            
            m_tempFlux = 0.;

            // (d)
            for (CFuint iDim = 0; iDim < m_dim; ++iDim)
            {
              const CFreal polyCoef = (*m_solPolyDerivAtSolPnts)[jSolIdx][iDim][m_pertSol]; 
          
              m_tempFlux += m_fluxJacobian[m_pertSide][m_pertSol][m_pertVar][iDim] * polyCoef;
            }
            
            acc.addValues(jSolIdx+pertSideTerm,m_pertSol+pertSideTerm,m_pertVar,&m_tempFlux[0]);
          }
          
          // add the discontinuous gradient part of the jacobian (m)
          for (CFuint kSolPnt = 0; kSolPnt < m_nbrSolSolDep; ++kSolPnt)
          {
            const CFuint kSolIdx = (*m_solSolDep)[m_pertSol][kSolPnt];
          
            // (i)
            for (CFuint jSol = 0; jSol < m_nbrSolSolDep; ++jSol)
            {
              const CFuint jSolIdx = (*m_solSolDep)[kSolIdx][jSol];
              
              m_tempFlux = 0.0;
                
              // (d)
              for (CFuint iDim = 0; iDim < m_dim; ++iDim)
              {
                const CFreal dl = (*m_solPolyDerivAtSolPnts)[jSolIdx][iDim][kSolIdx];
                
                /// llav jacob to state part //// should actually be added for all jSol
                m_tempFlux += m_llavFluxJacobian[m_pertSide][m_pertSol][m_pertVar][m_pertSide][jSolIdx][iDim] * dl;
                
                // (b)
                for (CFuint jDim = 0; jDim < m_dim; ++jDim)
                {
                  const CFreal dl_dqdu = dl * m_gradientStateJacobian[m_pertSide][kSolIdx][m_pertSide][m_pertSol][jDim];
                    
                  // (p)
                  for (CFuint iEq = 0; iEq < nbDepGradVar; ++iEq)
                  {
                    const CFuint var = m_varToGradVarDep[m_pertVar][iEq];
                    
                    const CFreal dl_dqdu_dudu = m_gradVarsToStateJacobian[m_pertSide][m_pertSol][m_pertVar][var] * dl_dqdu;
                      
                    m_tempFlux += m_gradientFluxJacobian[m_pertSide][kSolIdx][var][jDim][iDim] * dl_dqdu_dudu;
                  }
                  
//                  CFreal llavPart = m_solEpsilons[m_pertSide][kSolIdx] * m_neighbCellFluxProjVects[m_pertSide][iDim][kSolIdx][jDim];
//                  llavPart *= dl_dqdu;
//                  //if(m_cells[0]->getID()==1) CFLog(INFO, "pertSol: " << m_pertSol << ", pertVar: " << m_pertVar << ", llavJC: " << llavPart << "\n");
//                  // add part of analytical LLAV jacobian
//                  m_tempFlux[m_pertVar] += llavPart;
                }
              }
            
              acc.addValues(jSolIdx+pertSideTerm,m_pertSol+pertSideTerm,m_pertVar,&m_tempFlux[0]);
            }
          }
          
          // add the discontinuous part of the jacobian related to the flx pnt (f)
          for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFlxDep; ++iFlxPnt)
          {
            const CFuint flxIdx = (*m_solFlxDep)[m_pertSol][iFlxPnt];
            
            // (df)
            const CFuint dim = (*m_flxPntFlxDim)[flxIdx];
            
            m_temp = m_fluxJacobian[m_pertSide][m_pertSol][m_pertVar][dim] * (*m_solPolyValsAtFlxPnts)[flxIdx][m_pertSol];
            
            // add the second part of the discontinuous part of the jacobian (i)
            for (CFuint jSolPnt = 0; jSolPnt < m_nbrSolDep; ++jSolPnt)
            {
              const CFuint jSolIdx = (*m_flxSolDep)[flxIdx][jSolPnt];

              // get the divergence of the correction function
              const CFreal divh = m_corrFctDiv[jSolIdx][flxIdx];
                           
              m_tempFlux = -m_temp * divh;
            
              acc.addValues(jSolIdx+pertSideTerm,m_pertSol+pertSideTerm,m_pertVar,&m_tempFlux[0]);
            }
          }
          
          // add the second part of the discontinuous gradient part of the jacobian (m)
          for (CFuint kSolPnt = 0; kSolPnt < m_nbrSolSolDep; ++kSolPnt)
          {
            const CFuint kSolIdx = (*m_solSolDep)[m_pertSol][kSolPnt];
              
            for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFlxDep; ++iFlxPnt)
            {   
              const CFuint flxIdx = (*m_solFlxDep)[kSolIdx][iFlxPnt];
              
              // (df)
              const CFuint dim = (*m_flxPntFlxDim)[flxIdx];
              
              const CFreal l = (*m_solPolyValsAtFlxPnts)[flxIdx][kSolIdx];            
              
              // (i)
              for (CFuint jSolPnt = 0; jSolPnt < m_nbrSolDep; ++jSolPnt)
              {
                const CFuint jSolIdx = (*m_flxSolDep)[flxIdx][jSolPnt];
              
                m_tempFlux = 0.0;

                const CFreal divh_l = -m_corrFctDiv[jSolIdx][flxIdx] * l;
                
                /// llav jacob to state part //// actually should loop over all kSol
                m_tempFlux += m_llavFluxJacobian[m_pertSide][m_pertSol][m_pertVar][m_pertSide][jSolIdx][dim] * divh_l;
              
                // (b)
                for (CFuint jDim = 0; jDim < m_dim; ++jDim)
                {
                  const CFreal divh_l_dqdu = divh_l * m_gradientStateJacobian[m_pertSide][kSolIdx][m_pertSide][m_pertSol][jDim];
                  
                  // (p)
                  for (CFuint iEq = 0; iEq < nbDepGradVar; ++iEq)
                  {
                    const CFuint var = m_varToGradVarDep[m_pertVar][iEq];
                    
                    const CFreal divh_l_dqdu_dudu = divh_l_dqdu * m_gradVarsToStateJacobian[m_pertSide][m_pertSol][m_pertVar][var];             
                      
                    m_tempFlux += divh_l_dqdu_dudu * m_gradientFluxJacobian[m_pertSide][kSolIdx][var][jDim][dim];
                  }
                  
//                  CFreal llavPart = divh_l_dqdu * m_solEpsilons[m_pertSide][kSolIdx];
//                  
//                  llavPart *= m_neighbCellFluxProjVects[m_pertSide][dim][kSolIdx][jDim];
//                  //if(m_cells[0]->getID()==1) CFLog(INFO, "pertSol: " << m_pertSol << ", pertVar: " << m_pertVar << ", llavJF: " << llavPart << "\n");
//                  // add part of analytical LLAV Jacobian
//                  m_tempFlux[m_pertVar] += llavPart;
                }
                
                acc.addValues(jSolIdx+pertSideTerm,m_pertSol+pertSideTerm,m_pertVar,&m_tempFlux[0]);
              }
            }
          }
        }
      }
    }
  }
  
  // loop over left and right cell to add the riemann flux (face) part to the jacobian
  for (m_pertSide = 0; m_pertSide < 2; ++m_pertSide)
  {
    // variable for the other side
    const CFuint iOtherSide = m_pertSide == LEFT ? RIGHT : LEFT;
    
    // term depending on iSide
    const CFuint pertSideTerm = m_pertSide*m_nbrSolPnts;

    // term depending on iOtherSide
    const CFuint otherSideTerm = iOtherSide*m_nbrSolPnts;
    
    // loop over the variables in the state (k)
    for (m_pertVar = 0; m_pertVar < m_nbrEqs; ++m_pertVar)
    { 
      const CFuint nbDepGradVar = m_nbrVarToGradVarDep[m_pertVar];
        
      // loop over face flx pnts (f)
      for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
      {
        // local flux point indices in the left and right cell
        const CFuint flxPntIdxThis = (*m_faceFlxPntConnPerOrient)[m_orient][m_pertSide][iFlxPnt];
        const CFuint flxPntIdxOther = (*m_faceFlxPntConnPerOrient)[m_orient][iOtherSide][iFlxPnt];
        
        m_temp = m_riemannFluxJacobian[m_pertSide][iFlxPnt][m_pertVar]*m_faceJacobVecSizeFlxPnts[iFlxPnt][m_pertSide];
        m_tempOther = m_riemannFluxJacobian[m_pertSide][iFlxPnt][m_pertVar]*m_faceJacobVecSizeFlxPnts[iFlxPnt][iOtherSide];
        
        const CFreal halfFaceJacob = 0.5 * m_faceJacobVecSizeFlxPnts[iFlxPnt][m_pertSide];
        
        // loop over the states to perturb the states (l)
        for (m_pertSol = 0; m_pertSol < m_nbrSolDep; ++m_pertSol)
        {
          const CFuint pertSolIdx = (*m_flxSolDep)[flxPntIdxThis][m_pertSol];
          
          m_temp2 = m_temp * (*m_solPolyValsAtFlxPnts)[flxPntIdxThis][pertSolIdx];
          m_tempOther2 = m_tempOther * (*m_solPolyValsAtFlxPnts)[flxPntIdxThis][pertSolIdx];
            
          // add the second part of the discontinuous part of the jacobian (i)
          for (CFuint jSolPnt = 0; jSolPnt < m_nbrSolDep; ++jSolPnt)
          {
            const CFuint jSolIdxThis = (*m_flxSolDep)[flxPntIdxThis][jSolPnt];
            const CFuint jSolIdxOther = (*m_flxSolDep)[flxPntIdxOther][jSolPnt];

            // get the divergence of the correction function on this side
            CFreal divh = m_corrFctDiv[jSolIdxThis][flxPntIdxThis];
                          
            // add part on this side of face
            m_tempFlux = m_temp2 * divh;
              
            acc.addValues(jSolIdxThis+pertSideTerm,pertSolIdx+pertSideTerm,m_pertVar,&m_tempFlux[0]);
          
            // get the divergence of the correction function on other side
            divh = m_corrFctDiv[jSolIdxOther][flxPntIdxOther];
                        
            // add cross-cell part 
            m_tempFlux = m_tempOther2 * divh;   
              
            acc.addValues(jSolIdxOther+otherSideTerm,pertSolIdx+pertSideTerm,m_pertVar,&m_tempFlux[0]);
          }
        }
        
        // loop over the states to perturb the states (l) for LLAV to state part
        for (m_pertSol = 0; m_pertSol < m_nbrSolPnts; ++m_pertSol)
        { 
          m_temp2 = m_llavRiemannFluxJacobian[m_pertSide][m_pertSol][m_pertVar][iFlxPnt] * m_faceJacobVecSizeFlxPnts[iFlxPnt][m_pertSide];
          m_tempOther2 = m_llavRiemannFluxJacobian[m_pertSide][m_pertSol][m_pertVar][iFlxPnt] * m_faceJacobVecSizeFlxPnts[iFlxPnt][iOtherSide];
          
          // add the LLAV interface part of the jacobian (i)
          for (CFuint jSolPnt = 0; jSolPnt < m_nbrSolDep; ++jSolPnt)
          {
            const CFuint jSolIdxThis = (*m_flxSolDep)[flxPntIdxThis][jSolPnt];
            const CFuint jSolIdxOther = (*m_flxSolDep)[flxPntIdxOther][jSolPnt];

            // get the divergence of the correction function on this side
            CFreal divh = m_corrFctDiv[jSolIdxThis][flxPntIdxThis];
                          
            // add part on this side of face
            m_tempFlux = m_temp2 * divh;
              
            acc.addValues(jSolIdxThis+pertSideTerm,m_pertSol+pertSideTerm,m_pertVar,&m_tempFlux[0]);
          
            // get the divergence of the correction function on other side
            divh = m_corrFctDiv[jSolIdxOther][flxPntIdxOther];
                        
            // add cross-cell part 
            m_tempFlux = m_tempOther2 * divh;   
              
            acc.addValues(jSolIdxOther+otherSideTerm,m_pertSol+pertSideTerm,m_pertVar,&m_tempFlux[0]);
          }
        }
        
        for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
        { 
          m_needToAddSolPnt[iSol] = true;
        }
            
        // (i)
        for (CFuint jSolPnt = 0; jSolPnt < m_nbrSolDep; ++jSolPnt)
        {
          const CFuint jSolIdxThis = (*m_flxSolDep)[flxPntIdxThis][jSolPnt];
          const CFuint jSolIdxOther = (*m_flxSolDep)[flxPntIdxOther][jSolPnt];

          // get the divergence of the correction function on this side
          CFreal divh = m_corrFctDiv[jSolIdxThis][flxPntIdxThis];
            
          const CFreal divh_halfFaceJacob = divh * halfFaceJacob;
          
          // loop over the states to perturb the states (l)
          for (CFuint m_pertSol = 0; m_pertSol < m_nbrSolDep; ++m_pertSol)
          {   
            const CFuint pertSolIdx = (*m_flxSolDep)[flxPntIdxThis][m_pertSol];
            
            m_needToAddSolPnt[pertSolIdx] = false;
            
            if (m_addRiemannToGradJacob)
            {
            // add part on this side of face
            m_tempFlux = 0.0;

            // (m)
            for (CFuint kSolPnt = 0; kSolPnt < m_nbrSolDep; ++kSolPnt)
            {
              const CFuint kSolIdx = (*m_flxSolDep)[flxPntIdxThis][kSolPnt];
              const CFuint kSolIdxOther = (*m_flxSolDep)[flxPntIdxOther][kSolPnt];
              
              const CFreal divh_halfFaceJacob_l = divh_halfFaceJacob * (*m_solPolyValsAtFlxPnts)[flxPntIdxThis][kSolIdx];
              const CFreal divh_halfFaceJacob_lOther = divh_halfFaceJacob * (*m_solPolyValsAtFlxPnts)[flxPntIdxOther][kSolIdxOther];
              
              /// llav jacob to state part
              //m_tempFlux += m_llavRiemannFluxJacobian[m_pertSide][pertSolIdx][m_pertVar][iFlxPnt] * divh_halfFaceJacob_l;
              //m_tempFlux += m_llavRiemannFluxJacobian[m_pertSide][pertSolIdx][m_pertVar][iFlxPnt] * divh_halfFaceJacob_lOther;
                
              // (b)
              for (CFuint iDim = 0; iDim < m_dim; ++iDim)
              {
                const CFreal divh_halfFaceJacob_l_dqdu = divh_halfFaceJacob_l * m_gradientStateJacobian[m_pertSide][kSolIdx][m_pertSide][pertSolIdx][iDim];
                const CFreal divh_halfFaceJacob_l_dqduOther = divh_halfFaceJacob_lOther * m_gradientStateJacobian[iOtherSide][kSolIdxOther][m_pertSide][pertSolIdx][iDim];
                
                // (p)
                for (CFuint iEq = 0; iEq < nbDepGradVar; ++iEq)
                {
                  const CFuint var = m_varToGradVarDep[m_pertVar][iEq];
                  
                  m_temp =  m_riemannFluxGradJacobian[iFlxPnt][var][iDim] * m_gradVarsToStateJacobian[m_pertSide][pertSolIdx][m_pertVar][var];
                    
                  m_tempFlux += m_temp * divh_halfFaceJacob_l_dqdu;
                    
                  m_tempFlux += m_temp * divh_halfFaceJacob_l_dqduOther;
                }
                
//                const CFreal epsilon = 0.5*(m_epsilonLR[LEFT][iFlxPnt]+m_epsilonLR[RIGHT][iFlxPnt]);
//                
//                const CFreal llavPart = epsilon * m_unitNormalFlxPnts[iFlxPnt][iDim];
//                
//                // add part of analytical LLAV jacobian
//                m_tempFlux[m_pertVar] += llavPart * divh_halfFaceJacob_l_dqdu;    
//                m_tempFlux[m_pertVar] += llavPart * divh_halfFaceJacob_l_dqduOther;  
              }         
            }
              
            acc.addValues(jSolIdxThis+pertSideTerm,pertSolIdx+pertSideTerm,m_pertVar,&m_tempFlux[0]);
          
            // get the divergence of the correction function on other side
            divh = m_corrFctDiv[jSolIdxOther][flxPntIdxOther];
            
            const CFreal divh_halfFaceJacobOther = 0.5 * divh * m_faceJacobVecSizeFlxPnts[iFlxPnt][iOtherSide];

            // add cross-cell part 
            m_tempFlux = 0.0;

            // (m)
            for (CFuint kSolPnt = 0; kSolPnt < m_nbrSolDep; ++kSolPnt)
            {
              const CFuint kSolIdx = (*m_flxSolDep)[flxPntIdxThis][kSolPnt];
              const CFuint kSolIdxOther = (*m_flxSolDep)[flxPntIdxOther][kSolPnt];
              
              const CFreal divh_halfFaceJacobOther_lThis = divh_halfFaceJacobOther * (*m_solPolyValsAtFlxPnts)[flxPntIdxThis][kSolIdx];
              const CFreal divh_halfFaceJacobOther_lOther = divh_halfFaceJacobOther * (*m_solPolyValsAtFlxPnts)[flxPntIdxOther][kSolIdxOther];
              
              /// llav jacob to state part
              //m_tempFlux += m_llavRiemannFluxJacobian[m_pertSide][pertSolIdx][m_pertVar][iFlxPnt] * divh_halfFaceJacobOther_lThis;
              //m_tempFlux += m_llavRiemannFluxJacobian[m_pertSide][pertSolIdx][m_pertVar][iFlxPnt] * divh_halfFaceJacobOther_lOther;
                
              // (b)
              for (CFuint iDim = 0; iDim < m_dim; ++iDim)
              {
                const CFreal divh_halfFaceJacobOther_lThis_dqduThis = divh_halfFaceJacobOther_lThis * m_gradientStateJacobian[m_pertSide][kSolIdx][m_pertSide][pertSolIdx][iDim];
                const CFreal divh_halfFaceJacobOther_lOther_dqduOther = divh_halfFaceJacobOther_lOther * m_gradientStateJacobian[iOtherSide][kSolIdxOther][m_pertSide][pertSolIdx][iDim];
              
                // (p)
                for (CFuint iEq = 0; iEq < nbDepGradVar; ++iEq)
                {
                  const CFuint var = m_varToGradVarDep[m_pertVar][iEq];
                  
                  m_temp = m_riemannFluxGradJacobian[iFlxPnt][var][iDim] * m_gradVarsToStateJacobian[m_pertSide][pertSolIdx][m_pertVar][var];
                    
                  m_tempFlux += m_temp * divh_halfFaceJacobOther_lThis_dqduThis;
                    
                  m_tempFlux += m_temp * divh_halfFaceJacobOther_lOther_dqduOther;
                }
                
//                const CFreal epsilon = 0.5*(m_epsilonLR[LEFT][iFlxPnt]+m_epsilonLR[RIGHT][iFlxPnt]);
//                
//                const CFreal llavPart = epsilon * m_unitNormalFlxPnts[iFlxPnt][iDim];
//                
//                // add part of analytical LLAV jacobian
//                m_tempFlux[m_pertVar] += llavPart * divh_halfFaceJacobOther_lThis_dqduThis;    
//                m_tempFlux[m_pertVar] += llavPart * divh_halfFaceJacobOther_lOther_dqduOther;
              }         
            }
            
            acc.addValues(jSolIdxOther+otherSideTerm,pertSolIdx+pertSideTerm,m_pertVar,&m_tempFlux[0]);
          }
          }
          
          if (m_addRiemannToGradCrossCellJacob)
          {
          // loop over the states to perturb the states (l)
          for (CFuint pertSolIdx = 0; pertSolIdx < m_nbrSolPnts; ++pertSolIdx)
          {
            CFuint dependingKSol = 1000;
              
            if (m_needToAddSolPnt[pertSolIdx])
            {
              // add part on this side of face
              m_tempFlux = 0.0;

              // (m)
              for (CFuint kSolPnt = 0; kSolPnt < m_nbrSolDep; ++kSolPnt)
              {
                const CFuint kSolIdx = (*m_flxSolDep)[flxPntIdxThis][kSolPnt];
              
                for (CFuint lSol = 0; lSol < m_nbrSolSolDep; ++lSol)
                {
                  const CFuint lSolIdx = (*m_solSolDep)[pertSolIdx][lSol]; 
                
                  if (lSolIdx == kSolIdx)
                  {
                    dependingKSol = kSolIdx;
                    break;
                  }
                }
              }
  
              const CFreal divh_halfFaceJacob_l = divh_halfFaceJacob * (*m_solPolyValsAtFlxPnts)[flxPntIdxThis][dependingKSol];
              
              /// llav jacob to state part
              //m_tempFlux += m_llavRiemannFluxJacobian[m_pertSide][pertSolIdx][m_pertVar][iFlxPnt] * divh_halfFaceJacob_l;
                
              // (b)
              for (CFuint iDim = 0; iDim < m_dim; ++iDim)
              {
                const CFreal divh_halfFaceJacob_l_dqdu = divh_halfFaceJacob_l * m_gradientStateJacobian[m_pertSide][dependingKSol][m_pertSide][pertSolIdx][iDim];
                
                // (p)
                for (CFuint iEq = 0; iEq < nbDepGradVar; ++iEq)
                {
                  const CFuint var = m_varToGradVarDep[m_pertVar][iEq];
                  
                  const CFreal divh_halfFaceJacob_l_dqdu_dudu = m_gradVarsToStateJacobian[m_pertSide][pertSolIdx][m_pertVar][var] * divh_halfFaceJacob_l_dqdu;
                    
                  m_tempFlux += m_riemannFluxGradJacobian[iFlxPnt][var][iDim] * divh_halfFaceJacob_l_dqdu_dudu;
                }
                
//                const CFreal epsilon = 0.5*(m_epsilonLR[LEFT][iFlxPnt]+m_epsilonLR[RIGHT][iFlxPnt]);
//                
//                CFreal llavPart =  epsilon * m_unitNormalFlxPnts[iFlxPnt][iDim];
//                
//                llavPart *= divh_halfFaceJacob_l_dqdu;
//                
//                // add part of analytical LLAV jacobian
//                m_tempFlux[m_pertVar] += llavPart;     
              } 
              
              acc.addValues(jSolIdxThis+pertSideTerm,pertSolIdx+pertSideTerm,m_pertVar,&m_tempFlux[0]);
              
          
              // get the divergence of the correction function on other side
              divh = m_corrFctDiv[jSolIdxOther][flxPntIdxOther];
            
              const CFreal divh_halfFaceJacobOther = 0.5 * divh * m_faceJacobVecSizeFlxPnts[iFlxPnt][iOtherSide];

              // add cross-cell part 
              m_tempFlux = 0.0;
              
              const CFreal divh_halfFaceJacobOther_lThis = divh_halfFaceJacobOther * (*m_solPolyValsAtFlxPnts)[flxPntIdxThis][dependingKSol];
                
              /// llav jacob to state part
              //m_tempFlux += m_llavRiemannFluxJacobian[m_pertSide][pertSolIdx][m_pertVar][iFlxPnt] * divh_halfFaceJacobOther_lThis;
              
              // (b)
              for (CFuint iDim = 0; iDim < m_dim; ++iDim)
              {
                const CFreal divh_halfFaceJacobOther_lThis_dqduThis = divh_halfFaceJacobOther_lThis * m_gradientStateJacobian[m_pertSide][dependingKSol][m_pertSide][pertSolIdx][iDim];
              
                // (p)
                for (CFuint iEq = 0; iEq < nbDepGradVar; ++iEq)
                {
                  const CFuint var = m_varToGradVarDep[m_pertVar][iEq];
                  
                  const CFreal divh_halfFaceJacobOther_lThis_dqduThis_dudu = m_gradVarsToStateJacobian[m_pertSide][pertSolIdx][m_pertVar][var] * divh_halfFaceJacobOther_lThis_dqduThis;
                    
                  m_tempFlux += m_riemannFluxGradJacobian[iFlxPnt][var][iDim] * divh_halfFaceJacobOther_lThis_dqduThis_dudu;
                }
                
//                const CFreal epsilon = 0.5*(m_epsilonLR[LEFT][iFlxPnt]+m_epsilonLR[RIGHT][iFlxPnt]);
//                
//                CFreal llavPart = epsilon * m_unitNormalFlxPnts[iFlxPnt][iDim];
//                llavPart *= divh_halfFaceJacobOther_lThis_dqduThis;
//                
//                // add part of analytical LLAV jacobian
//                m_tempFlux[m_pertVar] += llavPart;    
              }         
            
              acc.addValues(jSolIdxOther+otherSideTerm,pertSolIdx+pertSideTerm,m_pertVar,&m_tempFlux[0]);
            }
          }
          }
        }
        
        //// add the cross-element gradient part
        if (m_addFluxToGradCrossCellJacob)
        {
        // loop over the states to perturb the states (l)
        for (m_pertSol = 0; m_pertSol < m_nbrSolDep; ++m_pertSol)
        {
          const CFuint pertSolIdx = (*m_flxSolDep)[flxPntIdxThis][m_pertSol];
  
          for (CFuint kSolPnt = 0; kSolPnt < m_nbrSolDep; ++kSolPnt)
          {
            const CFuint kSolIdxOther = (*m_flxSolDep)[flxPntIdxOther][kSolPnt];
            
            // add the first and second part of the discontinuous gradient part of the jacobian
            for (CFuint iInfluencedFlx = 0; iInfluencedFlx < m_nbrFlxDep; ++iInfluencedFlx)
            {
              const CFuint iInfluencedFlxIdx = (*m_solFlxDep)[kSolIdxOther][iInfluencedFlx];

              const CFuint dimOther = (*m_flxPntFlxDim)[iInfluencedFlxIdx];
              
              const CFreal lOther = (*m_solPolyValsAtFlxPnts)[iInfluencedFlxIdx][kSolIdxOther];
                
              for (CFuint jSolPnt = 0; jSolPnt < m_nbrSolSolDep; ++jSolPnt)
              {   
                const CFuint jSolIdx = (*m_solSolDep)[kSolIdxOther][jSolPnt];
                  
                m_tempFlux = 0.0;
                  
                // get the divergence of the correction function on this side
                const CFreal divh_lOther = -m_corrFctDiv[jSolIdx][iInfluencedFlxIdx] * lOther; 
                
                /// llav jacob to state part //// actually should go over all kSol
                m_tempFlux += m_llavFluxJacobian[m_pertSide][pertSolIdx][m_pertVar][iOtherSide][jSolIdx][dimOther] * divh_lOther;
              
                for (CFuint jDim = 0; jDim < m_dim; ++jDim)
                {
                  const CFreal divh_l_dqduOther = divh_lOther * m_gradientStateJacobian[iOtherSide][kSolIdxOther][m_pertSide][pertSolIdx][jDim];
                  
                  for (CFuint iEq = 0; iEq < nbDepGradVar; ++iEq)
                  {
                    const CFuint var = m_varToGradVarDep[m_pertVar][iEq];
                      
                    const CFreal divh_l_dqduOther_dudu = divh_l_dqduOther * m_gradVarsToStateJacobian[m_pertSide][pertSolIdx][m_pertVar][var];
                      
                    m_tempFlux += m_gradientFluxJacobian[iOtherSide][kSolIdxOther][var][jDim][dimOther] * divh_l_dqduOther_dudu;
                  }
                  
//                  CFreal llavPart = divh_l_dqduOther * m_solEpsilons[iOtherSide][kSolIdxOther];
//                  
//                  llavPart *= m_neighbCellFluxProjVects[iOtherSide][dimOther][kSolIdxOther][jDim];
//                  
//                  // add part of analytical LLAV Jacobian
//                  m_tempFlux[m_pertVar] += llavPart;
                }
                
                acc.addValues(jSolIdx+otherSideTerm,pertSolIdx+pertSideTerm,m_pertVar,&m_tempFlux[0]);
              }
            }
          }
        }
        
        // loop over the states to perturb the states (l)
        for (m_pertSol = 0; m_pertSol < m_nbrSolDep; ++m_pertSol)
        {
          const CFuint pertSolIdx = (*m_flxSolDep)[flxPntIdxThis][m_pertSol];
             
          for (CFuint kSolPnt = 0; kSolPnt < m_nbrSolDep; ++kSolPnt)
          {
            const CFuint kSolIdxOther = (*m_flxSolDep)[flxPntIdxOther][kSolPnt];
              
            for (CFuint jSolPnt = 0; jSolPnt < m_nbrSolSolDep; ++jSolPnt)
            {
              const CFuint jSolIdx = (*m_solSolDep)[kSolIdxOther][jSolPnt];
                
              m_tempFlux = 0.0;
                
              for (CFuint iDim = 0; iDim < m_dim; ++iDim)
              {
                const CFreal lOther = (*m_solPolyDerivAtSolPnts)[jSolIdx][iDim][kSolIdxOther];
                
                /// llav jacob to state part //// actually should go over all jSol
                m_tempFlux += m_llavFluxJacobian[m_pertSide][pertSolIdx][m_pertVar][iOtherSide][jSolIdx][iDim] * lOther;
                  
                for (CFuint jDim = 0; jDim < m_dim; ++jDim)
                { 
                  const CFreal l_dqduOther = lOther * m_gradientStateJacobian[iOtherSide][kSolIdxOther][m_pertSide][pertSolIdx][jDim];
                  
                  for (CFuint iEq = 0; iEq < nbDepGradVar; ++iEq)
                  {
                    const CFuint var = m_varToGradVarDep[m_pertVar][iEq];
                    
                    const CFreal l_dqduOther_dudu = m_gradVarsToStateJacobian[m_pertSide][pertSolIdx][m_pertVar][var] * l_dqduOther;
                      
                    m_tempFlux += m_gradientFluxJacobian[iOtherSide][kSolIdxOther][var][jDim][iDim] * l_dqduOther_dudu;
                  }
                  
//                  CFreal llavPart = l_dqduOther * m_solEpsilons[iOtherSide][kSolIdxOther];
//                  
//                  llavPart *= m_neighbCellFluxProjVects[iOtherSide][iDim][kSolIdxOther][jDim];
//                  
//                  // add part of analytical LLAV Jacobian
//                  m_tempFlux[m_pertVar] += llavPart;
                }
              }
                
              acc.addValues(jSolIdx+otherSideTerm,pertSolIdx+pertSideTerm,m_pertVar,&m_tempFlux[0]);
            }
          }
        }
        }
      }
    }
  }
  
//  if (m_cells[LEFT]->getID() == 1) 
//  {
//      //CFLog(INFO, "ACC: " << acc.getValue(0,4,3,3) << "\n");
//      //acc.printToScreen();
//  }
//  if (m_cells[RIGHT]->getID() == 1) 
//  {
//      //CFLog(INFO, "ACC: " << acc.getValue(4,4,3,3) << "\n");
//      //acc.printToScreen();
//  }

  if (getMethodData().doComputeJacobian())
  {
    // add the values to the jacobian matrix
    m_lss->getMatrix()->addValues(acc);
  }

  // reset to zero the entries in the block accumulator
  acc.reset();
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffLLAVJacobFluxReconstructionTurb::computeOneJacobDiffFaceTerm(const CFuint side)
{
  // get residual factor
  const CFreal resFactor = getMethodData().getResFactor();

  // dereference accumulator
  BlockAccumulator& acc = *m_acc;
  
  CFuint solIdx = 0;
  for (m_pertSide = 0; m_pertSide < 2; ++m_pertSide)
  {
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol, ++solIdx)
    {
      acc.setRowColIndex(solIdx,(*m_states[m_pertSide])[iSol]->getLocalID());
    }
  }

  //// compute the needed flux jacobians
  
  initJacobianComputation();
  
  computeCellFluxJacobianNum(resFactor);
  
  computeRiemannFluxJacobianNum(resFactor);
  
  computeFluxToGradJacobianNum(resFactor);
  
  if (m_addRiemannToGradJacob || m_addRiemannToGradCrossCellJacob) computeRiemannFluxToGradJacobianNum(resFactor);
  
  computeGradToStateJacobianAna();
  
  computeGradVarsToStateJacobianNum();
  
  computeEpsToStateJacobianAna();
  
  computeLLAVCellFluxJacobianAna(resFactor);
  
  computeLLAVRiemannFluxJacobianAna(resFactor);
  
  //// add the total jacobians to the system jacobian
  
  // loop over left and right cell to add the discontinuous (cell) part to the jacobian
  for (m_pertSide = 0; m_pertSide < 2; ++m_pertSide)
  {
    // make sure this is only done once per cell
    if (!m_cellFlags[m_cells[m_pertSide]->getID()]) 
    {
      // term depending on iSide
      const CFuint pertSideTerm = m_pertSide*m_nbrSolPnts;

      // loop over the states to which to derive (l)
      for (m_pertSol = 0; m_pertSol < m_nbrSolPnts; ++m_pertSol)
      {
        // loop over the variables in the state (k)
        for (m_pertVar = 0; m_pertVar < m_nbrEqs; ++m_pertVar)
        {
          const CFuint nbDepGradVar = m_nbrVarToGradVarDep[m_pertVar];
            
          // add the discontinuous part of the jacobian related to the sol pnt (i)
          for (CFuint jSolPnt = 0; jSolPnt < m_nbrSolSolDep; ++jSolPnt)
          {
            const CFuint jSolIdx = (*m_solSolDep)[m_pertSol][jSolPnt];
            
            m_tempFlux = 0.;

            // (d)
            for (CFuint iDim = 0; iDim < m_dim; ++iDim)
            {
              const CFreal polyCoef = (*m_solPolyDerivAtSolPnts)[jSolIdx][iDim][m_pertSol]; 
          
              m_tempFlux += m_fluxJacobian[m_pertSide][m_pertSol][m_pertVar][iDim] * polyCoef;
            }
            
            acc.addValues(jSolIdx+pertSideTerm,m_pertSol+pertSideTerm,m_pertVar,&m_tempFlux[0]);
          }
          
          // add the discontinuous gradient part of the jacobian (m)
          for (CFuint kSolPnt = 0; kSolPnt < m_nbrSolSolDep; ++kSolPnt)
          {
            const CFuint kSolIdx = (*m_solSolDep)[m_pertSol][kSolPnt];
          
            // (i)
            for (CFuint jSol = 0; jSol < m_nbrSolSolDep; ++jSol)
            {
              const CFuint jSolIdx = (*m_solSolDep)[kSolIdx][jSol];
              
              m_tempFlux = 0.0;
                
              // (d)
              for (CFuint iDim = 0; iDim < m_dim; ++iDim)
              {
                const CFreal dl = (*m_solPolyDerivAtSolPnts)[jSolIdx][iDim][kSolIdx];
                
                /// llav jacob to state part //// should actually be added for all jSol
                m_tempFlux += m_llavFluxJacobian[m_pertSide][m_pertSol][m_pertVar][m_pertSide][jSolIdx][iDim] * dl;
                
                // (b)
                for (CFuint jDim = 0; jDim < m_dim; ++jDim)
                {
                  const CFreal dl_dqdu = dl * m_gradientStateJacobian[m_pertSide][kSolIdx][m_pertSide][m_pertSol][jDim];
                    
                  // (p)
                  for (CFuint iEq = 0; iEq < nbDepGradVar; ++iEq)
                  {
                    const CFuint var = m_varToGradVarDep[m_pertVar][iEq];
                    
                    const CFreal dl_dqdu_dudu = m_gradVarsToStateJacobian[m_pertSide][m_pertSol][m_pertVar][var] * dl_dqdu;
                      
                    m_tempFlux += m_gradientFluxJacobian[m_pertSide][kSolIdx][var][jDim][iDim] * dl_dqdu_dudu;
                  }
                  
//                  CFreal llavPart = m_solEpsilons[m_pertSide][kSolIdx] * m_neighbCellFluxProjVects[m_pertSide][iDim][kSolIdx][jDim];
//                  llavPart *= dl_dqdu;
//                  //if(m_cells[0]->getID()==1) CFLog(INFO, "pertSol: " << m_pertSol << ", pertVar: " << m_pertVar << ", llavJC: " << llavPart << "\n");
//                  // add part of analytical LLAV jacobian
//                  m_tempFlux[m_pertVar] += llavPart;
                }
              }
            
              acc.addValues(jSolIdx+pertSideTerm,m_pertSol+pertSideTerm,m_pertVar,&m_tempFlux[0]);
            }
          }
          
          // add the discontinuous part of the jacobian related to the flx pnt (f)
          for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFlxDep; ++iFlxPnt)
          {
            const CFuint flxIdx = (*m_solFlxDep)[m_pertSol][iFlxPnt];
            
            // (df)
            const CFuint dim = (*m_flxPntFlxDim)[flxIdx];
            
            m_temp = m_fluxJacobian[m_pertSide][m_pertSol][m_pertVar][dim] * (*m_solPolyValsAtFlxPnts)[flxIdx][m_pertSol];
            
            // add the second part of the discontinuous part of the jacobian (i)
            for (CFuint jSolPnt = 0; jSolPnt < m_nbrSolDep; ++jSolPnt)
            {
              const CFuint jSolIdx = (*m_flxSolDep)[flxIdx][jSolPnt];

              // get the divergence of the correction function
              const CFreal divh = m_corrFctDiv[jSolIdx][flxIdx];
                           
              m_tempFlux = -m_temp * divh;
            
              acc.addValues(jSolIdx+pertSideTerm,m_pertSol+pertSideTerm,m_pertVar,&m_tempFlux[0]);
            }
          }
          
          // add the second part of the discontinuous gradient part of the jacobian (m)
          for (CFuint kSolPnt = 0; kSolPnt < m_nbrSolSolDep; ++kSolPnt)
          {
            const CFuint kSolIdx = (*m_solSolDep)[m_pertSol][kSolPnt];
              
            for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFlxDep; ++iFlxPnt)
            {   
              const CFuint flxIdx = (*m_solFlxDep)[kSolIdx][iFlxPnt];
              
              // (df)
              const CFuint dim = (*m_flxPntFlxDim)[flxIdx];
              
              const CFreal l = (*m_solPolyValsAtFlxPnts)[flxIdx][kSolIdx];            
              
              // (i)
              for (CFuint jSolPnt = 0; jSolPnt < m_nbrSolDep; ++jSolPnt)
              {
                const CFuint jSolIdx = (*m_flxSolDep)[flxIdx][jSolPnt];
              
                m_tempFlux = 0.0;

                const CFreal divh_l = -m_corrFctDiv[jSolIdx][flxIdx] * l;
                
                /// llav jacob to state part //// actually should loop over all kSol
                m_tempFlux += m_llavFluxJacobian[m_pertSide][m_pertSol][m_pertVar][m_pertSide][jSolIdx][dim] * divh_l;
              
                // (b)
                for (CFuint jDim = 0; jDim < m_dim; ++jDim)
                {
                  const CFreal divh_l_dqdu = divh_l * m_gradientStateJacobian[m_pertSide][kSolIdx][m_pertSide][m_pertSol][jDim];
                  
                  // (p)
                  for (CFuint iEq = 0; iEq < nbDepGradVar; ++iEq)
                  {
                    const CFuint var = m_varToGradVarDep[m_pertVar][iEq];
                    
                    const CFreal divh_l_dqdu_dudu = divh_l_dqdu * m_gradVarsToStateJacobian[m_pertSide][m_pertSol][m_pertVar][var];             
                      
                    m_tempFlux += divh_l_dqdu_dudu * m_gradientFluxJacobian[m_pertSide][kSolIdx][var][jDim][dim];
                  }
                  
//                  CFreal llavPart = divh_l_dqdu * m_solEpsilons[m_pertSide][kSolIdx];
//                  
//                  llavPart *= m_neighbCellFluxProjVects[m_pertSide][dim][kSolIdx][jDim];
//                  //if(m_cells[0]->getID()==1) CFLog(INFO, "pertSol: " << m_pertSol << ", pertVar: " << m_pertVar << ", llavJF: " << llavPart << "\n");
//                  // add part of analytical LLAV Jacobian
//                  m_tempFlux[m_pertVar] += llavPart;
                }
                
                acc.addValues(jSolIdx+pertSideTerm,m_pertSol+pertSideTerm,m_pertVar,&m_tempFlux[0]);
              }
            }
          }
        }
      }
    }
  }
  
  // loop over left and right cell to add the riemann flux (face) part to the jacobian
  for (m_pertSide = 0; m_pertSide < 2; ++m_pertSide)
  {
    // variable for the other side
    const CFuint iOtherSide = m_pertSide == LEFT ? RIGHT : LEFT;
    
    // term depending on iSide
    const CFuint pertSideTerm = m_pertSide*m_nbrSolPnts;

    // term depending on iOtherSide
    const CFuint otherSideTerm = iOtherSide*m_nbrSolPnts;
    
    // loop over the variables in the state (k)
    for (m_pertVar = 0; m_pertVar < m_nbrEqs; ++m_pertVar)
    { 
      const CFuint nbDepGradVar = m_nbrVarToGradVarDep[m_pertVar];
        
      // loop over face flx pnts (f)
      for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
      {
        // local flux point indices in the left and right cell
        const CFuint flxPntIdxThis = (*m_faceFlxPntConnPerOrient)[m_orient][m_pertSide][iFlxPnt];
        const CFuint flxPntIdxOther = (*m_faceFlxPntConnPerOrient)[m_orient][iOtherSide][iFlxPnt];
        
        m_temp = m_riemannFluxJacobian[m_pertSide][iFlxPnt][m_pertVar]*m_faceJacobVecSizeFlxPnts[iFlxPnt][m_pertSide];
        m_tempOther = m_riemannFluxJacobian[m_pertSide][iFlxPnt][m_pertVar]*m_faceJacobVecSizeFlxPnts[iFlxPnt][iOtherSide];
        
        const CFreal halfFaceJacob = 0.5 * m_faceJacobVecSizeFlxPnts[iFlxPnt][m_pertSide];
        
        // loop over the states to perturb the states (l)
        for (m_pertSol = 0; m_pertSol < m_nbrSolDep; ++m_pertSol)
        {
          const CFuint pertSolIdx = (*m_flxSolDep)[flxPntIdxThis][m_pertSol];
          
          m_temp2 = m_temp * (*m_solPolyValsAtFlxPnts)[flxPntIdxThis][pertSolIdx];
          m_tempOther2 = m_tempOther * (*m_solPolyValsAtFlxPnts)[flxPntIdxThis][pertSolIdx];
            
          // add the second part of the discontinuous part of the jacobian (i)
          for (CFuint jSolPnt = 0; jSolPnt < m_nbrSolDep; ++jSolPnt)
          {
            const CFuint jSolIdxThis = (*m_flxSolDep)[flxPntIdxThis][jSolPnt];
            const CFuint jSolIdxOther = (*m_flxSolDep)[flxPntIdxOther][jSolPnt];

            // get the divergence of the correction function on this side
            CFreal divh = m_corrFctDiv[jSolIdxThis][flxPntIdxThis];
                          
            // add part on this side of face
            m_tempFlux = m_temp2 * divh;
              
            acc.addValues(jSolIdxThis+pertSideTerm,pertSolIdx+pertSideTerm,m_pertVar,&m_tempFlux[0]);
          
            // get the divergence of the correction function on other side
            divh = m_corrFctDiv[jSolIdxOther][flxPntIdxOther];
                        
            // add cross-cell part 
            m_tempFlux = m_tempOther2 * divh;   
              
            acc.addValues(jSolIdxOther+otherSideTerm,pertSolIdx+pertSideTerm,m_pertVar,&m_tempFlux[0]);
          }
        }
        
        // loop over the states to perturb the states (l) for LLAV to state part
        for (m_pertSol = 0; m_pertSol < m_nbrSolPnts; ++m_pertSol)
        { 
          m_temp2 = m_llavRiemannFluxJacobian[m_pertSide][m_pertSol][m_pertVar][iFlxPnt] * m_faceJacobVecSizeFlxPnts[iFlxPnt][m_pertSide];
          m_tempOther2 = m_llavRiemannFluxJacobian[m_pertSide][m_pertSol][m_pertVar][iFlxPnt] * m_faceJacobVecSizeFlxPnts[iFlxPnt][iOtherSide];
          
          // add the LLAV interface part of the jacobian (i)
          for (CFuint jSolPnt = 0; jSolPnt < m_nbrSolDep; ++jSolPnt)
          {
            const CFuint jSolIdxThis = (*m_flxSolDep)[flxPntIdxThis][jSolPnt];
            const CFuint jSolIdxOther = (*m_flxSolDep)[flxPntIdxOther][jSolPnt];

            // get the divergence of the correction function on this side
            CFreal divh = m_corrFctDiv[jSolIdxThis][flxPntIdxThis];
                          
            // add part on this side of face
            m_tempFlux = m_temp2 * divh;
              
            acc.addValues(jSolIdxThis+pertSideTerm,m_pertSol+pertSideTerm,m_pertVar,&m_tempFlux[0]);
          
            // get the divergence of the correction function on other side
            divh = m_corrFctDiv[jSolIdxOther][flxPntIdxOther];
                        
            // add cross-cell part 
            m_tempFlux = m_tempOther2 * divh;   
              
            acc.addValues(jSolIdxOther+otherSideTerm,m_pertSol+pertSideTerm,m_pertVar,&m_tempFlux[0]);
          }
        }
        
        for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
        { 
          m_needToAddSolPnt[iSol] = true;
        }
            
        // (i)
        for (CFuint jSolPnt = 0; jSolPnt < m_nbrSolDep; ++jSolPnt)
        {
          const CFuint jSolIdxThis = (*m_flxSolDep)[flxPntIdxThis][jSolPnt];
          const CFuint jSolIdxOther = (*m_flxSolDep)[flxPntIdxOther][jSolPnt];

          // get the divergence of the correction function on this side
          CFreal divh = m_corrFctDiv[jSolIdxThis][flxPntIdxThis];
            
          const CFreal divh_halfFaceJacob = divh * halfFaceJacob;
          
          // loop over the states to perturb the states (l)
          for (CFuint m_pertSol = 0; m_pertSol < m_nbrSolDep; ++m_pertSol)
          {   
            const CFuint pertSolIdx = (*m_flxSolDep)[flxPntIdxThis][m_pertSol];
            
            m_needToAddSolPnt[pertSolIdx] = false;
            
            if (m_addRiemannToGradJacob)
            {
            // add part on this side of face
            m_tempFlux = 0.0;

            // (m)
            for (CFuint kSolPnt = 0; kSolPnt < m_nbrSolDep; ++kSolPnt)
            {
              const CFuint kSolIdx = (*m_flxSolDep)[flxPntIdxThis][kSolPnt];
              const CFuint kSolIdxOther = (*m_flxSolDep)[flxPntIdxOther][kSolPnt];
              
              const CFreal divh_halfFaceJacob_l = divh_halfFaceJacob * (*m_solPolyValsAtFlxPnts)[flxPntIdxThis][kSolIdx];
              const CFreal divh_halfFaceJacob_lOther = divh_halfFaceJacob * (*m_solPolyValsAtFlxPnts)[flxPntIdxOther][kSolIdxOther];
              
              /// llav jacob to state part
              //m_tempFlux += m_llavRiemannFluxJacobian[m_pertSide][pertSolIdx][m_pertVar][iFlxPnt] * divh_halfFaceJacob_l;
              //m_tempFlux += m_llavRiemannFluxJacobian[m_pertSide][pertSolIdx][m_pertVar][iFlxPnt] * divh_halfFaceJacob_lOther;
                
              // (b)
              for (CFuint iDim = 0; iDim < m_dim; ++iDim)
              {
                const CFreal divh_halfFaceJacob_l_dqdu = divh_halfFaceJacob_l * m_gradientStateJacobian[m_pertSide][kSolIdx][m_pertSide][pertSolIdx][iDim];
                const CFreal divh_halfFaceJacob_l_dqduOther = divh_halfFaceJacob_lOther * m_gradientStateJacobian[iOtherSide][kSolIdxOther][m_pertSide][pertSolIdx][iDim];
                
                // (p)
                for (CFuint iEq = 0; iEq < nbDepGradVar; ++iEq)
                {
                  const CFuint var = m_varToGradVarDep[m_pertVar][iEq];
                  
                  m_temp =  m_riemannFluxGradJacobian[iFlxPnt][var][iDim] * m_gradVarsToStateJacobian[m_pertSide][pertSolIdx][m_pertVar][var];
                    
                  m_tempFlux += m_temp * divh_halfFaceJacob_l_dqdu;
                    
                  m_tempFlux += m_temp * divh_halfFaceJacob_l_dqduOther;
                }
                
//                const CFreal epsilon = 0.5*(m_epsilonLR[LEFT][iFlxPnt]+m_epsilonLR[RIGHT][iFlxPnt]);
//                
//                const CFreal llavPart = epsilon * m_unitNormalFlxPnts[iFlxPnt][iDim];
//                
//                // add part of analytical LLAV jacobian
//                m_tempFlux[m_pertVar] += llavPart * divh_halfFaceJacob_l_dqdu;    
//                m_tempFlux[m_pertVar] += llavPart * divh_halfFaceJacob_l_dqduOther;  
              }         
            }
              
            acc.addValues(jSolIdxThis+pertSideTerm,pertSolIdx+pertSideTerm,m_pertVar,&m_tempFlux[0]);
          
            // get the divergence of the correction function on other side
            divh = m_corrFctDiv[jSolIdxOther][flxPntIdxOther];
            
            const CFreal divh_halfFaceJacobOther = 0.5 * divh * m_faceJacobVecSizeFlxPnts[iFlxPnt][iOtherSide];

            // add cross-cell part 
            m_tempFlux = 0.0;

            // (m)
            for (CFuint kSolPnt = 0; kSolPnt < m_nbrSolDep; ++kSolPnt)
            {
              const CFuint kSolIdx = (*m_flxSolDep)[flxPntIdxThis][kSolPnt];
              const CFuint kSolIdxOther = (*m_flxSolDep)[flxPntIdxOther][kSolPnt];
              
              const CFreal divh_halfFaceJacobOther_lThis = divh_halfFaceJacobOther * (*m_solPolyValsAtFlxPnts)[flxPntIdxThis][kSolIdx];
              const CFreal divh_halfFaceJacobOther_lOther = divh_halfFaceJacobOther * (*m_solPolyValsAtFlxPnts)[flxPntIdxOther][kSolIdxOther];
              
              /// llav jacob to state part
              //m_tempFlux += m_llavRiemannFluxJacobian[m_pertSide][pertSolIdx][m_pertVar][iFlxPnt] * divh_halfFaceJacobOther_lThis;
              //m_tempFlux += m_llavRiemannFluxJacobian[m_pertSide][pertSolIdx][m_pertVar][iFlxPnt] * divh_halfFaceJacobOther_lOther;
                
              // (b)
              for (CFuint iDim = 0; iDim < m_dim; ++iDim)
              {
                const CFreal divh_halfFaceJacobOther_lThis_dqduThis = divh_halfFaceJacobOther_lThis * m_gradientStateJacobian[m_pertSide][kSolIdx][m_pertSide][pertSolIdx][iDim];
                const CFreal divh_halfFaceJacobOther_lOther_dqduOther = divh_halfFaceJacobOther_lOther * m_gradientStateJacobian[iOtherSide][kSolIdxOther][m_pertSide][pertSolIdx][iDim];
              
                // (p)
                for (CFuint iEq = 0; iEq < nbDepGradVar; ++iEq)
                {
                  const CFuint var = m_varToGradVarDep[m_pertVar][iEq];
                  
                  m_temp = m_riemannFluxGradJacobian[iFlxPnt][var][iDim] * m_gradVarsToStateJacobian[m_pertSide][pertSolIdx][m_pertVar][var];
                    
                  m_tempFlux += m_temp * divh_halfFaceJacobOther_lThis_dqduThis;
                    
                  m_tempFlux += m_temp * divh_halfFaceJacobOther_lOther_dqduOther;
                }
                
//                const CFreal epsilon = 0.5*(m_epsilonLR[LEFT][iFlxPnt]+m_epsilonLR[RIGHT][iFlxPnt]);
//                
//                const CFreal llavPart = epsilon * m_unitNormalFlxPnts[iFlxPnt][iDim];
//                
//                // add part of analytical LLAV jacobian
//                m_tempFlux[m_pertVar] += llavPart * divh_halfFaceJacobOther_lThis_dqduThis;    
//                m_tempFlux[m_pertVar] += llavPart * divh_halfFaceJacobOther_lOther_dqduOther;
              }         
            }
            
            acc.addValues(jSolIdxOther+otherSideTerm,pertSolIdx+pertSideTerm,m_pertVar,&m_tempFlux[0]);
          }
          }
          
          if (m_addRiemannToGradCrossCellJacob)
          {
          // loop over the states to perturb the states (l)
          for (CFuint pertSolIdx = 0; pertSolIdx < m_nbrSolPnts; ++pertSolIdx)
          {
            CFuint dependingKSol = 1000;
              
            if (m_needToAddSolPnt[pertSolIdx])
            {
              // add part on this side of face
              m_tempFlux = 0.0;

              // (m)
              for (CFuint kSolPnt = 0; kSolPnt < m_nbrSolDep; ++kSolPnt)
              {
                const CFuint kSolIdx = (*m_flxSolDep)[flxPntIdxThis][kSolPnt];
              
                for (CFuint lSol = 0; lSol < m_nbrSolSolDep; ++lSol)
                {
                  const CFuint lSolIdx = (*m_solSolDep)[pertSolIdx][lSol]; 
                
                  if (lSolIdx == kSolIdx)
                  {
                    dependingKSol = kSolIdx;
                    break;
                  }
                }
              }
  
              const CFreal divh_halfFaceJacob_l = divh_halfFaceJacob * (*m_solPolyValsAtFlxPnts)[flxPntIdxThis][dependingKSol];
              
              /// llav jacob to state part
              //m_tempFlux += m_llavRiemannFluxJacobian[m_pertSide][pertSolIdx][m_pertVar][iFlxPnt] * divh_halfFaceJacob_l;
                
              // (b)
              for (CFuint iDim = 0; iDim < m_dim; ++iDim)
              {
                const CFreal divh_halfFaceJacob_l_dqdu = divh_halfFaceJacob_l * m_gradientStateJacobian[m_pertSide][dependingKSol][m_pertSide][pertSolIdx][iDim];
                
                // (p)
                for (CFuint iEq = 0; iEq < nbDepGradVar; ++iEq)
                {
                  const CFuint var = m_varToGradVarDep[m_pertVar][iEq];
                  
                  const CFreal divh_halfFaceJacob_l_dqdu_dudu = m_gradVarsToStateJacobian[m_pertSide][pertSolIdx][m_pertVar][var] * divh_halfFaceJacob_l_dqdu;
                    
                  m_tempFlux += m_riemannFluxGradJacobian[iFlxPnt][var][iDim] * divh_halfFaceJacob_l_dqdu_dudu;
                }
                
//                const CFreal epsilon = 0.5*(m_epsilonLR[LEFT][iFlxPnt]+m_epsilonLR[RIGHT][iFlxPnt]);
//                
//                CFreal llavPart =  epsilon * m_unitNormalFlxPnts[iFlxPnt][iDim];
//                
//                llavPart *= divh_halfFaceJacob_l_dqdu;
//                
//                // add part of analytical LLAV jacobian
//                m_tempFlux[m_pertVar] += llavPart;     
              } 
              
              acc.addValues(jSolIdxThis+pertSideTerm,pertSolIdx+pertSideTerm,m_pertVar,&m_tempFlux[0]);
              
          
              // get the divergence of the correction function on other side
              divh = m_corrFctDiv[jSolIdxOther][flxPntIdxOther];
            
              const CFreal divh_halfFaceJacobOther = 0.5 * divh * m_faceJacobVecSizeFlxPnts[iFlxPnt][iOtherSide];

              // add cross-cell part 
              m_tempFlux = 0.0;
              
              const CFreal divh_halfFaceJacobOther_lThis = divh_halfFaceJacobOther * (*m_solPolyValsAtFlxPnts)[flxPntIdxThis][dependingKSol];
                
              /// llav jacob to state part
              //m_tempFlux += m_llavRiemannFluxJacobian[m_pertSide][pertSolIdx][m_pertVar][iFlxPnt] * divh_halfFaceJacobOther_lThis;
              
              // (b)
              for (CFuint iDim = 0; iDim < m_dim; ++iDim)
              {
                const CFreal divh_halfFaceJacobOther_lThis_dqduThis = divh_halfFaceJacobOther_lThis * m_gradientStateJacobian[m_pertSide][dependingKSol][m_pertSide][pertSolIdx][iDim];
              
                // (p)
                for (CFuint iEq = 0; iEq < nbDepGradVar; ++iEq)
                {
                  const CFuint var = m_varToGradVarDep[m_pertVar][iEq];
                  
                  const CFreal divh_halfFaceJacobOther_lThis_dqduThis_dudu = m_gradVarsToStateJacobian[m_pertSide][pertSolIdx][m_pertVar][var] * divh_halfFaceJacobOther_lThis_dqduThis;
                    
                  m_tempFlux += m_riemannFluxGradJacobian[iFlxPnt][var][iDim] * divh_halfFaceJacobOther_lThis_dqduThis_dudu;
                }
                
//                const CFreal epsilon = 0.5*(m_epsilonLR[LEFT][iFlxPnt]+m_epsilonLR[RIGHT][iFlxPnt]);
//                
//                CFreal llavPart = epsilon * m_unitNormalFlxPnts[iFlxPnt][iDim];
//                llavPart *= divh_halfFaceJacobOther_lThis_dqduThis;
//                
//                // add part of analytical LLAV jacobian
//                m_tempFlux[m_pertVar] += llavPart;    
              }         
            
              acc.addValues(jSolIdxOther+otherSideTerm,pertSolIdx+pertSideTerm,m_pertVar,&m_tempFlux[0]);
            }
          }
          }
        }
        
        //// add the cross-element gradient part
        if (m_addFluxToGradCrossCellJacob)
        {
        // loop over the states to perturb the states (l)
        for (m_pertSol = 0; m_pertSol < m_nbrSolDep; ++m_pertSol)
        {
          const CFuint pertSolIdx = (*m_flxSolDep)[flxPntIdxThis][m_pertSol];
  
          for (CFuint kSolPnt = 0; kSolPnt < m_nbrSolDep; ++kSolPnt)
          {
            const CFuint kSolIdxOther = (*m_flxSolDep)[flxPntIdxOther][kSolPnt];
            
            // add the first and second part of the discontinuous gradient part of the jacobian
            for (CFuint iInfluencedFlx = 0; iInfluencedFlx < m_nbrFlxDep; ++iInfluencedFlx)
            {
              const CFuint iInfluencedFlxIdx = (*m_solFlxDep)[kSolIdxOther][iInfluencedFlx];

              const CFuint dimOther = (*m_flxPntFlxDim)[iInfluencedFlxIdx];
              
              const CFreal lOther = (*m_solPolyValsAtFlxPnts)[iInfluencedFlxIdx][kSolIdxOther];
                
              for (CFuint jSolPnt = 0; jSolPnt < m_nbrSolSolDep; ++jSolPnt)
              {   
                const CFuint jSolIdx = (*m_solSolDep)[kSolIdxOther][jSolPnt];
                  
                m_tempFlux = 0.0;
                  
                // get the divergence of the correction function on this side
                const CFreal divh_lOther = -m_corrFctDiv[jSolIdx][iInfluencedFlxIdx] * lOther; 
                
                /// llav jacob to state part //// actually should go over all kSol
                m_tempFlux += m_llavFluxJacobian[m_pertSide][pertSolIdx][m_pertVar][iOtherSide][jSolIdx][dimOther] * divh_lOther;
              
                for (CFuint jDim = 0; jDim < m_dim; ++jDim)
                {
                  const CFreal divh_l_dqduOther = divh_lOther * m_gradientStateJacobian[iOtherSide][kSolIdxOther][m_pertSide][pertSolIdx][jDim];
                  
                  for (CFuint iEq = 0; iEq < nbDepGradVar; ++iEq)
                  {
                    const CFuint var = m_varToGradVarDep[m_pertVar][iEq];
                      
                    const CFreal divh_l_dqduOther_dudu = divh_l_dqduOther * m_gradVarsToStateJacobian[m_pertSide][pertSolIdx][m_pertVar][var];
                      
                    m_tempFlux += m_gradientFluxJacobian[iOtherSide][kSolIdxOther][var][jDim][dimOther] * divh_l_dqduOther_dudu;
                  }
                  
//                  CFreal llavPart = divh_l_dqduOther * m_solEpsilons[iOtherSide][kSolIdxOther];
//                  
//                  llavPart *= m_neighbCellFluxProjVects[iOtherSide][dimOther][kSolIdxOther][jDim];
//                  
//                  // add part of analytical LLAV Jacobian
//                  m_tempFlux[m_pertVar] += llavPart;
                }
                
                acc.addValues(jSolIdx+otherSideTerm,pertSolIdx+pertSideTerm,m_pertVar,&m_tempFlux[0]);
              }
            }
          }
        }
        
        // loop over the states to perturb the states (l)
        for (m_pertSol = 0; m_pertSol < m_nbrSolDep; ++m_pertSol)
        {
          const CFuint pertSolIdx = (*m_flxSolDep)[flxPntIdxThis][m_pertSol];
             
          for (CFuint kSolPnt = 0; kSolPnt < m_nbrSolDep; ++kSolPnt)
          {
            const CFuint kSolIdxOther = (*m_flxSolDep)[flxPntIdxOther][kSolPnt];
              
            for (CFuint jSolPnt = 0; jSolPnt < m_nbrSolSolDep; ++jSolPnt)
            {
              const CFuint jSolIdx = (*m_solSolDep)[kSolIdxOther][jSolPnt];
                
              m_tempFlux = 0.0;
                
              for (CFuint iDim = 0; iDim < m_dim; ++iDim)
              {
                const CFreal lOther = (*m_solPolyDerivAtSolPnts)[jSolIdx][iDim][kSolIdxOther];
                
                /// llav jacob to state part //// actually should go over all jSol
                m_tempFlux += m_llavFluxJacobian[m_pertSide][pertSolIdx][m_pertVar][iOtherSide][jSolIdx][iDim] * lOther;
                  
                for (CFuint jDim = 0; jDim < m_dim; ++jDim)
                { 
                  const CFreal l_dqduOther = lOther * m_gradientStateJacobian[iOtherSide][kSolIdxOther][m_pertSide][pertSolIdx][jDim];
                  
                  for (CFuint iEq = 0; iEq < nbDepGradVar; ++iEq)
                  {
                    const CFuint var = m_varToGradVarDep[m_pertVar][iEq];
                    
                    const CFreal l_dqduOther_dudu = m_gradVarsToStateJacobian[m_pertSide][pertSolIdx][m_pertVar][var] * l_dqduOther;
                      
                    m_tempFlux += m_gradientFluxJacobian[iOtherSide][kSolIdxOther][var][jDim][iDim] * l_dqduOther_dudu;
                  }
                  
//                  CFreal llavPart = l_dqduOther * m_solEpsilons[iOtherSide][kSolIdxOther];
//                  
//                  llavPart *= m_neighbCellFluxProjVects[iOtherSide][iDim][kSolIdxOther][jDim];
//                  
//                  // add part of analytical LLAV Jacobian
//                  m_tempFlux[m_pertVar] += llavPart;
                }
              }
                
              acc.addValues(jSolIdx+otherSideTerm,pertSolIdx+pertSideTerm,m_pertVar,&m_tempFlux[0]);
            }
          }
        }
        }
      }
    }
  }
  
//  if (m_cells[LEFT]->getID() == 1) 
//  {
//      //CFLog(INFO, "ACC: " << acc.getValue(0,4,3,3) << "\n");
//      //acc.printToScreen();
//  }
//  if (m_cells[RIGHT]->getID() == 1) 
//  {
//      //CFLog(INFO, "ACC: " << acc.getValue(4,4,3,3) << "\n");
//      //acc.printToScreen();
//  }

  if (getMethodData().doComputeJacobian())
  {
    // add the values to the jacobian matrix
    m_lss->getMatrix()->addValues(acc);
  }

  // reset to zero the entries in the block accumulator
  acc.reset();
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffLLAVJacobFluxReconstructionTurb::setup()
{
  CFAUTOTRACE;
  
  // setup parent class
  ConvDiffLLAVJacobFluxReconstructionNS::setup();
  
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
