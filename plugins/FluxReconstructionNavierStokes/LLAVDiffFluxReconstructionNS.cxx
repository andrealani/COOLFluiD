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

#include "NavierStokes/EulerTerm.hh"

#include "NavierStokes/NavierStokesVarSet.hh"
#include "NavierStokes/EulerVarSet.hh"

#include "FluxReconstructionNavierStokes/LLAVDiffFluxReconstructionNS.hh"
#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"

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

MethodCommandProvider< LLAVDiffFluxReconstructionNS,
		       FluxReconstructionSolverData,
		       FluxReconstructionModule >
LLAVDiffFluxReconstructionNSFluxReconstructionProvider("LLAVDiffNS");

//////////////////////////////////////////////////////////////////////////////
  
LLAVDiffFluxReconstructionNS::LLAVDiffFluxReconstructionNS(const std::string& name) :
  LLAVDiffFluxReconstruction(name),
  m_gradsBackUp(),
  m_eulerVarSet(CFNULL),
  m_msEulerTerm(CFNULL),
  m_nbrSpecies(),
  m_pData(),
  m_pData2(),
  m_eulerVarSet2(CFNULL),
  m_tempGradTermL(),
  m_tempGradTermR(),
  m_diffusiveVarSet(CFNULL),
  m_tempStatesL(),
  m_tempStatesR(),
  m_cellGradsAV(),
  m_avgGradAV(),
  m_cellGradFlxPntAV(),
  m_dampCoeff()
  {
  }

//////////////////////////////////////////////////////////////////////////////

void LLAVDiffFluxReconstructionNS::configure ( Config::ConfigArgs& args )
{
  LLAVDiffFluxReconstruction::configure(args);
}  

//////////////////////////////////////////////////////////////////////////////

void LLAVDiffFluxReconstructionNS::prepareFluxComputation()
{
  const bool isPerturb = this->getMethodData().isPerturb();
  const CFuint iPerturbVar = this->getMethodData().iPerturbVar();

  m_diffusiveVarSet->setComposition(m_avgSol, isPerturb, iPerturbVar);
}

//////////////////////////////////////////////////////////////////////////////

void LLAVDiffFluxReconstructionNS::computeWaveSpeedUpdates(vector< CFreal >& waveSpeedUpd)
{
  // compute the wave speed updates for the neighbouring cells
  cf_assert(waveSpeedUpd.size() == 2);
  CFreal visc = 1.0;

  const CFreal dynVisc = m_diffusiveVarSet->getCurrDynViscosity();
  const CFreal factorPr = min(m_diffusiveVarSet->getModel().getPrandtl(),1.0); 
  cf_assert(factorPr>0.0);
  
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    waveSpeedUpd[iSide] = 0.0;
    for (CFuint iFlx = 0; iFlx < m_cellFlx[iSide].size(); ++iFlx)
    {
      const CFreal jacobXJacobXIntCoef = m_faceJacobVecAbsSizeFlxPnts[iFlx]*
                                 m_faceJacobVecAbsSizeFlxPnts[iFlx]*
                                   (*m_faceIntegrationCoefs)[iFlx]*
                                   m_cflConvDiffRatio;
      const CFreal rho = m_diffusiveVarSet->getDensity(*(m_cellStatesFlxPnt[iSide][iFlx]));
      visc = dynVisc/rho/factorPr;
      
      if (m_addUpdCoeff)
      {
        const CFreal epsilon = 0.5*(m_epsilonLR[LEFT][iFlx]+m_epsilonLR[RIGHT][iFlx]);
        const CFreal viscCoef = computeViscCoef(m_cellStatesFlxPnt[iSide][iFlx]);
        visc += epsilon*viscCoef;
      }
      
      // transform update states to physical data to calculate eigenvalues
      waveSpeedUpd[iSide] += visc*jacobXJacobXIntCoef/m_cellVolume[iSide];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void LLAVDiffFluxReconstructionNS::setFaceData(CFuint faceID)
{
  LLAVDiffFluxReconstruction::setFaceData(faceID);

  // get the gradients datahandle
  DataHandle< vector< RealVector > > gradientsAV = socket_gradientsAV.getDataHandle();

  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
    {
      const CFuint stateID = (*(m_states[iSide]))[iState]->getLocalID();
      m_cellGradsAV[iSide][iState] = &gradientsAV[stateID];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void LLAVDiffFluxReconstructionNS::setCellData()
{
  LLAVDiffFluxReconstruction::setCellData();

  // get the gradients datahandle
  DataHandle< vector< RealVector > > gradientsAV = socket_gradientsAV.getDataHandle();

  // get the grads in the current cell
  for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
  {
    const CFuint stateID = (*(m_cellStates))[iState]->getLocalID();
    m_cellGradsAV[0][iState] = &gradientsAV[stateID];
  }
}

//////////////////////////////////////////////////////////////////////////////

void LLAVDiffFluxReconstructionNS::computeInterfaceFlxCorrection()
{   
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_tempStatesL[iFlx] = m_cellStatesFlxPnt[LEFT][iFlx]->getData();
    m_tempStatesR[iFlx] = m_cellStatesFlxPnt[RIGHT][iFlx]->getData();
  }
  
  m_diffusiveVarSet->setGradientVars(m_tempStatesL,m_tempGradTermL,m_nbrFaceFlxPnts);
  m_diffusiveVarSet->setGradientVars(m_tempStatesR,m_tempGradTermR,m_nbrFaceFlxPnts);
  
  // Loop over the flux points to calculate FI
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  { 
    // get AV
    const CFreal epsilon = 0.5*(m_epsilonLR[LEFT][iFlxPnt]+m_epsilonLR[RIGHT][iFlxPnt]);
      
    // compute the average sol and grad to use the BR2 scheme
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      *(m_avgGrad[iVar]) = (*(m_cellGradFlxPnt[LEFT][iFlxPnt][iVar]) + *(m_cellGradFlxPnt[RIGHT][iFlxPnt][iVar]))/2.0;
      
      *(m_avgGradAV[iVar]) = (*(m_cellGradFlxPntAV[LEFT][iFlxPnt][iVar]) + *(m_cellGradFlxPntAV[RIGHT][iFlxPnt][iVar]))/2.0;
      
      m_avgSol[iVar] = ((*(m_cellStatesFlxPnt[LEFT][iFlxPnt]))[iVar] + (*(m_cellStatesFlxPnt[RIGHT][iFlxPnt]))[iVar])/2.0;
    }

    // damping factor
    const CFreal dampFactor = m_dampCoeff*m_faceInvCharLengths[iFlxPnt];

    // compute averaged (damped) gradients
    for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
    {
      // compute damping term
      const RealVector dGradVarXNormalAV = ((*m_cellStatesFlxPnt[LEFT][iFlxPnt])[iGrad] - (*m_cellStatesFlxPnt[RIGHT][iFlxPnt])[iGrad])*m_unitNormalFlxPnts[iFlxPnt];
      
      const RealVector dGradVarXNormal = (m_tempGradTermL(iGrad,iFlxPnt) - m_tempGradTermR(iGrad,iFlxPnt))*m_unitNormalFlxPnts[iFlxPnt];
      
      *m_avgGrad[iGrad] -= dampFactor*dGradVarXNormal;
      
      *m_avgGradAV[iGrad] -= dampFactor*dGradVarXNormalAV;
    }
    
    // prepare flux computation
    prepareFluxComputation();
     
    // compute the Riemann flux
    computeFlux(m_avgSol,m_avgGrad,m_unitNormalFlxPnts[iFlxPnt],0,m_flxPntRiemannFlux[iFlxPnt]);
    
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

CFreal LLAVDiffFluxReconstructionNS::computeViscCoef(RealVector* state)
{
  CFreal result = 1.0;
  
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void LLAVDiffFluxReconstructionNS::computeSmoothness()
{ 
  CFreal sNum = 0.0;
  
  CFreal sDenom = 0.0;
  
  // get datahandle
  DataHandle< CFreal > monPhysVar = socket_monPhysVar.getDataHandle();
  
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    CFreal stateP = 0.0;
    CFreal diffStatesPPMinOne = 0.0;
    
    if (m_monitoredPhysVar < m_pData.size())
    {
      RealVector statePMinOne = *((*m_cellStates)[iSol]->getData()) - m_statesPMinOne[iSol];
      
      const bool RhoivtTv = getMethodData().getUpdateVarStr() == "RhoivtTv";
      
      if (RhoivtTv)
      {
        m_eulerVarSet2->computePhysicalData(*((*m_cellStates)[iSol]),m_pData);
        m_eulerVarSet2->computePhysicalData(m_statesPMinOne[iSol],m_pData2);
      }
      else 
      {
	m_eulerVarSet->computePhysicalData(*((*m_cellStates)[iSol]),m_pData);
        m_eulerVarSet->computePhysicalData(m_statesPMinOne[iSol],m_pData2);
      }
      
      stateP = m_pData[m_monitoredPhysVar];
      diffStatesPPMinOne = stateP - m_pData2[m_monitoredPhysVar];
    
      monPhysVar[(((*m_cellStates)[iSol]))->getLocalID()] = stateP;
    }
    else
    {
      stateP = (*((*m_cellStates)[iSol]))[m_monitoredVar];
      diffStatesPPMinOne = stateP - m_statesPMinOne[iSol][m_monitoredVar];
    }
    
    sNum += diffStatesPPMinOne*diffStatesPPMinOne;
    sDenom += stateP*stateP;
  }
  const RealVector& coords = *((*m_cellNodes)[0]->getData());
  const CFreal d = sqrt(coords[0]*coords[0]+coords[1]*coords[1]);

  if (sNum <= MathTools::MathConsts::CFrealEps() || sDenom <= MathTools::MathConsts::CFrealEps() ) //  || d < 1.2 || coords[0] < 0.1
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

CFreal LLAVDiffFluxReconstructionNS::computePeclet()
{ 
  CFreal result = m_peclet;
  
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void LLAVDiffFluxReconstructionNS::computeFlxPntStatesAndGrads()
{
  // loop over flx pnts to extrapolate the states to the flux points
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {     
    // local flux point indices in the left and right cell
    const CFuint flxPntIdxL = (*m_faceFlxPntConnPerOrient)[m_orient][LEFT][iFlxPnt];
    const CFuint flxPntIdxR = (*m_faceFlxPntConnPerOrient)[m_orient][RIGHT][iFlxPnt];
    
    // reset states in flx pnt
    *(m_cellStatesFlxPnt[LEFT][iFlxPnt]) = 0.0;
    *(m_cellStatesFlxPnt[RIGHT][iFlxPnt]) = 0.0;

    // reset the grads in the flx pnts
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      *(m_cellGradFlxPnt[LEFT][iFlxPnt][iVar]) = 0.0;
      *(m_cellGradFlxPnt[RIGHT][iFlxPnt][iVar]) = 0.0;
      *(m_cellGradFlxPntAV[LEFT][iFlxPnt][iVar]) = 0.0;
      *(m_cellGradFlxPntAV[RIGHT][iFlxPnt][iVar]) = 0.0;
    }

    // extrapolate the left and right states to the flx pnts
    for (CFuint iSol = 0; iSol < m_nbrSolDep; ++iSol)
    {
      const CFuint solIdxL = (*m_flxSolDep)[flxPntIdxL][iSol];
      const CFuint solIdxR = (*m_flxSolDep)[flxPntIdxR][iSol];
 
      // add the contributions of the current sol pnt
      *(m_cellStatesFlxPnt[LEFT][iFlxPnt]) += (*m_solPolyValsAtFlxPnts)[flxPntIdxL][solIdxL]*(*((*(m_states[LEFT]))[solIdxL]));
      *(m_cellStatesFlxPnt[RIGHT][iFlxPnt]) += (*m_solPolyValsAtFlxPnts)[flxPntIdxR][solIdxR]*(*((*(m_states[RIGHT]))[solIdxR]));

      for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
      {
        *(m_cellGradFlxPnt[LEFT][iFlxPnt][iVar]) += (*m_solPolyValsAtFlxPnts)[flxPntIdxL][solIdxL]*((*(m_cellGrads[LEFT][solIdxL]))[iVar]);
	*(m_cellGradFlxPnt[RIGHT][iFlxPnt][iVar]) += (*m_solPolyValsAtFlxPnts)[flxPntIdxR][solIdxR]*((*(m_cellGrads[RIGHT][solIdxR]))[iVar]);
        *(m_cellGradFlxPntAV[LEFT][iFlxPnt][iVar]) += (*m_solPolyValsAtFlxPnts)[flxPntIdxL][solIdxL]*((*(m_cellGradsAV[LEFT][solIdxL]))[iVar]);
	*(m_cellGradFlxPntAV[RIGHT][iFlxPnt][iVar]) += (*m_solPolyValsAtFlxPnts)[flxPntIdxR][solIdxR]*((*(m_cellGradsAV[RIGHT][solIdxR]))[iVar]);
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void LLAVDiffFluxReconstructionNS::computeDivDiscontFlx(vector< RealVector >& residuals)
{
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
      *(m_tempGrad[iVar]) = (*(m_cellGrads[0][iSolPnt]))[iVar];
    }
    
    m_avgSol = *((*m_cellStates)[iSolPnt]->getData());

    // calculate the discontinuous flux projected on x, y, z-directions
    for (CFuint iDim = 0; iDim < m_dim+m_ndimplus; ++iDim)
    {
      // add diff part
      computeFlux(m_avgSol,m_tempGrad,m_cellFluxProjVects[iDim][iSolPnt],0,m_contFlx[iSolPnt][iDim]);

      // add artificial part
      for (CFuint iDim2 = 0; iDim2 < m_dim; ++iDim2)
      {
        for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
        {
          m_contFlx[iSolPnt][iDim][iVar] += m_solEpsilons[iSolPnt]*(((*(m_cellGradsAV[0][iSolPnt]))[iVar])[iDim2])*m_cellFluxProjVects[iDim][iSolPnt][iDim2];
        }
      }
    }

    // extrapolate to flx pnts
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
    residuals[iSolPnt] = 0.0;

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
          residuals[iSolPnt][iEq] += polyCoef*(m_contFlx[jSolIdx][iDir+m_ndimplus][iEq]);
	}
      }
    }
  }

  // add the contribution of the faces
  const CFuint nbrFaces = m_cell->nbNeighborGeos();

  for (CFuint iFace = 0; iFace < nbrFaces; ++iFace)
  {
    if (!((*m_isFaceOnBoundaryCell)[iFace]))
    {
      for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
      {
        const CFuint currFlxIdx = (*m_faceFlxPntConn)[iFace][iFlxPnt];

        for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolDep; ++iSolPnt)
        {
          const CFuint solIdx = (*m_flxSolDep)[currFlxIdx][iSolPnt];
          const CFreal divh = m_corrFctDiv[solIdx][currFlxIdx];
   
          // Fill in the corrections
          for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
          {
            residuals[solIdx][iVar] += -m_extrapolatedFluxes[currFlxIdx][iVar] * divh; 
          }
        }
      }
    }
    else
    {
      m_faceNodes = (*m_faces)[iFace]->getNodes();
      m_face = (*m_faces)[iFace];
      m_cellNodes = m_cell->getNodes();
      
      // get the datahandle of the update coefficients
      DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
  
      // compute flux point coordinates
      for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
      {
        m_flxPntCoords[iFlx] = m_face->computeCoordFromMappedCoord((*m_flxLocalCoords)[iFlx]);	
      }
          
      // compute face Jacobian vectors
      m_faceJacobVecs = m_face->computeFaceJacobDetVectorAtMappedCoords(*m_flxLocalCoords);
  
      // get face Jacobian vector sizes in the flux points
      DataHandle< vector< CFreal > > faceJacobVecSizeFaceFlxPnts = socket_faceJacobVecSizeFaceFlxPnts.getDataHandle();
  
      // Loop over flux points to compute the unit normals
      for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
      {
        // get face Jacobian vector size
        CFreal faceJacobVecAbsSizeFlxPnts = faceJacobVecSizeFaceFlxPnts[m_face->getID()][iFlxPnt];
	
	// set face Jacobian vector size with sign depending on mapped coordinate direction
        m_faceJacobVecSizeFlxPnts2[iFlxPnt] = faceJacobVecAbsSizeFlxPnts*((*m_faceLocalDir)[iFace]);
 
	// set unit normal vector
        m_unitNormalFlxPnts2[iFlxPnt] = (m_faceJacobVecs[iFlxPnt]/faceJacobVecAbsSizeFlxPnts);
      }
	
      for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
      {
        const CFuint currFlxIdx = (*m_faceFlxPntConn)[iFace][iFlxPnt];
    
        // reset the grads in the flx pnts
        for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
        {
          *(m_cellGradFlxPnt[0][iFlxPnt][iVar]) = 0.0;
        }
        
        *(m_cellStatesFlxPnt[0][iFlxPnt]) = 0.0;

        // loop over the sol pnts to compute the states and grads in the flx pnts
        for (CFuint iSol = 0; iSol < m_nbrSolDep; ++iSol)
        {
          const CFuint solIdx = (*m_flxSolDep)[currFlxIdx][iSol];
          const CFreal coeff = (*m_solPolyValsAtFlxPnts)[currFlxIdx][solIdx];

	  *(m_cellStatesFlxPnt[0][iFlxPnt]) += coeff*(*((*(m_cellStates))[solIdx]));
	  
          for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
          {
            *(m_cellGradFlxPnt[0][iFlxPnt][iVar]) += coeff*((*(m_cellGradsAV[0][solIdx]))[iVar]);
          }
        }
      }
	
      // compute ghost gradients
      if ((getMethodData().getUpdateVarStr() == "Cons" || getMethodData().getUpdateVarStr() == "RhoivtTv") && getMethodData().hasDiffTerm())
      {
	for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
        {
	  for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
          {
	    *(m_flxPntGhostGrads[iFlxPnt][iVar]) = *(m_cellGradFlxPnt[0][iFlxPnt][iVar]);
	  }
	}
      }
      else
      {
	(*m_bcStateComputers)[(*m_faceBCIdxCell)[iFace]]->computeGhostGradients(m_cellGradFlxPnt[0],m_flxPntGhostGrads,m_unitNormalFlxPnts2,m_flxPntCoords);
      }
	
      for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
      {
	const CFuint currFlxIdx = (*m_faceFlxPntConn)[iFace][iFlxPnt];
	
	CFreal epsilon = 0.0;
	
	// loop over the sol pnts to compute the states and grads in the flx pnts
        for (CFuint iNode = 0; iNode < m_faceNodes->size(); ++iNode)
        {
	  for (CFuint iNodeCell = 0; iNodeCell < m_nbrCornerNodes; ++iNodeCell)
          {
	    if ((*m_faceNodes)[iNode]->getLocalID() == (*m_cellNodes)[iNodeCell]->getLocalID())
	    {
              //const CFuint nodeIdx = (*m_faceNodes)[iNode]->getLocalID();
	      // get node local index
              const CFuint nodeIdx = (*m_cellNodesConn)(m_cell->getID(),iNodeCell);
	  
              epsilon += m_nodePolyValsAtFlxPnts[currFlxIdx][iNodeCell]*m_nodeEpsilons[nodeIdx]/m_nbNodeNeighbors[nodeIdx];
	    }
	  }
        }
	
	if (!m_jacob && m_addUpdCoeff)
	{
	  // adding updateCoeff
	  CFreal visc = 1.0;
  
          m_waveSpeedUpd[0] = 0.0;

          const CFreal jacobXJacobXIntCoef = m_faceJacobVecSizeFlxPnts2[iFlxPnt]*
                                             m_faceJacobVecSizeFlxPnts2[iFlxPnt]*
                                             (*m_faceIntegrationCoefs)[iFlxPnt]*
                                             m_cflConvDiffRatio;
          //const CFreal rho = (*(m_cellStatesFlxPnt[0][iFlxPnt]))[0];
          const CFreal viscCoef = computeViscCoef(m_cellStatesFlxPnt[0][iFlxPnt]);
          visc = epsilon*viscCoef;
      
          // transform update states to physical data to calculate eigenvalues
          m_waveSpeedUpd[0] += visc*jacobXJacobXIntCoef/m_cell->computeVolume();

          // loop over the sol pnts of both sides to update the wave speeds
          for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
          {
            const CFuint solID = (*m_cellStates)[iSol]->getLocalID();
            updateCoeff[solID] += m_waveSpeedUpd[0];
          }
          //if (m_waveSpeedUpd[0] > 10.0) CFLog(INFO, "wvspLLAVBnd: " << m_waveSpeedUpd[0] << "\n");
	}
	
        // compute the average sol and grad to use the BR2 scheme
        for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
        {
	  if (m_cell->getID() == 1092) CFLog(VERBOSE, "var: " << iVar << ", grad: " << *(m_cellGradFlxPnt[0][iFlxPnt][iVar]) << ", ghost: " << *(m_flxPntGhostGrads[iFlxPnt][iVar]) << "\n");
          *(m_avgGrad[iVar]) = (*(m_cellGradFlxPnt[0][iFlxPnt][iVar]) + *(m_flxPntGhostGrads[iFlxPnt][iVar]))/2.0;
        }
              
        m_flxPntRiemannFlux[iFlxPnt] = 0.0;
	      
        for (CFuint iDim = 0; iDim < m_dim; ++iDim)
        {
          for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
          {
            m_flxPntRiemannFlux[iFlxPnt][iVar] += epsilon*((*(m_avgGrad[iVar]))[iDim])*m_unitNormalFlxPnts2[iFlxPnt][iDim];
	    if (m_cell->getID() == 1092) CFLog(VERBOSE, "avgrad: " << (*(m_avgGrad[iVar]))[iDim] << "\n");
          }
        }

        // compute FI in the mapped coord frame
        m_cellFlx[0][iFlxPnt] = (m_flxPntRiemannFlux[iFlxPnt])*m_faceJacobVecSizeFlxPnts2[iFlxPnt]; 
	if (m_cell->getID() == 1092) CFLog(VERBOSE, "riemannunit: " << m_flxPntRiemannFlux[iFlxPnt] << "jacob: " << m_faceJacobVecSizeFlxPnts2[iFlxPnt] << "\n");
	
	for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolDep; ++iSolPnt)
        {  
          const CFuint solIdx = (*m_flxSolDep)[currFlxIdx][iSolPnt];
          const CFreal divh = m_corrFctDiv[solIdx][currFlxIdx];
   
          // Fill in the corrections
          for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
          {
            residuals[solIdx][iVar] += (m_cellFlx[0][iFlxPnt][iVar] - m_extrapolatedFluxes[currFlxIdx][iVar]) * divh; 
          }
        }
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void LLAVDiffFluxReconstructionNS::setup()
{
  CFAUTOTRACE;
  
  // setup parent class
  LLAVDiffFluxReconstruction::setup();
  
  const bool RhoivtTv = getMethodData().getUpdateVarStr() == "RhoivtTv";
  
  // get damping coeff
  m_dampCoeff = getMethodData().getDiffDampCoefficient();
  
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
  
  // if needed, get the MS 
  
  m_gradsBackUp.resize(2);
  m_gradsBackUp[LEFT].resize(m_nbrSolPnts);
  m_gradsBackUp[RIGHT].resize(m_nbrSolPnts);

  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  { 
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      m_gradsBackUp[LEFT][iSol].push_back(new RealVector(m_dim));
      m_gradsBackUp[RIGHT][iSol].push_back(new RealVector(m_dim));
    }
  }
  
  // get the diffusive varset
  m_diffusiveVarSet = (getMethodData().getDiffusiveVar()).d_castTo< NavierStokesVarSet >();
  
  m_tempGradTermL.resize(m_nbrEqs,m_nbrFaceFlxPnts);
  m_tempGradTermR.resize(m_nbrEqs,m_nbrFaceFlxPnts);
  
  m_tempStatesL.resize(m_nbrFaceFlxPnts);
  m_tempStatesR.resize(m_nbrFaceFlxPnts);
  
  m_cellGradsAV.resize(2);
  m_cellGradsAV[LEFT].resize(m_nbrSolPnts);
  m_cellGradsAV[RIGHT].resize(m_nbrSolPnts);
  
  m_cellGradFlxPntAV.resize(2);
  m_avgGradAV.reserve(m_nbrEqs);
  m_cellGradFlxPntAV[LEFT].resize(m_nbrFaceFlxPnts);
  m_cellGradFlxPntAV[RIGHT].resize(m_nbrFaceFlxPnts);
  
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  { 
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      m_cellGradFlxPntAV[LEFT][iFlx].push_back(new RealVector(m_dim));
      m_cellGradFlxPntAV[RIGHT][iFlx].push_back(new RealVector(m_dim));
    }
  }
  
  for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
  {
    m_avgGradAV[iVar] = new RealVector(m_dim);
  }
}

//////////////////////////////////////////////////////////////////////////////

void LLAVDiffFluxReconstructionNS::unsetup()
{
  CFAUTOTRACE;
  
  for (CFuint iSol= 0; iSol < m_nbrSolPnts; ++iSol)
  {
    for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
    {
      deletePtr(m_gradsBackUp[LEFT][iSol][iGrad]);  
      deletePtr(m_gradsBackUp[RIGHT][iSol][iGrad]);
    }
    m_gradsBackUp[LEFT][iSol].clear();
    m_gradsBackUp[RIGHT][iSol].clear();
  }
  for (CFuint iFlx= 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
    {
      deletePtr(m_cellGradFlxPntAV[LEFT][iFlx][iGrad]);  
      deletePtr(m_cellGradFlxPntAV[RIGHT][iFlx][iGrad]);
    }
    m_cellGradFlxPntAV[LEFT][iFlx].clear();
    m_cellGradFlxPntAV[RIGHT][iFlx].clear();
  }
  
  m_gradsBackUp[LEFT].clear();
  m_gradsBackUp[RIGHT].clear();
  m_gradsBackUp.clear();
  
  for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
  {
    deletePtr(m_avgGradAV[iVar]); 
  }
  m_avgGradAV.clear();
  m_cellGradFlxPntAV[LEFT].clear();
  m_cellGradFlxPntAV[RIGHT].clear();
  m_cellGradFlxPntAV.clear();
  
  // unsetup parent class
  LLAVDiffFluxReconstruction::unsetup();
}


//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

