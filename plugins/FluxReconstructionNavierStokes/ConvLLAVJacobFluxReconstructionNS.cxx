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

#include "FluxReconstructionNavierStokes/ConvLLAVJacobFluxReconstructionNS.hh"
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

MethodCommandProvider< ConvLLAVJacobFluxReconstructionNS,
		       FluxReconstructionSolverData,
		       FluxReconstructionNavierStokesModule >
convLLAVRHSJacobNSFluxReconstructionProvider("ConvLLAVRHSJacobNS");
  
//////////////////////////////////////////////////////////////////////////////
  
ConvLLAVJacobFluxReconstructionNS::ConvLLAVJacobFluxReconstructionNS(const std::string& name) :
  ConvLLAVJacobFluxReconstruction(name),
  m_eulerVarSet(CFNULL),
  m_msEulerTerm(CFNULL),
  m_nbrSpecies(),
  m_pData(),
  m_pData2(),
  m_eulerVarSet2(CFNULL)
{
}

//////////////////////////////////////////////////////////////////////////////

void ConvLLAVJacobFluxReconstructionNS::configure ( Config::ConfigArgs& args )
{
  ConvLLAVJacobFluxReconstruction::configure(args);
} 

//////////////////////////////////////////////////////////////////////////////

void ConvLLAVJacobFluxReconstructionNS::computeInterfaceFlxCorrection()
{
  // Loop over the flux points to calculate FI
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  { 
    // compute the average sol and grad to use the BR2 scheme
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      *(m_avgGradAV[iVar]) = (*(m_cellGradFlxPntAV[LEFT][iFlxPnt][iVar]) + *(m_cellGradFlxPntAV[RIGHT][iFlxPnt][iVar]))/2.0;
       
      m_avgSol[iVar] = ((*(m_cellStatesFlxPnt[LEFT][iFlxPnt]))[iVar] + (*(m_cellStatesFlxPnt[RIGHT][iFlxPnt]))[iVar])/2.0; 
    }

    // damping factor
    const CFreal dampFactor = m_dampCoeffDiff*m_faceInvCharLengths[iFlxPnt];

    // compute averaged (damped) gradients
    for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
    { 
      // compute damping term for LLAV
      const RealVector dGradVarXNormalAV = ((*(m_cellStatesFlxPnt[LEFT][iFlxPnt]))[iGrad] - (*(m_cellStatesFlxPnt[RIGHT][iFlxPnt]))[iGrad])*m_unitNormalFlxPnts[iFlxPnt];
      *m_avgGradAV[iGrad] -= dampFactor*dGradVarXNormalAV;
    }
    
    prepareFluxComputation();
     
    // compute the diff flux
    //computeFlux(m_avgSol,m_avgGrad,m_unitNormalFlxPnts[iFlxPnt],0,m_flxPntRiemannFluxDiff[iFlxPnt]);
    
    //m_flxPntRiemannFluxDiffConv[iFlxPnt] = m_flxPntRiemannFluxDiff[iFlxPnt];
    
    // compute the convective riemann flux
    m_flxPntRiemannFluxDiffConv[iFlxPnt] = -m_riemannFluxComputer->computeFlux(*(m_cellStatesFlxPnt[LEFT][iFlxPnt]),
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

void ConvLLAVJacobFluxReconstructionNS::computeRiemannFluxJacobianNum(const CFreal resFactor)
{
  CFLog(VERBOSE, "NS computeRiemannFluxJacobianNum\n");
      
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
        
//        // compute the average grad to use the BR2 scheme
//        for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
//        {
//          *(m_avgGrad[iVar]) = (*(m_cellGradFlxPnt[LEFT][iFlxPnt][iVar]) + *(m_cellGradFlxPnt[RIGHT][iFlxPnt][iVar]))/2.0;          
//        }

//        // compute averaged (damped) gradients
//        for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
//        {
//          // compute damping term
//          const RealVector dGradVarXNormal = (m_tempGradTermL(iGrad,iFlxPnt) - m_tempGradTermR(iGrad,iFlxPnt))*m_unitNormalFlxPnts[iFlxPnt];
//          *m_avgGrad[iGrad] -= dampFactor*dGradVarXNormal;
//        }
        
        // compute the average sol
        for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
        {        
          m_avgSol[iVar] = (pertState[iVar] + (*(m_cellStatesFlxPnt[iOtherSide][iFlxPnt]))[iVar])/2.0; 
        }
    
        prepareFluxComputation();
     
        // compute diffusive flux
        //computeFlux(m_avgSol,m_avgGrad,m_unitNormalFlxPnts[iFlxPnt],0,m_flxPntRiemannFluxPert[iFlxPnt]);
        
        if (m_pertSide == LEFT)
        {
          // compute the convective riemann flux
          m_flxPntRiemannFluxPert[iFlxPnt] = -m_riemannFluxComputer->computeFlux(pertState,
									              *(m_cellStatesFlxPnt[RIGHT][iFlxPnt]),
									              m_unitNormalFlxPnts[iFlxPnt]);
        }
        else
        {
          // compute the convective riemann flux
          m_flxPntRiemannFluxPert[iFlxPnt] = -m_riemannFluxComputer->computeFlux(*(m_cellStatesFlxPnt[LEFT][iFlxPnt]),
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
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvLLAVJacobFluxReconstructionNS::prepareFluxComputation()
{
  //const bool isPerturb = this->getMethodData().isPerturb();
  //const CFuint iPerturbVar = this->getMethodData().iPerturbVar();

  //m_diffusiveVarSetNS->setComposition(m_avgSol, isPerturb, iPerturbVar);
}

//////////////////////////////////////////////////////////////////////////////

void ConvLLAVJacobFluxReconstructionNS::computeGradients()
{
  // get the sol pnt normals
  DataHandle< CFreal > solPntNormals = socket_solPntNormals.getDataHandle();
  
  // get the volumes
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
  
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
        
      // Loop over  variables
      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      {
        // Loop over gradient directions
        for (CFuint iDir = 0; iDir < m_dim; ++iDir)
        {
          for (CFuint jDir = 0; jDir < m_dim; ++jDir)
          {
	    // project the state on a normal and reuse a RealVector variable of the class to store
	    m_projectedCorrL[jDir] = (*((*m_cellStates)[iSolPnt]))[iEq] * solPntNormals[solID*(m_dim+m_ndimplus)*m_dim+(iDir+m_ndimplus)*m_dim+jDir];
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

void ConvLLAVJacobFluxReconstructionNS::computeGradientFaceCorrections()
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

      // compute the face corrections to the gradients
      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      {
        const CFreal avgSol = ((*(m_cellStatesFlxPnt[LEFT][iFlx]))[iEq]+(*(m_cellStatesFlxPnt[RIGHT][iFlx]))[iEq])/2.0;
	m_projectedCorrL = (avgSol-(*(m_cellStatesFlxPnt[LEFT][iFlx]))[iEq])*m_faceJacobVecSizeFlxPnts[iFlx][LEFT]*m_unitNormalFlxPnts[iFlx];
	m_projectedCorrR = (avgSol-(*(m_cellStatesFlxPnt[RIGHT][iFlx]))[iEq])*m_faceJacobVecSizeFlxPnts[iFlx][RIGHT]*m_unitNormalFlxPnts[iFlx];

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

void ConvLLAVJacobFluxReconstructionNS::computeSmoothness()
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

void ConvLLAVJacobFluxReconstructionNS::computeSmoothness(const CFuint side)
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

void ConvLLAVJacobFluxReconstructionNS::computeEpsilon0()
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
  //if(m_cell->getID()==1) CFLog(INFO, "wvspd: " << h*wavespeed << ", P: " << peclet << ", dx: " << m_subcellRes << "\n");
  m_epsilon0 = max(h*wavespeed*(2.0/peclet - m_subcellRes/peclet),0.0);
}

//////////////////////////////////////////////////////////////////////////////

void ConvLLAVJacobFluxReconstructionNS::computeEpsilon0(const CFuint side)
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

void ConvLLAVJacobFluxReconstructionNS::computeEpsToStateJacobianAna()
{
  CFLog(VERBOSE, "NS computeEpsToStateJacobianAna\n");
    
  for (m_pertSide = 0; m_pertSide < 2; ++m_pertSide)
  {
    for (m_pertSol = 0; m_pertSol < m_nbrSolPnts; ++m_pertSol)
    {
      m_epsJacobian[m_pertSide][m_pertSol] = 0.0;
    }
  }
  
  const bool RhoivtTv = getMethodData().getUpdateVarStr() == "RhoivtTv";
  
  // get the cell volumes
  DataHandle< CFreal > cellVolumes = socket_cellVolumes.getDataHandle();
  
  const CFreal oneOverDim = 1./m_dim;
  
  const CFreal peclet = computePeclet();
  
  for (m_pertSide = 0; m_pertSide < 2; ++m_pertSide)
  { 
    CFreal h_f_S;
    
    const CFreal h = pow(cellVolumes[m_cells[m_pertSide]->getID()],oneOverDim);
    
    const CFreal h_f = h*(2.0/peclet - m_subcellRes/peclet);
    
    computeProjStates(m_statesPMinOne, m_pertSide);
    
    computeSmoothness(m_pertSide);
    
    CFreal sBefore;
    
    if (m_s < m_s0 - m_kappa)
    {
      h_f_S = 0.0;
      
      sBefore = 0.0;
    }
    else if (m_s > m_s0 + m_kappa)
    { 
      h_f_S = h_f/m_nbrSolPnts;
      
      sBefore = 1.0;
    }
    else
    {   
      sBefore = 0.5*(1.0 + sin(0.5*MathTools::MathConsts::CFrealPi()*(m_s-m_s0)/m_kappa));
              
      h_f_S = h_f/m_nbrSolPnts*sBefore;
    }
    
    CFreal wallFactor = 1.0;
    
    if (m_useWallCutOff)
    {
      // Get the wall distance
      DataHandle< CFreal > wallDist = socket_wallDistance.getDataHandle();
  
      CFreal centroidDistance = 0.0;
      
      for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
      {
        const CFuint stateID = (*m_states[m_pertSide])[iSol]->getLocalID();
        centroidDistance += wallDist[stateID];
      }
    
      centroidDistance /= m_nbrSolPnts;
    
      if (centroidDistance < m_wallCutOff) 
      {
        if (centroidDistance < 0.5*m_wallCutOff)
        {
          wallFactor = 0.0; 
        }
        else
        {
          wallFactor = 0.5*(1.0 + sin(0.5*MathTools::MathConsts::CFrealPi()*(centroidDistance-0.75*m_wallCutOff)/(0.25*m_wallCutOff)));
        }
      }
    }
    
    sBefore *= wallFactor;
      
    for (m_pertSol = 0; m_pertSol < m_nbrSolPnts; ++m_pertSol)
    {
      // dereference state
      State& pertState = *(*m_states[m_pertSide])[m_pertSol];
      
      if (RhoivtTv)
      {
        m_eulerVarSet2->computePhysicalData(*((*m_states[m_pertSide])[m_pertSol]),m_pData);
      }
      else
      {
        m_eulerVarSet->computePhysicalData(*((*m_states[m_pertSide])[m_pertSol]),m_pData);
      }
        
      const CFreal lambda = m_pData[EulerTerm::V] + m_pData[EulerTerm::A];
      cf_assert(m_pData[EulerTerm::V] > 0.0 && m_pData[EulerTerm::A] > 0.0);
          
      const CFreal h_f_lambda = h_f * lambda;
      
      for (m_pertVar = 0; m_pertVar < m_nbrEqs; ++m_pertVar)
      {
        const CFreal uBefore = pertState[m_pertVar];
        
        // perturb physical variable in state
        m_numJacob->perturb(m_pertVar,pertState[m_pertVar]);
        
        const CFreal eps = pertState[m_pertVar] - uBefore;
        
        if (RhoivtTv)
        {
          m_eulerVarSet2->computePhysicalData(*((*m_states[m_pertSide])[m_pertSol]),m_pData);
        }
        else
        {
          m_eulerVarSet->computePhysicalData(*((*m_states[m_pertSide])[m_pertSol]),m_pData);
        }
        
        const CFreal lambdaPert = m_pData[EulerTerm::V] + m_pData[EulerTerm::A];
          
        computeProjStates(m_statesPMinOne, m_pertSide);
    
        computeSmoothness(m_pertSide);
        
        CFreal sPert;
        
        if (m_s < m_s0 - m_kappa)
        { 
          sPert = 0.0;
        }
        else if (m_s > m_s0 + m_kappa)
        { 
          sPert = 1.0;
        }
        else
        {   
          sPert = 0.5*(1.0 + sin(0.5*MathTools::MathConsts::CFrealPi()*(m_s-m_s0)/m_kappa));
        }
        
        sPert *= wallFactor;

        // restore physical variable in state
        m_numJacob->restore(pertState[m_pertVar]);  
          
        m_epsJacobian[m_pertSide][m_pertSol][m_pertVar] += h_f_S*(lambdaPert - lambda)/eps;
          
        m_epsJacobian[m_pertSide][m_pertSol][m_pertVar] += h_f_lambda*(sPert - sBefore)/eps;
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvLLAVJacobFluxReconstructionNS::setup()
{
  CFAUTOTRACE;
  
  // setup parent class
  ConvLLAVJacobFluxReconstruction::setup();
  
  // get the diffusive varset
  //m_diffusiveVarSetNS = (getMethodData().getDiffusiveVar()).d_castTo< NavierStokesVarSet >();
  
  m_updateToSolutionVecTrans = getMethodData().getUpdateToSolutionVecTrans();
  
  m_updateToSolutionVecTrans->setup(2);
  
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
