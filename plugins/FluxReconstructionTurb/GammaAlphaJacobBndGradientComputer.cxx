// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/MethodCommandProvider.hh"

#include "Framework/CFSide.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/BaseTerm.hh"

#include "MathTools/MathFunctions.hh"

#include "FluxReconstructionTurb/GammaAlphaJacobBndGradientComputer.hh"
#include "FluxReconstructionNavierStokes/FluxReconstructionNavierStokes.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"
#include "NavierStokes/NavierStokesVarSet.hh"

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

MethodCommandProvider< GammaAlphaJacobBndGradientComputer,
		       FluxReconstructionSolverData,
		       FluxReconstructionNavierStokesModule >
GammaAlphaJacobBndGradientComputerProvider("ConvBndCorrectionsRHSJacobGammaAlpha");
  
//////////////////////////////////////////////////////////////////////////////
  
GammaAlphaJacobBndGradientComputer::GammaAlphaJacobBndGradientComputer(const std::string& name) :
  NSJacobBndGradientComputer(name),
  socket_volumes("volumes"),
  socket_solPntNormals("solPntNormals"),
  m_gradDummy(),
  m_uGrad(),
  m_vGrad(),
  m_wGrad(),
  m_cellFluxProjVects(),
  m_solSolDep(CFNULL),
  m_nbrSolSolDep(),
  m_solPolyDerivAtSolPnts(CFNULL),
  m_tempGradTermIntCell(),
  m_tempStatesIntCell()
{
}

//////////////////////////////////////////////////////////////////////////////

void GammaAlphaJacobBndGradientComputer::setup()
{
  CFAUTOTRACE;
  
  // setup parent class
  NSJacobBndGradientComputer::setup();
  
  m_uGrad.resize(m_dim);
  m_vGrad.resize(m_dim);
  m_wGrad.resize(m_dim);
  
  m_cellFluxProjVects.resize(m_dim);
  for (CFuint iDim = 0; iDim < m_dim; ++iDim)
  {
    m_cellFluxProjVects[iDim].resize(m_nbrSolPnts);
    for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
    {
      m_cellFluxProjVects[iDim][iSolPnt].resize(m_dim);
    }
  }
  
  // get the local FR data
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
  cf_assert(frLocalData.size() > 0);
  // for now, there should be only one type of element
  cf_assert(frLocalData.size() == 1);
  
  // get the coefs for derivation of the states in the sol pnts
  m_solPolyDerivAtSolPnts = frLocalData[0]->getCoefSolPolyDerivInSolPnts();
  
  m_solSolDep = frLocalData[0]->getSolPntSolDependency();
  
  m_nbrSolSolDep = ((*m_solSolDep)[0]).size();
  
  m_tempGradTermIntCell.resize(m_nbrEqs,m_nbrSolPnts);

  m_tempStatesIntCell.resize(m_nbrSolPnts);
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > >
GammaAlphaJacobBndGradientComputer::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result = NSJacobBndGradientComputer::needsSockets();
  result.push_back(&socket_solPntNormals);
  result.push_back(&socket_volumes);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void GammaAlphaJacobBndGradientComputer::computeFlxPntStates()
{ 
  // Loop over flux points to extrapolate the states to the flux points
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {
    // reset the extrapolated states
    *(m_cellStatesFlxPnt[iFlxPnt]) = 0.0;
    
    // get current flx pnt idx
    const CFuint currFlxIdx = (*m_faceFlxPntConn)[m_orient][iFlxPnt];
    
    // extrapolate the states to current flx pnt
    for (CFuint iSol = 0; iSol < m_nbrSolDep; ++iSol)
    {
      const CFuint solIdx = (*m_flxSolDep)[currFlxIdx][iSol];

      *(m_cellStatesFlxPnt[iFlxPnt]) += (*m_solPolyValsAtFlxPnts)[currFlxIdx][solIdx]*(*((*m_cellStates)[solIdx]));
    }
    
    m_bcStateComputer->setTransitionCriterion(iFlxPnt,false); 
  }
  
  // pre-compute the ghost states assuming no transition
  m_bcStateComputer->computeGhostStates(m_cellStatesFlxPnt,m_flxPntGhostSol,m_unitNormalFlxPnts,m_flxPntCoords);

  // get the sol pnt normals
  DataHandle< CFreal > solPntNormals = socket_solPntNormals.getDataHandle();

  // store the sol pnt normals
  for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
  {
    const CFuint solID = (*(m_cellStates))[iState]->getLocalID();
      
    for (CFuint iDim = 0; iDim < m_dim; ++iDim)
    {
      for (CFuint jDim = 0; jDim < m_dim; ++jDim)
      {
        m_cellFluxProjVects[iDim][iState][jDim] = solPntNormals[solID*m_dim*m_dim+iDim*m_dim+jDim];
      }
    }
  }

  SafePtr< NavierStokesVarSet > navierStokesVarSet = m_diffusiveVarSet.d_castTo< NavierStokesVarSet >();
  
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_tempStates[iFlx] = m_cellStatesFlxPnt[iFlx]->getData();
    m_tempStatesGhost[iFlx] = m_flxPntGhostSol[iFlx]->getData();
  }
  
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    m_tempStatesIntCell[iSol] = (*m_cellStates)[iSol]->getData();
    
    //set the grad updates to 0 
    for (CFuint iDim = 0; iDim < m_dim; ++iDim)
    {
      m_gradUpdates[iSol][1+iDim] = 0.0;
    }
  }
  
  m_diffusiveVarSet->setGradientVars(m_tempStates,m_tempGradTerm,m_nbrFaceFlxPnts);
  m_diffusiveVarSet->setGradientVars(m_tempStatesGhost,m_tempGradTermGhost,m_nbrFaceFlxPnts);
  m_diffusiveVarSet->setGradientVars(m_tempStatesIntCell,m_tempGradTermIntCell,m_nbrSolPnts);

  // compute the face corrections to the gradients
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    const CFuint flxIdx = (*m_faceFlxPntConn)[m_orient][iFlx];
    // Loop over  variables
    for (CFuint iEq = 1; iEq < 1+m_dim; ++iEq)
    {
      const CFreal avgSol = (m_tempGradTerm(iEq,iFlx)+m_tempGradTermGhost(iEq,iFlx))/2.0;
      m_projectedCorr = (avgSol-m_tempGradTerm(iEq,iFlx))*m_faceJacobVecSizeFlxPnts[iFlx]*m_unitNormalFlxPnts[iFlx];

      // Loop over solution pnts to calculate the grad updates
      for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolDep; ++iSolPnt)
      {
        const CFuint iSolIdx = (*m_flxSolDep)[flxIdx][iSolPnt];

	/// @todo Check if this is also OK for triangles!!
	m_gradUpdates[iSolIdx][iEq] += m_projectedCorr*m_corrFctDiv[iSolIdx][flxIdx];
      }
    }
  }

  // Loop over solution pnts to calculate the grad updates
  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
  {
    // Loop over  variables
    for (CFuint iEq = 1; iEq < 1+m_dim; ++iEq)
    {
      // Loop over gradient directions
      for (CFuint iDir = 0; iDir < m_dim; ++iDir)
      {

        m_projectedCorr = m_tempGradTermIntCell(iEq,iSolPnt) * m_cellFluxProjVects[iDir][iSolPnt];

        // Loop over solution pnts to count factor of all sol pnt polys
        for (CFuint jSolPnt = 0; jSolPnt < m_nbrSolSolDep; ++jSolPnt)
        { 
          const CFuint jSolIdx = (*m_solSolDep)[iSolPnt][jSolPnt];

          // compute the grad updates
          m_gradUpdates[jSolIdx][iEq] += (*m_solPolyDerivAtSolPnts)[jSolIdx][iDir][iSolPnt]*m_projectedCorr;
	}
      }
    }
  }

  // get the volumes
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();

  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    // get state ID
    const CFuint solID = (*m_cellStates)[iSol]->getLocalID();

    // inverse Jacobian determinant
    const CFreal invJacobDet = 1.0/volumes[solID];

    // update gradients
    for (CFuint iGrad = 1; iGrad < 1+m_dim; ++iGrad)
    {
      m_gradUpdates[iSol][iGrad] *= invJacobDet;
    }
  }

  // Loop over flux points to check the transition criterion
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {
    // reset the extrapolated gradients
    m_uGrad = 0.0;
    m_vGrad = 0.0;
    m_wGrad = 0.0;
    
    // get current flx pnt idx
    const CFuint currFlxIdx = (*m_faceFlxPntConn)[m_orient][iFlxPnt];
    
    // extrapolate the states to current flx pnt
    if (m_dim == 2)
    {
    for (CFuint iSol = 0; iSol < m_nbrSolDep; ++iSol)
    {
      const CFuint solIdx = (*m_flxSolDep)[currFlxIdx][iSol];

      m_uGrad += (*m_solPolyValsAtFlxPnts)[currFlxIdx][solIdx]*m_gradUpdates[solIdx][1];
      m_vGrad += (*m_solPolyValsAtFlxPnts)[currFlxIdx][solIdx]*m_gradUpdates[solIdx][2];
    }
    }
    else
    {
    for (CFuint iSol = 0; iSol < m_nbrSolDep; ++iSol)
    {
      const CFuint solIdx = (*m_flxSolDep)[currFlxIdx][iSol];

      m_uGrad += (*m_solPolyValsAtFlxPnts)[currFlxIdx][solIdx]*m_gradUpdates[solIdx][1];
      m_vGrad += (*m_solPolyValsAtFlxPnts)[currFlxIdx][solIdx]*m_gradUpdates[solIdx][2];
      m_wGrad += (*m_solPolyValsAtFlxPnts)[currFlxIdx][solIdx]*m_gradUpdates[solIdx][3];
    }
    }
    
    const CFreal rho = navierStokesVarSet->getDensity(*(m_cellStatesFlxPnt[iFlxPnt]));
    
    const CFreal muTot = navierStokesVarSet->getDynViscosity(*(m_cellStatesFlxPnt[iFlxPnt]), m_gradDummy);
    
//    const CFreal tau = muTot*(m_unitNormalFlxPnts[iFlxPnt][YY]*MathFunctions::innerProd(m_uGrad,m_unitNormalFlxPnts[iFlxPnt]) -
//                       m_unitNormalFlxPnts[iFlxPnt][XX]*MathFunctions::innerProd(m_vGrad,m_unitNormalFlxPnts[iFlxPnt]));
    
    CFreal tau = 0.0;
    
    if (m_dim == 2)
    {
      tau = muTot*(m_unitNormalFlxPnts[iFlxPnt][YY]*MathFunctions::innerProd(m_uGrad,m_unitNormalFlxPnts[iFlxPnt]) -
            m_unitNormalFlxPnts[iFlxPnt][XX]*MathFunctions::innerProd(m_vGrad,m_unitNormalFlxPnts[iFlxPnt]));
    }
    else if (fabs(m_unitNormalFlxPnts[iFlxPnt][ZZ]) <= fabs(m_unitNormalFlxPnts[iFlxPnt][XX]))
    {
      const CFreal tauT1 = muTot*(m_unitNormalFlxPnts[iFlxPnt][YY]*MathFunctions::innerProd(m_uGrad,m_unitNormalFlxPnts[iFlxPnt]) -
                           m_unitNormalFlxPnts[iFlxPnt][XX]*MathFunctions::innerProd(m_vGrad,m_unitNormalFlxPnts[iFlxPnt]));
      
      const CFreal tauT2 = muTot*(m_unitNormalFlxPnts[iFlxPnt][XX]*m_unitNormalFlxPnts[iFlxPnt][ZZ]*MathFunctions::innerProd(m_uGrad,m_unitNormalFlxPnts[iFlxPnt]) +
                           m_unitNormalFlxPnts[iFlxPnt][YY]*m_unitNormalFlxPnts[iFlxPnt][ZZ]*MathFunctions::innerProd(m_vGrad,m_unitNormalFlxPnts[iFlxPnt]) - 
                           (m_unitNormalFlxPnts[iFlxPnt][XX]*m_unitNormalFlxPnts[iFlxPnt][XX]+m_unitNormalFlxPnts[iFlxPnt][YY]*m_unitNormalFlxPnts[iFlxPnt][YY])*MathFunctions::innerProd(m_wGrad,m_unitNormalFlxPnts[iFlxPnt]));
      tau = sqrt(tauT1*tauT1+tauT2*tauT2);
    }
    else
    {
      const CFreal tauT1 = muTot*(-m_unitNormalFlxPnts[iFlxPnt][ZZ]*MathFunctions::innerProd(m_vGrad,m_unitNormalFlxPnts[iFlxPnt]) +
                           m_unitNormalFlxPnts[iFlxPnt][YY]*MathFunctions::innerProd(m_wGrad,m_unitNormalFlxPnts[iFlxPnt]));
      
      const CFreal tauT2 = muTot*(-m_unitNormalFlxPnts[iFlxPnt][XX]*m_unitNormalFlxPnts[iFlxPnt][YY]*MathFunctions::innerProd(m_vGrad,m_unitNormalFlxPnts[iFlxPnt]) -
                           m_unitNormalFlxPnts[iFlxPnt][XX]*m_unitNormalFlxPnts[iFlxPnt][ZZ]*MathFunctions::innerProd(m_wGrad,m_unitNormalFlxPnts[iFlxPnt]) + 
                           (m_unitNormalFlxPnts[iFlxPnt][YY]*m_unitNormalFlxPnts[iFlxPnt][YY]+m_unitNormalFlxPnts[iFlxPnt][ZZ]*m_unitNormalFlxPnts[iFlxPnt][ZZ])*MathFunctions::innerProd(m_uGrad,m_unitNormalFlxPnts[iFlxPnt]));
      tau = sqrt(tauT1*tauT1+tauT2*tauT2);
    }
    
    const CFreal tauCrit = tau/sqrt(rho*muTot);
    
    if (tauCrit <= (*(m_cellStatesFlxPnt[iFlxPnt]))[5+m_dim])
    {
      m_bcStateComputer->setTransitionCriterion(iFlxPnt,true);
    }
    else
    {
      m_bcStateComputer->setTransitionCriterion(iFlxPnt,false); 
    }
  }

  // compute ghost states
  m_bcStateComputer->computeGhostStates(m_cellStatesFlxPnt,m_flxPntGhostSol,m_unitNormalFlxPnts,m_flxPntCoords);
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

