// Copyright (C) 2016 KU Leuven, Belgium
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

#include "FluxReconstructionNavierStokes/NSBndGradientComputer.hh"
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

MethodCommandProvider< NSBndGradientComputer,
		       FluxReconstructionSolverData,
		       FluxReconstructionNavierStokesModule >
NsBndGradientComputerProvider("ConvBndCorrectionsRHSNS");

MethodCommandProvider< NSBndGradientComputer,
		       FluxReconstructionSolverData,
		       FluxReconstructionNavierStokesModule >
NSBndGradientComputerTurbProvider("ConvBndCorrectionsRHSTurb");
  
//////////////////////////////////////////////////////////////////////////////
  
NSBndGradientComputer::NSBndGradientComputer(const std::string& name) :
  ConvBndCorrectionsRHSFluxReconstruction(name),
  m_diffusiveVarSet(CFNULL),
  m_tempGradTerm(),
  m_tempGradTermGhost(),
  m_tempStates(),
  m_tempStatesGhost()
{
}

//////////////////////////////////////////////////////////////////////////////

void NSBndGradientComputer::computeGradientBndFaceCorrections()
{ 
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_tempStates[iFlx] = m_cellStatesFlxPnt[iFlx]->getData();
    m_tempStatesGhost[iFlx] = m_flxPntGhostSol[iFlx]->getData();
  }
  
  m_diffusiveVarSet->setGradientVars(m_tempStates,m_tempGradTerm,m_nbrFaceFlxPnts);
  m_diffusiveVarSet->setGradientVars(m_tempStatesGhost,m_tempGradTermGhost,m_nbrFaceFlxPnts);
  
  // Loop over solution pnts to reset the grad updates
  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
  {
    // Loop over  variables
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      //set the grad updates to 0 
      m_gradUpdates[iSolPnt][iEq] = 0.0;
    }
  }
      
  // compute the face corrections to the gradients
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    const CFuint flxIdx = (*m_faceFlxPntConn)[m_orient][iFlx];

    m_nbrSolDep = ((*m_flxSolDep)[flxIdx]).size();
    // Loop over  variables
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
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
  
  // get the gradients
  DataHandle< vector< RealVector > > gradients = socket_gradients.getDataHandle();

  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    // get state ID
    const CFuint solID = (*m_cellStates)[iSol]->getLocalID();

    // update gradients
    for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
    {
      gradients[solID][iGrad] += m_gradUpdates[iSol][iGrad];
     
    }
  }
  
  // if needed, compute the gradients for the artificial viscosity
  if ((getMethodData().getUpdateVarStr() == "Cons") && getMethodData().hasArtificialViscosity())
  {
    // Loop over solution pnts to reset the grad updates
    for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
    {
      // Loop over  variables
      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      {
        //set the grad updates to 0 
        m_gradUpdates[iSolPnt][iEq] = 0.0;
      }
    }
      
    // compute the face corrections to the gradients
    for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
    {
      const CFuint flxIdx = (*m_faceFlxPntConn)[m_orient][iFlx];
      
      m_nbrSolDep = ((*m_flxSolDep)[flxIdx]).size();
      // Loop over  variables
      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      {
        const CFreal avgSol = ((*(m_cellStatesFlxPnt[iFlx]))[iEq]+(*(m_flxPntGhostSol[iFlx]))[iEq])/2.0;
	m_projectedCorr = (avgSol-(*(m_cellStatesFlxPnt[iFlx]))[iEq])*m_faceJacobVecSizeFlxPnts[iFlx]*m_unitNormalFlxPnts[iFlx];

        // Loop over solution pnts to calculate the grad updates
        for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolDep; ++iSolPnt)
        {
          const CFuint iSolIdx = (*m_flxSolDep)[flxIdx][iSolPnt];

	  /// @todo Check if this is also OK for triangles!!
	  m_gradUpdates[iSolIdx][iEq] += m_projectedCorr*m_corrFctDiv[iSolIdx][flxIdx];
        }
      }
    }
  
    // get the gradients
    DataHandle< vector< RealVector > > gradientsAV = socket_gradientsAV.getDataHandle();

    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      // get state ID
      const CFuint solID = (*m_cellStates)[iSol]->getLocalID();

      // update gradients
      for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
      {
        gradientsAV[solID][iGrad] += m_gradUpdates[iSol][iGrad];
      }
    }
  }
  else if ((getMethodData().getUpdateVarStr() == "Puvt" || getMethodData().getUpdateVarStr() == "RhoivtTv") && getMethodData().hasArtificialViscosity())
  {
    // Loop over solution pnts to calculate the grad updates
    for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
    {
      // Loop over  variables
      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      {
        //set the grad updates to 0 
        m_gradUpdates[iSolPnt][iEq] = 0.0;
      }
    }
      
    // compute the face corrections to the gradients
    for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
    {
      const CFuint flxIdx = (*m_faceFlxPntConn)[m_orient][iFlx];

      const RealVector transformedState = static_cast<RealVector&>(*m_updateToSolutionVecTrans->transform(m_cellStatesFlxPnt[iFlx]));
      const RealVector transformedGhostState = static_cast<RealVector&>(*m_updateToSolutionVecTrans->transform(m_flxPntGhostSol[iFlx]));
      
      m_nbrSolDep = ((*m_flxSolDep)[flxIdx]).size();
      // Loop over  variables
      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      {
	const CFreal avgSol = (transformedState[iEq]+transformedGhostState[iEq])/2.0;
	m_projectedCorr = (avgSol-transformedState[iEq])*m_faceJacobVecSizeFlxPnts[iFlx]*m_unitNormalFlxPnts[iFlx];
	  
        // Loop over solution pnts to calculate the grad updates
        for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolDep; ++iSolPnt)
        {
          const CFuint iSolIdx = (*m_flxSolDep)[flxIdx][iSolPnt];

	  /// @todo Check if this is also OK for triangles!!
	  m_gradUpdates[iSolIdx][iEq] += m_projectedCorr*m_corrFctDiv[iSolIdx][flxIdx];
        }
      }
    }
  
    // get the gradients
    DataHandle< vector< RealVector > > gradientsAV = socket_gradientsAV.getDataHandle();

    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      // get state ID
      const CFuint solID = (*m_cellStates)[iSol]->getLocalID();

      // update gradients
      for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
      {
        gradientsAV[solID][iGrad] += m_gradUpdates[iSol][iGrad];
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void NSBndGradientComputer::setup()
{
  CFAUTOTRACE;
  
  // setup parent class
  ConvBndCorrectionsRHSFluxReconstruction::setup();
  
  // get the diffusive varset
  m_diffusiveVarSet = (getMethodData().getDiffusiveVar()).d_castTo< NavierStokesVarSet >();
  
  m_updateToSolutionVecTrans = getMethodData().getUpdateToSolutionVecTrans();
  
  m_updateToSolutionVecTrans->setup(2);
  
  m_tempGradTerm.resize(m_nbrEqs,m_nbrFaceFlxPnts);
  m_tempGradTermGhost.resize(m_nbrEqs,m_nbrFaceFlxPnts);
  
  m_tempStates.resize(m_nbrFaceFlxPnts);
  m_tempStatesGhost.resize(m_nbrFaceFlxPnts);
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

