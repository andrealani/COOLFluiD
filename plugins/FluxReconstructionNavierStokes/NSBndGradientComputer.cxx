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
  
//////////////////////////////////////////////////////////////////////////////
  
NSBndGradientComputer::NSBndGradientComputer(const std::string& name) :
  ConvBndCorrectionsRHSFluxReconstruction(name),
  m_diffusiveVarSet(CFNULL)
{
}

//////////////////////////////////////////////////////////////////////////////

void NSBndGradientComputer::computeGradientBndFaceCorrections()
{ 
  //ConvBndCorrectionsRHSFluxReconstruction::computeGradientBndFaceCorrections();
  // get the diffusive varset
  m_diffusiveVarSet = getMethodData().getDiffusiveVar();
  
  SafePtr< NavierStokesVarSet > navierStokesVarSet = m_diffusiveVarSet.d_castTo< NavierStokesVarSet >();

  RealMatrix tempGradTerm(m_nbrEqs,m_nbrFaceFlxPnts);
  RealMatrix tempGhostGradTerm(m_nbrEqs,m_nbrFaceFlxPnts);
  vector< RealVector* > tempStates;
  vector< RealVector* > tempGhostStates;
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    tempStates.push_back(m_cellStatesFlxPnt[iFlx]->getData());
    tempGhostStates.push_back(m_flxPntGhostSol[iFlx]->getData());
    //if (m_intCell->getID() == 191)
    //{
      CFLog(VERBOSE,"BndState " << iFlx << ": " << *(m_cellStatesFlxPnt[iFlx]->getData()) << "\n");
      CFLog(VERBOSE,"GstState " << iFlx << ": " << *(m_flxPntGhostSol[iFlx]->getData()) << "\n");
      //}
  }
  
  navierStokesVarSet->setGradientVars(tempStates,tempGradTerm,m_nbrFaceFlxPnts);
  navierStokesVarSet->setGradientVars(tempGhostStates,tempGhostGradTerm,m_nbrFaceFlxPnts);
  
  // Loop over solution pnts to calculate the grad updates
  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
  {
    // Loop over  variables
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      //set the grad updates to 0 
      m_gradUpdates[iSolPnt][iEq] = 0.0;
      
      // compute the face corrections to the gradients
      for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
      {
        const CFreal avgSol = (tempGradTerm(iEq,iFlx)+tempGhostGradTerm(iEq,iFlx))/2.0;
	const RealVector projectedCorr = (avgSol-tempGradTerm(iEq,iFlx))*m_faceJacobVecSizeFlxPnts[iFlx]*m_unitNormalFlxPnts[iFlx];
	if (m_intCell->getID() == 191 && iEq == 1)
    {
      CFLog(VERBOSE,"projecteddCorr Bnd of flx pnt " << (*m_faceFlxPntConn)[m_orient][iFlx] << " for state " << iSolPnt << ": " << projectedCorr << "\n");
      CFLog(VERBOSE,"deriv h: " << m_corrFctDiv[iSolPnt][(*m_faceFlxPntConn)[m_orient][iFlx]] << "\n");
    }
	/// @todo Check if this is also OK for triangles!!
	m_gradUpdates[iSolPnt][iEq] += projectedCorr*m_corrFctDiv[iSolPnt][(*m_faceFlxPntConn)[m_orient][iFlx]];
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
  if (getMethodData().getUpdateVarStr() == "Cons" && getMethodData().hasArtificialViscosity())
  {
    // Loop over solution pnts to calculate the grad updates
    for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
    {
      // Loop over  variables
      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      {
        //set the grad updates to 0 
        m_gradUpdates[iSolPnt][iEq] = 0.0;
      
        // compute the face corrections to the gradients
        for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
        {
          const CFreal avgSol = ((*(m_cellStatesFlxPnt[iFlx]))[iEq]+(*(m_flxPntGhostSol[iFlx]))[iEq])/2.0;
	  const RealVector projectedCorr = (avgSol-(*(m_cellStatesFlxPnt[iFlx]))[iEq])*m_faceJacobVecSizeFlxPnts[iFlx]*m_unitNormalFlxPnts[iFlx];
	
	  /// @todo Check if this is also OK for triangles!!
	  m_gradUpdates[iSolPnt][iEq] += projectedCorr*m_corrFctDiv[iSolPnt][(*m_faceFlxPntConn)[m_orient][iFlx]];
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
      
        // compute the face corrections to the gradients
        for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
        {
	  const RealVector transformedState = static_cast<RealVector&>(*m_updateToSolutionVecTrans->transform(m_cellStatesFlxPnt[iFlx]));
	  const RealVector transformedGhostState = static_cast<RealVector&>(*m_updateToSolutionVecTrans->transform(m_flxPntGhostSol[iFlx]));
	  const CFreal avgSol = (transformedState[iEq]+transformedGhostState[iEq])/2.0;
	  const RealVector projectedCorr = (avgSol-transformedState[iEq])*m_faceJacobVecSizeFlxPnts[iFlx]*m_unitNormalFlxPnts[iFlx];
	
	  /// @todo Check if this is also OK for triangles!!
	  m_gradUpdates[iSolPnt][iEq] += projectedCorr*m_corrFctDiv[iSolPnt][(*m_faceFlxPntConn)[m_orient][iFlx]];
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
  
  m_updateToSolutionVecTrans = getMethodData().getUpdateToSolutionVecTrans();
  
  m_updateToSolutionVecTrans->setup(2);
  
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

