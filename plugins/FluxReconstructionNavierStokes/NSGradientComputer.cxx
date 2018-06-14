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

#include "FluxReconstructionNavierStokes/NSGradientComputer.hh"
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

MethodCommandProvider< NSGradientComputer,
		       FluxReconstructionSolverData,
		       FluxReconstructionNavierStokesModule >
nsGradientComputerProvider("ConvRHSNS");
  
//////////////////////////////////////////////////////////////////////////////
  
NSGradientComputer::NSGradientComputer(const std::string& name) :
  ConvRHSFluxReconstruction(name),
  m_diffusiveVarSet(CFNULL)
{
  addConfigOptionsTo(this);
}
  
//////////////////////////////////////////////////////////////////////////////

void NSGradientComputer::computeGradients()
{
  // get the diffusive varset
  m_diffusiveVarSet = getMethodData().getDiffusiveVar();
  
  SafePtr< NavierStokesVarSet > navierStokesVarSet = m_diffusiveVarSet.d_castTo< NavierStokesVarSet >();
  
  RealMatrix tempGradTerm(m_nbrEqs,m_nbrSolPnts);
  vector< RealVector* > tempStates;
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    tempStates.push_back((*m_cellStates)[iSol]->getData());
    if (m_cell->getID() == 191)
    {
      CFLog(VERBOSE,"State " << iSol << ": " << *((*m_cellStates)[iSol]->getData()) << "\n");
    }
  }
  
  navierStokesVarSet->setGradientVars(tempStates,tempGradTerm,m_nbrSolPnts);
  
  // Loop over solution pnts to calculate the grad updates
  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
  {
    // Loop over  variables
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      //set the grad updates to 0 
      m_gradUpdates[0][iSolPnt][iEq] = 0.0;

      // Loop over gradient directions
      for (CFuint iDir = 0; iDir < m_dim; ++iDir)
      {
        // Loop over solution pnts to count factor of all sol pnt polys
        for (CFuint jSolPnt = 0; jSolPnt < m_nbrSolPnts; ++jSolPnt)
        {
	  const RealVector projectedState = tempGradTerm(iEq,jSolPnt) * m_cellFluxProjVects[iDir][jSolPnt];
	  
	  if (m_cell->getID() == 191)
    {
      CFLog(VERBOSE,"Projected State " << jSolPnt << " on dir " << iDir << ": " << projectedState << " of eq " << iEq << "\n");
    }
	  
          // compute the grad updates
          m_gradUpdates[0][iSolPnt][iEq] += (*m_solPolyDerivAtSolPnts)[iSolPnt][iDir][jSolPnt]*projectedState;
	}
      }
    }
  }
  
  // get the gradients
  DataHandle< vector< RealVector > > gradients = socket_gradients.getDataHandle();

  // get jacobian determinants at solution points
  const std::valarray<CFreal> jacobDet =
      m_cell->computeGeometricShapeFunctionJacobianDeterminant(*m_solPntsLocalCoords);

  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    // get state ID
    const CFuint solID = (*m_cellStates)[iSol]->getLocalID();

    // inverse Jacobian determinant
    const CFreal invJacobDet = 1.0/jacobDet[iSol];

    // update gradients
    for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
    {
      gradients[solID][iGrad] += m_gradUpdates[0][iSol][iGrad];
      gradients[solID][iGrad] *= invJacobDet;
      if(m_cell->getID() == 191)
      {
	CFLog(VERBOSE, "Vol gradient updates " << iGrad << " of  " << iSol << ": " << m_gradUpdates[0][iSol][iGrad] << "\n");
      }
    }
  }
  
  // if needed, compute the gradients for the artificial viscosity
  if ((getMethodData().getUpdateVarStr() == "Cons") && getMethodData().hasArtificialViscosity())
  {
    // Loop over solution pnts to calculate the grad updates
    for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
    {
      // Loop over  variables
      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      {
        //set the grad updates to 0 
        m_gradUpdates[0][iSolPnt][iEq] = 0.0;

        // Loop over gradient directions
        for (CFuint iDir = 0; iDir < m_dim; ++iDir)
        {
          // Loop over solution pnts to count factor of all sol pnt polys
          for (CFuint jSolPnt = 0; jSolPnt < m_nbrSolPnts; ++jSolPnt)
          {
	    const RealVector projectedState = (*((*m_cellStates)[jSolPnt]))[iEq] * m_cellFluxProjVects[iDir][jSolPnt];
	  
            // compute the grad updates
            m_gradUpdates[0][iSolPnt][iEq] += (*m_solPolyDerivAtSolPnts)[iSolPnt][iDir][jSolPnt]*projectedState;
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
      const CFreal invJacobDet = 1.0/jacobDet[iSol];

      // update gradientsAV
      for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
      {
        gradientsAV[solID][iGrad] += m_gradUpdates[0][iSol][iGrad];
        gradientsAV[solID][iGrad] *= invJacobDet;
      }
    }
  }
  else if ( (getMethodData().getUpdateVarStr() == "Puvt" || getMethodData().getUpdateVarStr() == "RhoivtTv") && getMethodData().hasArtificialViscosity())
  { 
    // Loop over solution pnts to calculate the grad updates
    for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
    {
      // Loop over  variables
      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      {
        //set the grad updates to 0 
        m_gradUpdates[0][iSolPnt][iEq] = 0.0;

        // Loop over gradient directions
        for (CFuint iDir = 0; iDir < m_dim; ++iDir)
        {
          // Loop over solution pnts to count factor of all sol pnt polys
          for (CFuint jSolPnt = 0; jSolPnt < m_nbrSolPnts; ++jSolPnt)
          {
	    const RealVector transformedState = static_cast<RealVector&>(*m_updateToSolutionVecTrans->transform((*m_cellStates)[jSolPnt]));
	    
	    const RealVector projectedState = transformedState[iEq] * m_cellFluxProjVects[iDir][jSolPnt];
	  
            // compute the grad updates
            m_gradUpdates[0][iSolPnt][iEq] += (*m_solPolyDerivAtSolPnts)[iSolPnt][iDir][jSolPnt]*projectedState;
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
      const CFreal invJacobDet = 1.0/jacobDet[iSol];

      // update gradientsAV
      for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
      {
        gradientsAV[solID][iGrad] += m_gradUpdates[0][iSol][iGrad];
        gradientsAV[solID][iGrad] *= invJacobDet;
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void NSGradientComputer::computeGradientFaceCorrections()
{
  //ConvRHSFluxReconstruction::computeGradientFaceCorrections();
  // get the diffusive varset
  m_diffusiveVarSet = getMethodData().getDiffusiveVar();

  SafePtr< NavierStokesVarSet > navierStokesVarSet = m_diffusiveVarSet.d_castTo< NavierStokesVarSet >();

  RealMatrix tempGradTermL(m_nbrEqs,m_nbrFaceFlxPnts);
  RealMatrix tempGradTermR(m_nbrEqs,m_nbrFaceFlxPnts);
  vector< vector< RealVector* > > tempStates;
  tempStates.resize(2);
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    tempStates[LEFT].push_back(m_cellStatesFlxPnt[LEFT][iFlx]->getData());
    tempStates[RIGHT].push_back(m_cellStatesFlxPnt[RIGHT][iFlx]->getData());
  }
  
  navierStokesVarSet->setGradientVars(tempStates[LEFT],tempGradTermL,m_nbrFaceFlxPnts);
  navierStokesVarSet->setGradientVars(tempStates[RIGHT],tempGradTermR,m_nbrFaceFlxPnts);
  
  // Loop over solution pnts to calculate the grad updates
  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
  {
    // Loop over  variables
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      //set the grad updates to 0 
      m_gradUpdates[LEFT][iSolPnt][iEq] = 0.0;
      m_gradUpdates[RIGHT][iSolPnt][iEq] = 0.0;
      
      // compute the face corrections to the gradients
      for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
      {
        const CFreal avgSol = (tempGradTermL(iEq,iFlx)+tempGradTermR(iEq,iFlx))/2.0;
	const RealVector projectedCorrL = (avgSol-tempGradTermL(iEq,iFlx))*m_faceJacobVecSizeFlxPnts[iFlx][LEFT]*m_unitNormalFlxPnts[iFlx];
	const RealVector projectedCorrR = (avgSol-tempGradTermR(iEq,iFlx))*m_faceJacobVecSizeFlxPnts[iFlx][RIGHT]*m_unitNormalFlxPnts[iFlx];
	/// @todo Check if this is also OK for triangles!!
	m_gradUpdates[LEFT][iSolPnt][iEq] += projectedCorrL*m_corrFctDiv[iSolPnt][(*m_faceFlxPntConnPerOrient)[m_orient][LEFT][iFlx]];
	m_gradUpdates[RIGHT][iSolPnt][iEq] += projectedCorrR*m_corrFctDiv[iSolPnt][(*m_faceFlxPntConnPerOrient)[m_orient][RIGHT][iFlx]];
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
  
  // if needed, compute the gradients for the artificial viscosity
  if ((getMethodData().getUpdateVarStr() == "Cons") && getMethodData().hasArtificialViscosity())
  {
    // Loop over solution pnts to calculate the grad updates
    for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
    {
      // Loop over  variables
      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      {
        //set the grad updates to 0 
        m_gradUpdates[LEFT][iSolPnt][iEq] = 0.0;
        m_gradUpdates[RIGHT][iSolPnt][iEq] = 0.0;
      
        // compute the face corrections to the gradients
        for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
        {
          const CFreal avgSol = ((*(m_cellStatesFlxPnt[LEFT][iFlx]))[iEq]+(*(m_cellStatesFlxPnt[RIGHT][iFlx]))[iEq])/2.0;
	  const RealVector projectedCorrL = (avgSol-(*(m_cellStatesFlxPnt[LEFT][iFlx]))[iEq])*m_faceJacobVecSizeFlxPnts[iFlx][LEFT]*m_unitNormalFlxPnts[iFlx];
	  const RealVector projectedCorrR = (avgSol-(*(m_cellStatesFlxPnt[RIGHT][iFlx]))[iEq])*m_faceJacobVecSizeFlxPnts[iFlx][RIGHT]*m_unitNormalFlxPnts[iFlx];
	  /// @todo Check if this is also OK for triangles!!
	  m_gradUpdates[LEFT][iSolPnt][iEq] += projectedCorrL*m_corrFctDiv[iSolPnt][(*m_faceFlxPntConnPerOrient)[m_orient][LEFT][iFlx]];
	  m_gradUpdates[RIGHT][iSolPnt][iEq] += projectedCorrR*m_corrFctDiv[iSolPnt][(*m_faceFlxPntConnPerOrient)[m_orient][RIGHT][iFlx]];
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
  else if ((getMethodData().getUpdateVarStr() == "Puvt" || getMethodData().getUpdateVarStr() == "RhoivtTv") && getMethodData().hasArtificialViscosity())
  {
    // Loop over solution pnts to calculate the grad updates
    for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
    {
      // Loop over  variables
      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      {
        //set the grad updates to 0 
        m_gradUpdates[LEFT][iSolPnt][iEq] = 0.0;
        m_gradUpdates[RIGHT][iSolPnt][iEq] = 0.0;
      
        // compute the face corrections to the gradients
        for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
        {
	  const RealVector transformededStateL = static_cast<RealVector&>(*m_updateToSolutionVecTrans->transform(m_cellStatesFlxPnt[LEFT][iFlx]));
	  const RealVector transformededStateR = static_cast<RealVector&>(*m_updateToSolutionVecTrans->transform(m_cellStatesFlxPnt[RIGHT][iFlx]));
	  const CFreal avgSol = (transformededStateL[iEq]+transformededStateR[iEq])/2.0;
	  const RealVector projectedCorrL = (avgSol-transformededStateL[iEq])*m_faceJacobVecSizeFlxPnts[iFlx][LEFT]*m_unitNormalFlxPnts[iFlx];
	  const RealVector projectedCorrR = (avgSol-transformededStateR[iEq])*m_faceJacobVecSizeFlxPnts[iFlx][RIGHT]*m_unitNormalFlxPnts[iFlx];
	  /// @todo Check if this is also OK for triangles!!
	  m_gradUpdates[LEFT][iSolPnt][iEq] += projectedCorrL*m_corrFctDiv[iSolPnt][(*m_faceFlxPntConnPerOrient)[m_orient][LEFT][iFlx]];
	  m_gradUpdates[RIGHT][iSolPnt][iEq] += projectedCorrR*m_corrFctDiv[iSolPnt][(*m_faceFlxPntConnPerOrient)[m_orient][RIGHT][iFlx]];
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
}

//////////////////////////////////////////////////////////////////////////////

void NSGradientComputer::setup()
{
  CFAUTOTRACE;
  
  // setup parent class
  ConvRHSFluxReconstruction::setup();
  
  m_updateToSolutionVecTrans = getMethodData().getUpdateToSolutionVecTrans();
  
  m_updateToSolutionVecTrans->setup(2);
  
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

