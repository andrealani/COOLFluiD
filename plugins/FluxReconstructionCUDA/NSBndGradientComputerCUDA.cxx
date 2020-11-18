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

#include "FluxReconstructionCUDA/NSBndGradientComputerCUDA.hh"
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

MethodCommandProvider< NSBndGradientComputerCUDA,
		       FluxReconstructionSolverData,
		       FluxReconstructionNavierStokesModule >
NsBndGradientComputerCUDAProvider("ConvBndCorrectionsRHSNSCUDA");
  
//////////////////////////////////////////////////////////////////////////////
  
NSBndGradientComputerCUDA::NSBndGradientComputerCUDA(const std::string& name) :
  ConvBndCorrectionsRHSFluxReconstruction(name),
  socket_gradientsCUDA("gradientsCUDA"),
  socket_gradientsAVCUDA("gradientsAVCUDA"),
  m_diffusiveVarSet(CFNULL),
  m_tempGradTerm(),
  m_tempGradTermGhost(),
  m_tempStates(),
  m_tempStatesGhost()
{
}

//////////////////////////////////////////////////////////////////////////////

void NSBndGradientComputerCUDA::computeGradientBndFaceCorrections()
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
  DataHandle< CFreal > gradients = socket_gradientsCUDA.getDataHandle();

  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    // get state ID
    const CFuint solID = (*m_cellStates)[iSol]->getLocalID();

    // update gradients
    for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
    {
      for (CFuint iDim = 0; iDim < m_dim; ++iDim)
      {
        gradients[solID*m_nbrEqs*m_dim+iGrad*m_dim+iDim] += m_gradUpdates[iSol][iGrad][iDim];
      } 
    }
  }
  
  // if needed, compute the gradients for the artificial viscosity
  if (getMethodData().hasArtificialViscosity())
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
    DataHandle< CFreal > gradientsAV = socket_gradientsAVCUDA.getDataHandle();

    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      // get state ID
      const CFuint solID = (*m_cellStates)[iSol]->getLocalID();

      // update gradients
      for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
      {
        for (CFuint iDim = 0; iDim < m_dim; ++iDim)
        {
          gradientsAV[solID*m_nbrEqs*m_dim+iGrad*m_dim+iDim] += m_gradUpdates[iSol][iGrad][iDim];
        }
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> >
NSBndGradientComputerCUDA::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result = ConvBndCorrectionsRHSFluxReconstruction::needsSockets();

  result.push_back(&socket_gradientsCUDA);
  result.push_back(&socket_gradientsAVCUDA);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void NSBndGradientComputerCUDA::setup()
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

