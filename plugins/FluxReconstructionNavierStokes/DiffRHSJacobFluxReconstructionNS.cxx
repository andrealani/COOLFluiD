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

#include "FluxReconstructionNavierStokes/DiffRHSJacobFluxReconstructionNS.hh"
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

MethodCommandProvider< DiffRHSJacobFluxReconstructionNS,
		       FluxReconstructionSolverData,
		       FluxReconstructionNavierStokesModule >
diffRHSJacobNSFluxReconstructionProvider("DiffRHSJacobNS");
  
//////////////////////////////////////////////////////////////////////////////
  
DiffRHSJacobFluxReconstructionNS::DiffRHSJacobFluxReconstructionNS(const std::string& name) :
  DiffRHSJacobFluxReconstruction(name)
{
  addConfigOptionsTo(this);
}
  
//////////////////////////////////////////////////////////////////////////////

void DiffRHSJacobFluxReconstructionNS::computeWaveSpeedUpdates(vector< CFreal >& waveSpeedUpd)
{
  // compute the wave speed updates for the neighbouring cells
  cf_assert(waveSpeedUpd.size() == 2);
  CFreal visc = 1.0;

  SafePtr< NavierStokesVarSet > navierStokesVarSet = m_diffusiveVarSet.d_castTo< NavierStokesVarSet >();
  const CFreal dynVisc = navierStokesVarSet->getCurrDynViscosity();
  
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    waveSpeedUpd[iSide] = 0.0;
    for (CFuint iFlx = 0; iFlx < m_cellFlx[iSide].size(); ++iFlx)
    {
      const CFreal jacobXJacobXIntCoef = m_faceJacobVecAbsSizeFlxPnts[iFlx]*
                                 m_faceJacobVecAbsSizeFlxPnts[iFlx]*
                                   (*m_faceIntegrationCoefs)[iFlx]*
                                   m_cflConvDiffRatio;
      const CFreal rho = navierStokesVarSet->getDensity(*(m_cellStatesFlxPnt[iSide][iFlx]));
      visc = dynVisc/rho;
      
      // transform update states to physical data to calculate eigenvalues
      waveSpeedUpd[iSide] += visc*jacobXJacobXIntCoef/m_cellVolume[iSide];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSJacobFluxReconstructionNS::computeBndGradTerms(RealMatrix& gradTerm, RealMatrix& ghostGradTerm)
{ 
  SafePtr< NavierStokesVarSet > navierStokesVarSet = m_diffusiveVarSet.d_castTo< NavierStokesVarSet >();

  vector< RealVector* > tempStates;
  vector< RealVector* > tempGhostStates;
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    tempStates.push_back(m_cellStatesFlxPnt[0][iFlx]->getData());
    tempGhostStates.push_back(m_flxPntGhostSol[iFlx]->getData());
  }
  
  navierStokesVarSet->setGradientVars(tempStates,gradTerm,m_nbrFaceFlxPnts);
  navierStokesVarSet->setGradientVars(tempGhostStates,ghostGradTerm,m_nbrFaceFlxPnts);
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSJacobFluxReconstructionNS::computeCellGradTerm(RealMatrix& gradTerm)
{ 
  SafePtr< NavierStokesVarSet > navierStokesVarSet = m_diffusiveVarSet.d_castTo< NavierStokesVarSet >();
  
  vector< RealVector* > tempStates;
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    tempStates.push_back((*m_cellStates)[iSol]->getData());
  }
  
  navierStokesVarSet->setGradientVars(tempStates,gradTerm,m_nbrSolPnts);
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSJacobFluxReconstructionNS::computeFaceGradTerms(RealMatrix& gradTermL, RealMatrix& gradTermR)
{
  SafePtr< NavierStokesVarSet > navierStokesVarSet = m_diffusiveVarSet.d_castTo< NavierStokesVarSet >();

  vector< vector< RealVector* > > tempStates;
  tempStates.resize(2);
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    tempStates[LEFT].push_back(m_cellStatesFlxPnt[LEFT][iFlx]->getData());
    tempStates[RIGHT].push_back(m_cellStatesFlxPnt[RIGHT][iFlx]->getData());
  }
  
  navierStokesVarSet->setGradientVars(tempStates[LEFT],gradTermL,m_nbrFaceFlxPnts);
  navierStokesVarSet->setGradientVars(tempStates[RIGHT],gradTermR,m_nbrFaceFlxPnts);
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

