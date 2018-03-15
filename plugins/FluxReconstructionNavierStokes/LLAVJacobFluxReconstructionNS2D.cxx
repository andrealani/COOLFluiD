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

#include "NavierStokes/Euler2DVarSet.hh"
#include "NavierStokes/EulerTerm.hh"

#include "FluxReconstructionNavierStokes/LLAVJacobFluxReconstructionNS2D.hh"
#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"

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

MethodCommandProvider< LLAVJacobFluxReconstructionNS2D,
		       FluxReconstructionSolverData,
		       FluxReconstructionModule >
LLAVJacobFluxReconstructionNS2DFluxReconstructionProvider("LLAVJacobNS2D");

//////////////////////////////////////////////////////////////////////////////
  
LLAVJacobFluxReconstructionNS2D::LLAVJacobFluxReconstructionNS2D(const std::string& name) :
  LLAVJacobFluxReconstruction(name),
  m_eulerVarSet(CFNULL),
  m_pData()
  {
  }

//////////////////////////////////////////////////////////////////////////////

void LLAVJacobFluxReconstructionNS2D::configure ( Config::ConfigArgs& args )
{
  LLAVJacobFluxReconstruction::configure(args);
}  

//////////////////////////////////////////////////////////////////////////////

void LLAVJacobFluxReconstructionNS2D::computeBndGradTerms(RealMatrix& gradTerm, RealMatrix& ghostGradTerm)
{ 
  vector< RealVector* > tempStates;
  vector< RealVector* > tempGhostStates;
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    tempStates.push_back(m_cellStatesFlxPnt[0][iFlx]->getData());
    tempGhostStates.push_back(m_flxPntGhostSol[iFlx]->getData());
  }
  
  setNecGradientVars(tempStates,gradTerm,m_nbrFaceFlxPnts);
  setNecGradientVars(tempGhostStates,ghostGradTerm,m_nbrFaceFlxPnts);
}

//////////////////////////////////////////////////////////////////////////////

void LLAVJacobFluxReconstructionNS2D::computeCellGradTerm(RealMatrix& gradTerm)
{ 
  vector< RealVector* > tempStates;
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    tempStates.push_back((*m_cellStates)[iSol]->getData());
  }
  
  setNecGradientVars(tempStates,gradTerm,m_nbrSolPnts);
}

//////////////////////////////////////////////////////////////////////////////

void LLAVJacobFluxReconstructionNS2D::computeFaceGradTerms(RealMatrix& gradTermL, RealMatrix& gradTermR)
{
  vector< vector< RealVector* > > tempStates;
  tempStates.resize(2);
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    tempStates[LEFT].push_back(m_cellStatesFlxPnt[LEFT][iFlx]->getData());
    tempStates[RIGHT].push_back(m_cellStatesFlxPnt[RIGHT][iFlx]->getData());
  }
  
  setNecGradientVars(tempStates[LEFT],gradTermL,m_nbrFaceFlxPnts);
  setNecGradientVars(tempStates[RIGHT],gradTermR,m_nbrFaceFlxPnts);
}

//////////////////////////////////////////////////////////////////////////////

void LLAVJacobFluxReconstructionNS2D::setNecGradientVars(const std::vector<RealVector*>& states, RealMatrix& values, const CFuint stateSize)
{
  const CFreal R = m_eulerVarSet->getModel()->getR();
  const CFreal ovCv = (m_eulerVarSet->getModel()->getGamma() - 1.)/R;
  
  for (CFuint i = 0; i < stateSize; ++i) 
  {
    const RealVector& state = *states[i];
    const CFreal ovRho = 1./state[0]; 
    values(1,i) = state[1]*ovRho;
    values(2,i) = state[2]*ovRho;
    
    const CFreal V2 = values(1,i)*values(1,i) + values(2,i)*values(2,i);
    values(3,i) = (state[3]*ovRho - 0.5*V2)*ovCv;
    values(0,i) = R*state[0]*values(3,i);
  }
}

//////////////////////////////////////////////////////////////////////////////

void LLAVJacobFluxReconstructionNS2D::setup()
{
  CFAUTOTRACE;
  LLAVJacobFluxReconstruction::setup();
  
  // get Euler 2D varset
  m_eulerVarSet = getMethodData().getUpdateVar().d_castTo<Euler2DVarSet>();
}

//////////////////////////////////////////////////////////////////////////////

void LLAVJacobFluxReconstructionNS2D::unsetup()
{
  CFAUTOTRACE;
  
  LLAVJacobFluxReconstruction::unsetup();
}


//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

