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

#include "NavierStokes/NavierStokesVarSet.hh"

#include "FluxReconstructionNavierStokes/LLAVFluxReconstructionNS.hh"
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

MethodCommandProvider< LLAVFluxReconstructionNS,
		       FluxReconstructionSolverData,
		       FluxReconstructionModule >
LLAVFluxReconstructionNSFluxReconstructionProvider("LLAVNS");

//////////////////////////////////////////////////////////////////////////////
  
LLAVFluxReconstructionNS::LLAVFluxReconstructionNS(const std::string& name) :
  LLAVFluxReconstruction(name),
  m_diffVarSet(CFNULL),
  m_pData(),
  m_gradsBackUp()
  {
  }

//////////////////////////////////////////////////////////////////////////////

void LLAVFluxReconstructionNS::configure ( Config::ConfigArgs& args )
{
  LLAVFluxReconstruction::configure(args);
}  

//////////////////////////////////////////////////////////////////////////////

void LLAVFluxReconstructionNS::setFaceData(CFuint faceID)
{
  LLAVFluxReconstruction::setFaceData(faceID);

  m_diffVarSet = getMethodData().getDiffusiveVar();

  SafePtr< NavierStokesVarSet > navierStokesVarSet = m_diffVarSet.d_castTo< NavierStokesVarSet >();

  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    vector< vector< RealVector* > > grads;
    vector< RealVector* > tempStates;
    vector< vector< RealVector > > temp;
    temp.resize(m_nbrSolPnts);
    grads.resize(m_nbrSolPnts);

    // make a back up of the grads
    for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
    {
      tempStates.push_back((*(m_states[iSide]))[iState]->getData());

      temp[iState] = *(m_cellGrads[iSide][iState]);
      
      grads[iState].resize(m_nbrEqs);

      for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
      {
        cf_assert(temp[iState].size() == m_nbrEqs);
        grads[iState][iVar] = & (temp[iState][iVar]);
      }
    }

    navierStokesVarSet->setStateGradients(tempStates,grads,m_gradsBackUp[iSide],m_nbrSolPnts);

    // make store the correct gradients
    for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
    {    
      for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
      { 
        (*m_cellGrads[iSide][iState])[iVar] = *(m_gradsBackUp[iSide][iState][iVar]);
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void LLAVFluxReconstructionNS::setCellData()
{
  LLAVFluxReconstruction::setCellData();
  
  m_diffVarSet = getMethodData().getDiffusiveVar();
  
  SafePtr< NavierStokesVarSet > navierStokesVarSet = m_diffVarSet.d_castTo< NavierStokesVarSet >();
  
  vector< vector< RealVector* > > grads;
  vector< RealVector* > tempStates;
  vector< vector< RealVector > > temp;
  temp.resize(m_nbrSolPnts);
  grads.resize(m_nbrSolPnts);
  
  // make a back up of the grads
  for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
  {
    tempStates.push_back((*m_cellStates)[iState]->getData());
    
    temp[iState] = *(m_cellGrads[0][iState]);
      
    grads[iState].resize(m_nbrEqs);
    
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      cf_assert(temp[iState].size() == m_nbrEqs);
      grads[iState][iVar] = & (temp[iState][iVar]);
    }
  }
  
  navierStokesVarSet->setStateGradients(tempStates,grads,m_gradsBackUp[0],m_nbrSolPnts);
    
  // make store the correct gradients
  for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
  {    
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    { 
      (*m_cellGrads[0][iState])[iVar] = *(m_gradsBackUp[0][iState][iVar]);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void LLAVFluxReconstructionNS::setup()
{
  CFAUTOTRACE;
  
  LLAVFluxReconstruction::setup();
  
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
}

//////////////////////////////////////////////////////////////////////////////

void LLAVFluxReconstructionNS::unsetup()
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
  m_gradsBackUp[LEFT].clear();
  m_gradsBackUp[RIGHT].clear();
  m_gradsBackUp.clear();
  
  LLAVFluxReconstruction::unsetup();
}


//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

