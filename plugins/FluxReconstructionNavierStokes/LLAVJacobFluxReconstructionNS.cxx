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

#include "FluxReconstructionNavierStokes/LLAVJacobFluxReconstructionNS.hh"
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

MethodCommandProvider< LLAVJacobFluxReconstructionNS,
		       FluxReconstructionSolverData,
		       FluxReconstructionModule >
LLAVJacobFluxReconstructionNSFluxReconstructionProvider("LLAVJacobNS");

//////////////////////////////////////////////////////////////////////////////
  
LLAVJacobFluxReconstructionNS::LLAVJacobFluxReconstructionNS(const std::string& name) :
  LLAVJacobFluxReconstruction(name)
  {
  }

//////////////////////////////////////////////////////////////////////////////

void LLAVJacobFluxReconstructionNS::configure ( Config::ConfigArgs& args )
{
  LLAVJacobFluxReconstruction::configure(args);
}  

//////////////////////////////////////////////////////////////////////////////

void LLAVJacobFluxReconstructionNS::setFaceData(CFuint faceID)
{
  LLAVJacobFluxReconstruction::setFaceData(faceID);
  
  // get the gradients datahandle
  DataHandle< vector< RealVector > > gradientsAV = socket_gradientsAV.getDataHandle();

  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
    {
      const CFuint stateID = (*(m_states[iSide]))[iState]->getLocalID();
      m_cellGrads[iSide][iState] = &gradientsAV[stateID];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void LLAVJacobFluxReconstructionNS::setCellData()
{
  LLAVJacobFluxReconstruction::setCellData();
  
  // get the gradients datahandle
  DataHandle< vector< RealVector > > gradientsAV = socket_gradientsAV.getDataHandle();

  // get the grads in the current cell
  for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
  {
    const CFuint stateID = (*(m_cellStates))[iState]->getLocalID();
    m_cellGrads[0][iState] = &gradientsAV[stateID];
  }
}

//////////////////////////////////////////////////////////////////////////////

void LLAVJacobFluxReconstructionNS::setFaceNeighbourGradients()
{
  // get the gradients datahandle
  DataHandle< vector< RealVector > > gradientsAV = socket_gradientsAV.getDataHandle();

  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    // get neighbour states of other faces
    const CFuint nbrOtherFaces = m_otherFaceLocalIdxs[iSide].size();
    for (CFuint iFace = 0; iFace < nbrOtherFaces; ++iFace)
    {
      // get local face index
      const CFuint faceIdx = m_otherFaceLocalIdxs[iSide][iFace];

      // get neigbouring states
      const CFuint nbrFaceNeighbours = (*m_isFaceOnBoundary[iSide])[faceIdx] ? 1 : 2;
      for (CFuint iSide2 = 0; iSide2 < nbrFaceNeighbours; ++iSide2)
      {
        // get number of states
        const CFuint nbrStates = m_faceNghbrStates[iSide][iFace][iSide2]->size();

        // resize m_faceNghbrGrads[iSide][iFace][iSide2]
        m_faceNghbrGrads[iSide][iFace][iSide2].resize(nbrStates);

        // set the gradients
        for (CFuint iState = 0; iState < nbrStates; ++iState)
        {
          const CFuint stateID = (*m_faceNghbrStates[iSide][iFace][iSide2])[iState]->getLocalID();
          m_faceNghbrGrads[iSide][iFace][iSide2][iState] = &gradientsAV[stateID];
        }
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void LLAVJacobFluxReconstructionNS::setup()
{
  CFAUTOTRACE;
  
  // setup parent class
  LLAVJacobFluxReconstruction::setup();
}

//////////////////////////////////////////////////////////////////////////////

void LLAVJacobFluxReconstructionNS::unsetup()
{
  CFAUTOTRACE;
  
  // unsetup parent class
  LLAVJacobFluxReconstruction::unsetup();
}


//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

