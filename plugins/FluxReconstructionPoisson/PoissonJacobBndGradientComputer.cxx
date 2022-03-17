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

#include "FluxReconstructionPoisson/PoissonJacobBndGradientComputer.hh"
#include "FluxReconstructionPoisson/FluxReconstructionPoisson.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::Poisson;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< PoissonJacobBndGradientComputer,
		       FluxReconstructionSolverData,
		       FluxReconstructionPoissonModule >
PoissonJacobBndGradientComputerProvider("ConvBndCorrectionsRHSJacobPoisson");
  
//////////////////////////////////////////////////////////////////////////////
  
PoissonJacobBndGradientComputer::PoissonJacobBndGradientComputer(const std::string& name) :
  ConvBndCorrectionsRHSJacobFluxReconstruction(name),
  m_diffVarSetPoisson(CFNULL),
  m_tempGradTerm(),
  m_tempGradTermGhost(),
  m_tempStates(),
  m_tempStatesGhost()
{
}

//////////////////////////////////////////////////////////////////////////////

void PoissonJacobBndGradientComputer::computeWaveSpeedUpdates(CFreal& waveSpeedUpd)
{
  // reset the wave speed update
  waveSpeedUpd = 0.0;
}

//////////////////////////////////////////////////////////////////////////////

void PoissonJacobBndGradientComputer::executeOnTrs()
{
  CFAUTOTRACE;

  // get InnerCells TopologicalRegionSet
  SafePtr<TopologicalRegionSet> cellTrs = MeshDataStack::getActive()->getTrs("InnerCells");

  // get current bnd face TRS
  SafePtr<TopologicalRegionSet> faceTrs = getCurrentTRS();
  
  CFLog(VERBOSE,"ConvBndCorrectionRHSJacobFluxReconstruction::executeOnTRS: " << faceTrs->getName() << "\n");

  // get bndFacesStartIdxs from FluxReconstructionMethodData
  map< std::string , vector< vector< CFuint > > >&
    bndFacesStartIdxsPerTRS = getMethodData().getBndFacesStartIdxs();
  vector< vector< CFuint > > bndFacesStartIdxs = bndFacesStartIdxsPerTRS[faceTrs->getName()];

  // number of face orientations (should be the same for all TRs)
  cf_assert(bndFacesStartIdxs.size() != 0);
  const CFuint nbOrients = bndFacesStartIdxs[0].size()-1;

  // number of TRs
  const CFuint nbTRs = faceTrs->getNbTRs();
  cf_assert(bndFacesStartIdxs.size() == nbTRs);

  // get the geodata of the face builder and set the TRSs
  FaceToCellGEBuilder::GeoData& geoData = m_faceBuilder->getDataGE();
  geoData.cellsTRS = cellTrs;
  geoData.facesTRS = faceTrs;
  geoData.isBoundary = true;
  
  // boolean telling whether there is a diffusive term
  const bool hasDiffTerm = getMethodData().hasDiffTerm() || getMethodData().hasArtificialViscosity();
  m_bcStateComputer->preProcess();
  // loop over TRs
  for (CFuint iTR = 0; iTR < nbTRs; ++iTR)
  {
    // loop over different orientations
    for (m_orient = 0; m_orient < nbOrients; ++m_orient)
    {
      CFLog(VERBOSE,"m_orient: " << m_orient << "\n");
      
      // select the correct flx pnts on the face out of all cell flx pnts for the current orient
      for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
      {
        m_flxPntsLocalCoords[iFlx] = (*m_allCellFlxPnts)[(*m_faceFlxPntConn)[m_orient][iFlx]];
      }
      
      // start and stop index of the faces with this orientation
      const CFuint startFaceIdx = bndFacesStartIdxs[iTR][m_orient  ];
      const CFuint stopFaceIdx  = bndFacesStartIdxs[iTR][m_orient+1];

      // loop over faces with this orientation
      for (CFuint faceID = startFaceIdx; faceID < stopFaceIdx; ++faceID)
      {
        // build the face GeometricEntity
        geoData.idx = faceID;
        m_face = m_faceBuilder->buildGE();

        // get the neighbouring cell
        m_intCell = m_face->getNeighborGeo(0);
	
	// get the states in the neighbouring cell
        m_cellStates = m_intCell->getStates();
	
        CFLog(VERBOSE,"cellID: " << m_intCell->getID() << "\n");
	CFLog(VERBOSE,"coord state 0: " << (((*m_cellStates)[0])->getCoordinates()) << "\n");

        // if cell is parallel updatable or the gradients have to be computed, compute the necessary data
//        if ((*m_cellStates)[0]->isParUpdatable() || hasDiffTerm)
//        {  
	  // set the bnd face data
	  setBndFaceData(m_face->getID());//faceID
	  
	  // compute the perturbed states and ghost states in the flx pnts
          computeFlxPntStates();
//	}
	
//	// if the cell is parallel updatable, compute the flx correction
//	if ((*m_cellStates)[0]->isParUpdatable())
//	{
//	  // compute FI-FD
//          computeInterfaceFlxCorrection();
//	  
//          // compute the wave speed updates
//          computeWaveSpeedUpdates(m_waveSpeedUpd);
//      
//          // update the wave speeds
//          updateWaveSpeed();
//       
//	  // compute the correction -(FI-FD)divh of the bnd face for each sol pnt
//          computeCorrection(m_corrections);
//	  
//	  // update the rhs
//          updateRHS();
//	}
	  
	// if there is a diffusive term, compute the gradients
//        if (hasDiffTerm)
//        {
          computeGradientBndFaceCorrections();
//        }
        
//        const CFuint iter = SubSystemStatusStack::getActive()->getNbIter();
//    
//        const CFuint iterFreeze = getMethodData().getFreezeJacobIter();
//    
//        const CFuint interval = iter - iterFreeze;
//      
//        if (!getMethodData().freezeJacob() || iter < iterFreeze || interval % getMethodData().getFreezeJacobInterval() == 0)
//        {
//	
//	  // if the cell is parallel updatable, compute the contribution to the numerical jacobian
//	  if ((*m_cellStates)[0]->isParUpdatable())
//	  {
//	    // compute the convective boundary flux correction contribution to the jacobian
//	    computeJacobConvBndCorrection();
//          }
//        }
        
        // release the face
        m_faceBuilder->releaseGE();
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void PoissonJacobBndGradientComputer::computeGradientBndFaceCorrections()
{ 
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_tempStates[iFlx] = m_cellStatesFlxPnt[iFlx]->getData();
    m_tempStatesGhost[iFlx] = m_flxPntGhostSol[iFlx]->getData();
  }
  
  m_diffVarSetPoisson->setGradientVars(m_tempStates,m_tempGradTerm,m_nbrFaceFlxPnts);
  m_diffVarSetPoisson->setGradientVars(m_tempStatesGhost,m_tempGradTermGhost,m_nbrFaceFlxPnts);
  
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

void PoissonJacobBndGradientComputer::setup()
{
  CFAUTOTRACE;
  
  // setup parent class
  ConvBndCorrectionsRHSJacobFluxReconstruction::setup();
  
  m_diffVarSetPoisson = (getMethodData().getDiffusiveVar()).d_castTo< Physics::Poisson::PoissonDiffVarSet >();
  cf_assert(m_diffVarSetPoisson.isNotNull());
  
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

