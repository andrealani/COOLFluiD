// Copyright (C) 2022 KU Leuven CmPA, Belgium
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
  ConvBndCorrectionsRHSJacobFluxReconstruction(name)
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
  CFuint nbOrients = bndFacesStartIdxs[0].size()-1;

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
    nbOrients = bndFacesStartIdxs[iTR].size()-1;

    // loop over different orientations
    for (m_orient = 0; m_orient < nbOrients; ++m_orient)
    {
      CFLog(VERBOSE,"m_orient: " << m_orient << "\n");
      
      // Reset the value of m_nbrFaceFlxPnts in case it is not the same for all faces (Prism)
      m_nbrFaceFlxPnts=(*m_faceFlxPntConn)[m_orient].size();

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

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

