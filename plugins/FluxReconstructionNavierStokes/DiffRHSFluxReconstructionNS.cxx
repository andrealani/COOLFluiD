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

#include "FluxReconstructionNavierStokes/DiffRHSFluxReconstructionNS.hh"
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

MethodCommandProvider< DiffRHSFluxReconstructionNS,
		       FluxReconstructionSolverData,
		       FluxReconstructionNavierStokesModule >
diffRHSFluxReconstructionProvider("DiffRHS");
  
//////////////////////////////////////////////////////////////////////////////
  
DiffRHSFluxReconstructionNS::DiffRHSFluxReconstructionNS(const std::string& name) :
  DiffRHSFluxReconstruction(name)
{
  addConfigOptionsTo(this);
}
  
//////////////////////////////////////////////////////////////////////////////

void DiffRHSFluxReconstructionNS::computeWaveSpeedUpdates(vector< CFreal >& waveSpeedUpd)
{
  // compute the wave speed updates for the neighbouring cells
  cf_assert(waveSpeedUpd.size() == 2);
  CFreal visc = 1.0;
  /// @todo needs to be changed for non-NS
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

void DiffRHSFluxReconstructionNS::addPartitionFacesCorrection()
{
  // get InnerCells TopologicalRegionSet
  SafePtr<TopologicalRegionSet> cellTrs = MeshDataStack::getActive()->getTrs("InnerCells");

  // get current QuadFreeBCFluxReconstruction TRS
  SafePtr<TopologicalRegionSet> faceTrs = MeshDataStack::getActive()->getTrs("PartitionFaces");

  // get the partition face start indexes
  vector< CFuint >& partitionFacesStartIdxs = getMethodData().getPartitionFacesStartIdxs();

  // get number of face orientations
  const CFuint nbrFaceOrients = partitionFacesStartIdxs.size()-1;

  // get the geodata of the face builder and set the TRSs
  FaceToCellGEBuilder::GeoData& geoData = m_faceBuilder->getDataGE();
  geoData.cellsTRS = cellTrs;
  geoData.facesTRS = faceTrs;
  geoData.isBoundary = true;
  
  // get the local FR data
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
  
  // get solution polynomial values at nodes
  vector< vector< CFreal > > solPolyValsAtNodes
	= frLocalData[0]->getSolPolyValsAtNode(*m_flxPntsLocalCoords);
  
  // loop over different orientations
  for (m_orient = 0; m_orient < nbrFaceOrients; ++m_orient)
  {
    CFLog(VERBOSE, "Partition Orient = " << m_orient << "\n");
    // start and stop index of the faces with this orientation
    const CFuint faceStartIdx = partitionFacesStartIdxs[m_orient  ];
    const CFuint faceStopIdx  = partitionFacesStartIdxs[m_orient+1];

    // loop over faces with this orientation
    for (CFuint faceID = faceStartIdx; faceID < faceStopIdx; ++faceID)
    {
      // build the face GeometricEntity
      geoData.idx = faceID;
      m_face = m_faceBuilder->buildGE();
      
      // get the neighbouring cells
      m_cells[0] = m_face->getNeighborGeo(0);

      // get the states in the neighbouring cells
      m_states[0] = m_cells[0]->getStates();
      
      // compute volume
      m_cellVolume[0] = m_cells[0]->computeVolume();
      
      // compute flux point coordinates
      SafePtr< vector<RealVector> > flxLocalCoords = frLocalData[0]->getFaceFlxPntsFaceLocalCoords();
  
      // Loop over flux points to extrapolate the states to the flux points and get the 
      // discontinuous normal flux in the flux points and the Riemann flux
      for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
      {
        // get face Jacobian vector sizes in the flux points
        DataHandle< vector< CFreal > >
        faceJacobVecSizeFaceFlxPnts = socket_faceJacobVecSizeFaceFlxPnts.getDataHandle();
  
        // get face Jacobian vector size
        m_faceJacobVecAbsSizeFlxPnts[iFlxPnt] = faceJacobVecSizeFaceFlxPnts[m_face->getID()][iFlxPnt];//faceID
      }
      
      // Loop over flux points to extrapolate the states to the flux points and get the 
      // discontinuous normal flux in the flux points and the Riemann flux
      for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
      {
        State* tempVector = new State(solPolyValsAtNodes[iFlxPnt][0]*(*((*(m_states[0]))[0]->getData())));

        m_cellStatesFlxPnt[0][iFlxPnt] = tempVector;
  
        for (CFuint iSol = 1; iSol < m_nbrSolPnts; ++iSol)
        {
          CFLog(DEBUG_MIN,"cellStates: " << *((*m_cellStates)[iSol]->getData()) << "\n");
          *(m_cellStatesFlxPnt[0][iFlxPnt]) = (State) (*(m_cellStatesFlxPnt[0][iFlxPnt]->getData()) + 
			  solPolyValsAtNodes[iFlxPnt][iSol]*(*((*(m_states[0]))[iSol]->getData())));
        }
      }
           
      CFreal waveSpeedUpd = 0.0;
      CFreal visc = 1.0;
      /// @todo needs to be changed for non-NS
      SafePtr< NavierStokesVarSet > navierStokesVarSet = m_diffusiveVarSet.d_castTo< NavierStokesVarSet >();
      const CFreal dynVisc = navierStokesVarSet->getCurrDynViscosity();
      
      for (CFuint iFlx = 0; iFlx < m_cellFlx.size(); ++iFlx)
      {  
	const CFreal jacobXJacobXIntCoef = m_faceJacobVecAbsSizeFlxPnts[iFlx]*
                                           m_faceJacobVecAbsSizeFlxPnts[iFlx]*
                                           (*m_faceIntegrationCoefs)[iFlx]*
                                           m_cflConvDiffRatio;
        const CFreal rho = navierStokesVarSet->getDensity(*(m_cellStatesFlxPnt[0][iFlx]));
        visc = dynVisc/rho;
      
        // transform update states to physical data to calculate eigenvalues
        waveSpeedUpd += visc*jacobXJacobXIntCoef/m_cellVolume[0];
				   
      }

      // get the datahandle of the update coefficients
      DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();

      for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
      {
        const CFuint solID = (*(m_states[0]))[iSol]->getLocalID();
        updateCoeff[solID] += waveSpeedUpd;
      }

      // release the GeometricEntity
      m_faceBuilder->releaseGE();

    }
  }
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

