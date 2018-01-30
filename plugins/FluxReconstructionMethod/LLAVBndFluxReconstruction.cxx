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

#include "FluxReconstructionMethod/LLAVBndFluxReconstruction.hh"
#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {
    
//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< LLAVBndFluxReconstruction,
		       FluxReconstructionSolverData,
		       FluxReconstructionModule >
LLAVBndFluxReconstructionFluxReconstructionProvider("LLAVBndRHS");

//////////////////////////////////////////////////////////////////////////////
  
LLAVBndFluxReconstruction::LLAVBndFluxReconstruction(const std::string& name) :
  DiffBndCorrectionsRHSFluxReconstruction(name),
  m_updateVarSet(CFNULL),
  m_order(),
  m_transformationMatrix(),
  m_statesPMinOne(),
  m_epsilon(),
  m_s0(),
  m_s(),
  m_epsilon0(),
  m_kappa()
  {
    addConfigOptionsTo(this);
  }
  
  
//////////////////////////////////////////////////////////////////////////////

void LLAVBndFluxReconstruction::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

void LLAVBndFluxReconstruction::configure ( Config::ConfigArgs& args )
{
  FluxReconstructionSolverCom::configure(args);
}  

//////////////////////////////////////////////////////////////////////////////

void LLAVBndFluxReconstruction::executeOnTrs()
{
  CFAUTOTRACE;
  
  CFLog(VERBOSE, "LLAVBndFluxReconstruction::executeOnTrs()\n");
  
  // get InnerCells TopologicalRegionSet
  SafePtr<TopologicalRegionSet> cellTrs = MeshDataStack::getActive()->getTrs("InnerCells");

  // get current QuadFreeBCFluxReconstruction TRS
  SafePtr<TopologicalRegionSet> faceTrs = getCurrentTRS();

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

        // if cell is parallel updatable or the gradients need to be computed, compute the correction flux
        if ((*m_cellStates)[0]->isParUpdatable())
        {
	  // set the face ID in the BCStateComputer
	  m_bcStateComputer->setFaceID(m_face->getID());
	  
	  // set the bnd face data
	  setBndFaceData(m_face->getID());//faceID 
	  
	  // compute the states and ghost states in the flx pnts
	  computeFlxPntStates();

	  // compute the interface flux
          computeInterfaceFlxCorrection();

          // compute the wave speed updates
          //computeWaveSpeedUpdates(m_waveSpeedUpd);

          // update the wave speeds
          //updateWaveSpeed();

	  // compute the correction -(FI)divh of the bnd face for each sol pnt
          computeCorrection(m_corrections);

	  // update the rhs
          updateRHS();
        } 
        
        // release the face
        m_faceBuilder->releaseGE();
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void LLAVBndFluxReconstruction::computeInterfaceFlxCorrection()
{ 
  // compute the riemann flux in the flx pnts
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      *(m_avgGrad[iVar]) = (*(m_flxPntGhostGrads[iFlxPnt][iVar]) + *(m_cellGradFlxPnt[iFlxPnt][iVar]))/2.0;
      //*(m_avgGrad[iVar]) = *(m_cellGradFlxPnt[iFlxPnt][iVar]);
    }

    // compute FI
    m_flxPntRiemannFlux[iFlxPnt] = 0.0;
    
    for (CFuint iDim = 0; iDim < m_dim; ++iDim)
    {
      for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
      {
        m_flxPntRiemannFlux[iFlxPnt][iVar] += m_epsilon*((*(m_avgGrad[iVar]))[iDim])*m_unitNormalFlxPnts[iFlxPnt][iDim];
      }
    }
    
    // compute FI in the local frame
    m_cellFlx[iFlxPnt] = (m_flxPntRiemannFlux[iFlxPnt])*m_faceJacobVecSizeFlxPnts[iFlxPnt];
  }
}

//////////////////////////////////////////////////////////////////////////////

void LLAVBndFluxReconstruction::setBndFaceData(CFuint faceID)
{
  DiffBndCorrectionsRHSFluxReconstruction::setBndFaceData(faceID);
  
  computeEpsilon();
}

//////////////////////////////////////////////////////////////////////////////

void LLAVBndFluxReconstruction::computeProjStates(std::vector< RealVector >& projStates)
{
  for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
  {
    RealVector temp(projStates.size());
    
    for (CFuint iSol = 0; iSol < projStates.size(); ++iSol)
    {
      temp[iSol] = (*((*m_cellStates)[iSol]))[iEq];
    }

    RealVector tempProj(projStates.size());

    tempProj = m_transformationMatrix*temp;

    for (CFuint iSol = 0; iSol < projStates.size(); ++iSol)
    {
      projStates[iSol][iEq] = tempProj[iSol];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void LLAVBndFluxReconstruction::computeEpsilon()
{
  computeEpsilon0();
  
  computeSmoothness();
  
  m_kappa = 5.0;
  
  if (m_s < m_s0 - m_kappa)
  {
    m_epsilon = 0.0;
  }
  else if (m_s > m_s0 + m_kappa)
  {
    m_epsilon = m_epsilon0;
  }
  else
  {
    m_epsilon = m_epsilon0*0.5*(1.0 + sin(0.5*MathTools::MathConsts::CFrealPi()*(m_s-m_s0)/m_kappa));
  }
}

//////////////////////////////////////////////////////////////////////////////

void LLAVBndFluxReconstruction::computeEpsilon0()
{ 
  // get the datahandle of the update coefficients
  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
  
  const CFreal wavespeed = updateCoeff[(*m_cellStates)[0]->getLocalID()];
  
  const CFreal deltaKsi = 1.0/(m_order+2.0);
  
  const CFreal Pe = 4.0;//2.0 P1: 20
  
  m_epsilon0 = wavespeed*(2.0/Pe - deltaKsi/Pe);

}

//////////////////////////////////////////////////////////////////////////////

void LLAVBndFluxReconstruction::computeSmoothness()
{ 
  CFreal sNum = 0.0;
  
  CFreal sDenom = 0.0;
  
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    CFreal stateP = (*((*m_cellStates)[iSol]))[0];
    CFreal diffStatesPPMinOne = stateP - m_statesPMinOne[iSol][0];
    sNum += diffStatesPPMinOne*diffStatesPPMinOne;
    sDenom += stateP*stateP;
  }
  if (sNum <= MathTools::MathConsts::CFrealEps() || sDenom <= MathTools::MathConsts::CFrealEps())
  {
    m_s = -100.0;
  }
  else
  {
    m_s = log10(sNum/sDenom);
  }
}

//////////////////////////////////////////////////////////////////////////////

void LLAVBndFluxReconstruction::setup()
{
  CFAUTOTRACE;
  DiffBndCorrectionsRHSFluxReconstruction::setup();

  // get the update varset
  m_updateVarSet = getMethodData().getUpdateVar();
  
  // get the local spectral FD data
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
  cf_assert(frLocalData.size() > 0);
  // for now, there should be only one type of element
  cf_assert(frLocalData.size() == 1);
  
  const CFPolyOrder::Type order = frLocalData[0]->getPolyOrder();
  
  m_order = static_cast<CFuint>(order);

  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    RealVector temp(m_nbrEqs);
    temp = 0.0;
    m_statesPMinOne.push_back(temp);
  }
  
  SafePtr<RealMatrix> vdm = frLocalData[0]->getVandermondeMatrix();
  
  SafePtr<RealMatrix> vdmInv = frLocalData[0]->getVandermondeMatrixInv();
  
  RealMatrix temp(m_nbrSolPnts,m_nbrSolPnts);
  temp = 0.0;
  for (CFuint idx = 0; idx < (m_order)*(m_order); ++idx)
  {
    temp(idx,idx) = 1.0;
  }
  
  m_transformationMatrix.resize(m_nbrSolPnts,m_nbrSolPnts);
  
  m_transformationMatrix = (*vdm)*temp*(*vdmInv);
  
  m_s0 = -3.0*log10(static_cast<CFreal>(m_order));
}

//////////////////////////////////////////////////////////////////////////////

void LLAVBndFluxReconstruction::unsetup()
{
  CFAUTOTRACE;
  
  DiffBndCorrectionsRHSFluxReconstruction::unsetup();
}


//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

