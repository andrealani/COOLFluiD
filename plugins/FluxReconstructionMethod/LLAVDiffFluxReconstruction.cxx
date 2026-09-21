// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/MethodCommandProvider.hh"

#include "Framework/CFSide.hh"
#include "Framework/MeshData.hh"
#include "Framework/BaseTerm.hh"

#include "MathTools/MathFunctions.hh"

#include "FluxReconstructionMethod/LLAVDiffFluxReconstruction.hh"
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

MethodCommandProvider< LLAVDiffFluxReconstruction,
		       FluxReconstructionSolverData,
		       FluxReconstructionModule >
LLAVDiffFluxReconstructionFluxReconstructionProvider("LLAVDiff");

//////////////////////////////////////////////////////////////////////////////
  
LLAVDiffFluxReconstruction::LLAVDiffFluxReconstruction(const std::string& name) :
  LLAVFluxReconstruction(name),
  m_physFlxPntRiemannFlux(),
  m_physDivContFlx(),
  m_cellGradsPtrsBackUp(),
  m_pData()
  {
    addConfigOptionsTo(this);
  }
  
  
//////////////////////////////////////////////////////////////////////////////

void LLAVDiffFluxReconstruction::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

void LLAVDiffFluxReconstruction::configure ( Config::ConfigArgs& args )
{
  FluxReconstructionSolverCom::configure(args);
}  

//////////////////////////////////////////////////////////////////////////////

void LLAVDiffFluxReconstruction::execute()
{
  LLAVFluxReconstruction::execute();
}

//////////////////////////////////////////////////////////////////////////////

void LLAVDiffFluxReconstruction::computeInterfaceFlxCorrection()
{
  // physical common flux
  DiffRHSFluxReconstruction::computeInterfaceFlxCorrection();

  m_physFlxPntRiemannFlux = m_flxPntRiemannFlux;

  // artificial viscosity common flux
  LLAVFluxReconstruction::computeInterfaceFlxCorrection();

  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_flxPntRiemannFlux[iFlx] += m_physFlxPntRiemannFlux[iFlx];

    // compute FI in the mapped coord frame
    for (CFuint iSide = 0; iSide < 2; ++iSide)
    {
      m_cellFlx[iSide][iFlx] = m_flxPntRiemannFlux[iFlx]*m_faceJacobVecSizeFlxPnts[iFlx][iSide];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void LLAVDiffFluxReconstruction::computeWaveSpeedUpdates(vector< CFreal >& waveSpeedUpd)
{
  // compute the wave speed updates for the neighbouring cells
  cf_assert(waveSpeedUpd.size() == 2);
  CFreal visc = 1.0;
  
  // get the correct m_faceIntegrationCoefs depending on the face type (only applicable for Prism for now) @todo should be updated for hybrid grid
  if (m_dim>2)
  {
    // get face geo
    const CFGeoShape::Type geo = m_face->getShape(); 

    if (geo == CFGeoShape::TRIAG) // triag face
    {
      //(*m_faceIntegrationCoefs).resize(m_nbrFaceFlxPnts);
      (m_faceIntegrationCoefs) = &(*m_faceIntegrationCoefsPerType)[0];
    }
    else  // quad face
    {
      //(*m_faceIntegrationCoefs).resize(m_nbrFaceFlxPnts);
      (m_faceIntegrationCoefs) = &(*m_faceIntegrationCoefsPerType)[1];
    } 
  }

  if(m_addUpdCoeff)
  {
    for (CFuint iSide = 0; iSide < 2; ++iSide)
    {
      waveSpeedUpd[iSide] = 0.0;
      //for (CFuint iFlx = 0; iFlx < m_cellFlx[iSide].size(); ++iFlx)
      for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
      {
        const CFreal jacobXJacobXIntCoef = m_faceJacobVecAbsSizeFlxPnts[iFlx]*
                                         m_faceJacobVecAbsSizeFlxPnts[iFlx]*
                                         (*m_faceIntegrationCoefs)[iFlx]*
                                         m_cflConvDiffRatio;
        //const CFreal rho = (*(m_cellStatesFlxPnt[iSide][iFlx]))[0];
        const CFreal epsilon = 0.5*(m_epsilonLR[LEFT][iFlx]+m_epsilonLR[RIGHT][iFlx]);
        const CFreal viscCoef = computeViscCoef(m_cellStatesFlxPnt[iSide][iFlx]);
        visc = epsilon*viscCoef;
      
        // transform update states to physical data to calculate eigenvalues
        waveSpeedUpd[iSide] += visc*jacobXJacobXIntCoef/m_cellVolume[iSide];
      }
    } 
  }
}

//////////////////////////////////////////////////////////////////////////////

void LLAVDiffFluxReconstruction::computeDivDiscontFlx(vector< RealVector >& residuals)
{
  // physical volume term, its boundary flux is computed by the diffusive boundary command
  DiffRHSFluxReconstruction::computeDivDiscontFlx(residuals);

  m_physDivContFlx = residuals;

  // switch the gradients of the cell to the artificial viscosity gradients
  m_cellGradsPtrsBackUp = m_cellGrads[LEFT];

  DataHandle< vector< RealVector > > gradientsAV = socket_gradientsAV.getDataHandle();

  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    m_cellGrads[LEFT][iSol] = &gradientsAV[(*m_cellStates)[iSol]->getLocalID()];
  }

  // artificial viscosity volume term and boundary flux
  LLAVFluxReconstruction::computeDivDiscontFlx(residuals);

  m_cellGrads[LEFT] = m_cellGradsPtrsBackUp;

  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    residuals[iSol] += m_physDivContFlx[iSol];
  }
}

//////////////////////////////////////////////////////////////////////////////

void LLAVDiffFluxReconstruction::setup()
{
  CFAUTOTRACE;
  
  // setup parent class
  LLAVFluxReconstruction::setup();

  m_physFlxPntRiemannFlux.resize(m_nbrFaceFlxPnts);
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_physFlxPntRiemannFlux[iFlx].resize(m_nbrEqs);
  }

  m_physDivContFlx.resize(m_nbrSolPnts);
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    m_physDivContFlx[iSol].resize(m_nbrEqs);
  }

  m_cellGradsPtrsBackUp.resize(m_nbrSolPnts);
}

//////////////////////////////////////////////////////////////////////////////

void LLAVDiffFluxReconstruction::unsetup()
{
  CFAUTOTRACE;
  
  // unsetup parent class
  LLAVFluxReconstruction::unsetup();
}


//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

