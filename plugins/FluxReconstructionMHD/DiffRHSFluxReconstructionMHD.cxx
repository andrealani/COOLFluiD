// Copyright (C) 2022 KU Leuven CmPA, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/MethodCommandProvider.hh"

#include "Framework/CFSide.hh"
#include "Framework/MeshData.hh"
#include "Framework/BaseTerm.hh"

#include "MathTools/MathFunctions.hh"

#include "FluxReconstructionMHD/DiffRHSFluxReconstructionMHD.hh"
#include "FluxReconstructionMHD/FluxReconstructionMHD.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"
#include "MHD/MHDProjectionDiffVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::MHD;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< DiffRHSFluxReconstructionMHD,
		       FluxReconstructionSolverData,
		       FluxReconstructionMHDModule >
diffRHSMHDFluxReconstructionProvider("DiffRHSMHD");
  
//////////////////////////////////////////////////////////////////////////////
  
DiffRHSFluxReconstructionMHD::DiffRHSFluxReconstructionMHD(const std::string& name) :
  DiffRHSFluxReconstruction(name),
  m_tempGradTermL(),
  m_tempGradTermR(),
  m_diffusiveVarSet(CFNULL),
  m_tempStatesL(),
  m_tempStatesR()
{
  //addConfigOptionsTo(this);
}
  
//////////////////////////////////////////////////////////////////////////////

void DiffRHSFluxReconstructionMHD::computeWaveSpeedUpdates(vector< CFreal >& waveSpeedUpd)
{
  // compute the wave speed updates for the neighbouring cells
  cf_assert(waveSpeedUpd.size() == 2);
  CFreal visc = 1.0;

  const CFreal dynVisc = m_diffusiveVarSet->getCurrDynViscosity();
  
  const CFreal factorPr = 1.0;//min(m_diffusiveVarSet->getModel().getPrandtl(),1.0);
  cf_assert(factorPr>0.0);
  
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
      const CFreal rho = m_diffusiveVarSet->getDensity(*(m_cellStatesFlxPnt[iSide][iFlx]));
      visc = dynVisc/rho/factorPr;
      
      // transform update states to physical data to calculate eigenvalues
      waveSpeedUpd[iSide] += visc*jacobXJacobXIntCoef/m_cellVolume[iSide];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSFluxReconstructionMHD::prepareFluxComputation()
{
  const bool isPerturb = this->getMethodData().isPerturb();
  const CFuint iPerturbVar = this->getMethodData().iPerturbVar();

  m_diffusiveVarSet->setComposition(m_avgSol, isPerturb, iPerturbVar);
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSFluxReconstructionMHD::setup()
{
  CFAUTOTRACE;
  
  // setup parent class
  DiffRHSFluxReconstruction::setup();
  
  // get the diffusive varset
  m_diffusiveVarSet = (getMethodData().getDiffusiveVar()).d_castTo< MHDProjectionDiffVarSet >();
  
  m_tempGradTermL.resize(m_nbrEqs,m_nbrFaceFlxPnts);
  m_tempGradTermR.resize(m_nbrEqs,m_nbrFaceFlxPnts);
  
  m_tempStatesL.resize(m_nbrFaceFlxPnts);
  m_tempStatesR.resize(m_nbrFaceFlxPnts);
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

