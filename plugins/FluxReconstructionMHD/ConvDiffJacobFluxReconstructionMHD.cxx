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

#include "FluxReconstructionMHD/ConvDiffJacobFluxReconstructionMHD.hh"
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

MethodCommandProvider< ConvDiffJacobFluxReconstructionMHD,
		       FluxReconstructionSolverData,
		       FluxReconstructionMHDModule >
convDiffRHSJacobMHDFluxReconstructionProvider("ConvDiffRHSJacobMHD");
  
//////////////////////////////////////////////////////////////////////////////
  
ConvDiffJacobFluxReconstructionMHD::ConvDiffJacobFluxReconstructionMHD(const std::string& name) :
  ConvDiffJacobFluxReconstruction(name),
  m_diffusiveVarSetMHD(CFNULL),
  m_nbrSpecies(),
  m_pData(),
  m_pData2(),
  m_unpertGradVars(),
  m_pertGradVars()
{
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffJacobFluxReconstructionMHD::configure ( Config::ConfigArgs& args )
{
  ConvDiffJacobFluxReconstruction::configure(args);
} 
  
//////////////////////////////////////////////////////////////////////////////

void ConvDiffJacobFluxReconstructionMHD::computeWaveSpeedUpdates(vector< CFreal >& waveSpeedUpd)
{
  // compute the wave speed updates for the neighbouring cells
  cf_assert(waveSpeedUpd.size() == 2);
  
  // here convective and artificial parts are added!
  ConvDiffJacobFluxReconstruction::computeWaveSpeedUpdates(waveSpeedUpd);
          
  // now add diffusive part
  CFreal visc = 1.0;

  const CFreal dynVisc = m_diffusiveVarSetMHD->getCurrDynViscosity();
  
  const CFreal factorPr = 1.0;//min(m_diffusiveVarSetMHD->getModel().getPrandtl(),1.0);
  cf_assert(factorPr>0.0);
    
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    //for (CFuint iFlx = 0; iFlx < m_cellFlx[iSide].size(); ++iFlx)
    for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
    {
      const CFreal jacobXJacobXIntCoef = m_faceJacobVecAbsSizeFlxPnts[iFlx]*
                                 m_faceJacobVecAbsSizeFlxPnts[iFlx]*
                                   (*m_faceIntegrationCoefs)[iFlx]*
                                   m_cflConvDiffRatio;
      const CFreal rho = m_diffusiveVarSetMHD->getDensity(*(m_cellStatesFlxPnt[iSide][iFlx]));
      visc = dynVisc/rho/factorPr;
      
      // transform update states to physical data to calculate eigenvalues
      waveSpeedUpd[iSide] += visc*jacobXJacobXIntCoef/m_cellVolume[iSide];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffJacobFluxReconstructionMHD::prepareFluxComputation()
{
  const bool isPerturb = this->getMethodData().isPerturb();
  const CFuint iPerturbVar = this->getMethodData().iPerturbVar();

  m_diffusiveVarSetMHD->setComposition(m_avgSol, isPerturb, iPerturbVar);
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffJacobFluxReconstructionMHD::setup()
{
  CFAUTOTRACE;
  
  // setup parent class
  ConvDiffJacobFluxReconstruction::setup();
  
  // get the diffusive varset
  m_diffusiveVarSetMHD = (getMethodData().getDiffusiveVar()).d_castTo< MHDProjectionDiffVarSet >();
  
  m_updateToSolutionVecTrans = getMethodData().getUpdateToSolutionVecTrans();
  
  m_updateToSolutionVecTrans->setup(2);
  
  m_unpertGradVars.resize(m_nbrEqs);
  m_pertGradVars.resize(m_nbrEqs);
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD
