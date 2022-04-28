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

#include "FluxReconstructionMultiFluidMHD/DiffRHSFluxReconstructionMFMHD.hh"
#include "FluxReconstructionMultiFluidMHD/FluxReconstructionMultiFluidMHD.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"
#include "MultiFluidMHD/DiffMFMHDVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::MultiFluidMHD;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< DiffRHSFluxReconstructionMFMHD,
		       FluxReconstructionSolverData,
		       FluxReconstructionMultiFluidMHDModule >
diffRHSMFMHDFluxReconstructionProvider("DiffRHSMFMHD");
  
//////////////////////////////////////////////////////////////////////////////
  
DiffRHSFluxReconstructionMFMHD::DiffRHSFluxReconstructionMFMHD(const std::string& name) :
  DiffRHSFluxReconstruction(name)
{
}
  
//////////////////////////////////////////////////////////////////////////////

void DiffRHSFluxReconstructionMFMHD::computeWaveSpeedUpdates(vector< CFreal >& waveSpeedUpd)
{
  // compute the wave speed updates for the neighbouring cells
  cf_assert(waveSpeedUpd.size() == 2);
  CFreal visc = 1.0;

  SafePtr< DiffMFMHDVarSet > diffMFMHDVarSet = m_diffusiveVarSet.d_castTo< DiffMFMHDVarSet >();
  const CFreal dynVisc = diffMFMHDVarSet->getCurrDynViscosity();
  
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    waveSpeedUpd[iSide] = 0.0;
    for (CFuint iFlx = 0; iFlx < m_cellFlx[iSide].size(); ++iFlx)
    {
      const CFreal jacobXJacobXIntCoef = m_faceJacobVecAbsSizeFlxPnts[iFlx]*
                                 m_faceJacobVecAbsSizeFlxPnts[iFlx]*
                                   (*m_faceIntegrationCoefs)[iFlx]*
                                   m_cflConvDiffRatio;
      const CFreal rho = diffMFMHDVarSet->getDensity(*(m_cellStatesFlxPnt[iSide][iFlx]));
      visc = dynVisc/rho;
      
      // transform update states to physical data to calculate eigenvalues
      waveSpeedUpd[iSide] += visc*jacobXJacobXIntCoef/m_cellVolume[iSide];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSFluxReconstructionMFMHD::prepareFluxComputation()
{
  const bool isPerturb = this->getMethodData().isPerturb();
  const CFuint iPerturbVar = this->getMethodData().iPerturbVar();
  SafePtr< DiffMFMHDVarSet > diffMFMHDVarSet = m_diffusiveVarSet.d_castTo< DiffMFMHDVarSet >();
  diffMFMHDVarSet->setComposition(m_avgSol, isPerturb, iPerturbVar);
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

