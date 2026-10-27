// Copyright (C) 2026 KU Leuven CmPA, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/MethodCommandProvider.hh"

#include "FluxReconstructionMHD/MHDJacobBlendingCommands.hh"
#include "FluxReconstructionMHD/FluxReconstructionMHD.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< MHDJacobGradientComputerBlending,
                       FluxReconstructionSolverData,
                       FluxReconstructionMHDModule >
MHDJacobGradientComputerBlendingProvider("ConvRHSJacobMHDBlending");

MethodCommandProvider< MHDJacobBndGradientComputerBlending,
                       FluxReconstructionSolverData,
                       FluxReconstructionMHDModule >
MHDJacobBndGradientComputerBlendingProvider("ConvBndCorrectionsRHSJacobMHDBlending");

//////////////////////////////////////////////////////////////////////////////

MHDJacobGradientComputerBlending::MHDJacobGradientComputerBlending(const std::string& name) :
  ConvRHSJacobFluxReconstructionBlending(name)
{
}

//////////////////////////////////////////////////////////////////////////////

MHDJacobBndGradientComputerBlending::MHDJacobBndGradientComputerBlending(const std::string& name) :
  ConvBndCorrectionsRHSJacobFluxReconstructionBlending(name)
{
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD
