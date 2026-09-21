// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/MethodCommandProvider.hh"

#include "FluxReconstructionNavierStokes/NSJacobBndGradientComputerSubcellBlending.hh"
#include "FluxReconstructionNavierStokes/FluxReconstructionNavierStokes.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< NSJacobBndGradientComputerSubcellBlending,
		       FluxReconstructionSolverData,
		       FluxReconstructionNavierStokesModule >
NSJacobBndGradientComputerSubcellBlendingProvider("ConvBndCorrectionsRHSJacobNSSubcellBlending");

MethodCommandProvider< NSJacobBndGradientComputerSubcellBlending,
		       FluxReconstructionSolverData,
		       FluxReconstructionNavierStokesModule >
NSJacobBndGradientComputerSubcellBlendingTurbProvider("ConvBndCorrectionsRHSJacobTurbSubcellBlending");

//////////////////////////////////////////////////////////////////////////////
  
NSJacobBndGradientComputerSubcellBlending::NSJacobBndGradientComputerSubcellBlending(const std::string& name) :
  ConvBndCorrectionsRHSJacobFluxReconstructionSubcellBlending(name)
{
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD
