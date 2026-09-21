// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/MethodCommandProvider.hh"

#include "FluxReconstructionNavierStokes/NSJacobGradientComputerBlending.hh"
#include "FluxReconstructionNavierStokes/FluxReconstructionNavierStokes.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< NSJacobGradientComputerBlending,
		       FluxReconstructionSolverData,
		       FluxReconstructionNavierStokesModule >
NSJacobGradientComputerBlendingProvider("ConvRHSJacobNSBlending");

MethodCommandProvider< NSJacobGradientComputerBlending,
		       FluxReconstructionSolverData,
		       FluxReconstructionNavierStokesModule >
NSJacobGradientComputerBlendingTurbProvider("ConvRHSJacobTurbBlending");
  
//////////////////////////////////////////////////////////////////////////////
  
NSJacobGradientComputerBlending::NSJacobGradientComputerBlending(const std::string& name) :
  ConvRHSJacobFluxReconstructionBlending(name)
{
  addConfigOptionsTo(this);
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD
