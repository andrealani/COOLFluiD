// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/MethodCommandProvider.hh"

#include "FluxReconstructionNavierStokes/NSGradientComputerSubcellBlending.hh"
#include "FluxReconstructionNavierStokes/FluxReconstructionNavierStokes.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< NSGradientComputerSubcellBlending,
		       FluxReconstructionSolverData,
		       FluxReconstructionNavierStokesModule >
nsGradientComputerSubcellBlendingProvider("ConvRHSNSSubcellBlending");
  
//////////////////////////////////////////////////////////////////////////////
  
NSGradientComputerSubcellBlending::NSGradientComputerSubcellBlending(const std::string& name) :
  ConvRHSFluxReconstructionSubcellBlending(name)
{
  // no addConfigOptionsTo(this) here: this class adds no options of its own, and
  // calling it would re-register the base class options (FaceFluxBlending) and
  // throw DuplicateNameException at configure time.
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD
