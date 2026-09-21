// Copyright (C) 2019 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/MethodCommandProvider.hh"

#include "FluxReconstructionTurb/DiffBndCorrectionsRHSJacobFluxReconstructionGammaAlpha.hh"
#include "FluxReconstructionNavierStokes/FluxReconstructionNavierStokes.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< DiffBndCorrectionsRHSJacobFluxReconstructionGammaAlpha,
		       FluxReconstructionSolverData,
		       FluxReconstructionNavierStokesModule >
diffBndCorrectionsRHSJacobGammaAlphaFluxReconstructionProvider("DiffBndCorrectionsRHSJacobGammaAlpha");

//////////////////////////////////////////////////////////////////////////////

DiffBndCorrectionsRHSJacobFluxReconstructionGammaAlpha::DiffBndCorrectionsRHSJacobFluxReconstructionGammaAlpha(const std::string& name) :
  DiffBndCorrectionsRHSJacobFluxReconstructionTurb(name)
{
}

//////////////////////////////////////////////////////////////////////////////

DiffBndCorrectionsRHSJacobFluxReconstructionGammaAlpha::~DiffBndCorrectionsRHSJacobFluxReconstructionGammaAlpha()
{
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluxReconstructionMethod

} // namespace COOLFluiD
