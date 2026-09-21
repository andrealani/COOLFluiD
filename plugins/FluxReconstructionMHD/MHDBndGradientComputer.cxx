// Copyright (C) 2022 KU Leuven CmPA, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/MethodCommandProvider.hh"

#include "Framework/CFSide.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/BaseTerm.hh"

#include "MathTools/MathFunctions.hh"

#include "FluxReconstructionMHD/MHDBndGradientComputer.hh"
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

MethodCommandProvider< MHDBndGradientComputer,
		       FluxReconstructionSolverData,
		       FluxReconstructionMHDModule >
MhdBndGradientComputerProvider("ConvBndCorrectionsRHSMHD");
  
//////////////////////////////////////////////////////////////////////////////
  
MHDBndGradientComputer::MHDBndGradientComputer(const std::string& name) :
  ConvBndCorrectionsRHSFluxReconstruction(name)
{
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

