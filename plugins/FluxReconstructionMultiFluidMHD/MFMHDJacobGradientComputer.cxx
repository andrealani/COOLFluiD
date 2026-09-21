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

#include "FluxReconstructionMultiFluidMHD/MFMHDJacobGradientComputer.hh"
#include "FluxReconstructionMultiFluidMHD/FluxReconstructionMultiFluidMHD.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"
#include "MultiFluidMHD/MultiFluidMHDVarSet.hh"
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

MethodCommandProvider< MFMHDJacobGradientComputer,
		       FluxReconstructionSolverData,
		       FluxReconstructionMultiFluidMHDModule >
MFMHDJacobGradientComputerProvider("ConvRHSJacobMFMHD");
  
//////////////////////////////////////////////////////////////////////////////
  
MFMHDJacobGradientComputer::MFMHDJacobGradientComputer(const std::string& name) :
  ConvRHSJacobFluxReconstruction(name)
{
  addConfigOptionsTo(this);
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

