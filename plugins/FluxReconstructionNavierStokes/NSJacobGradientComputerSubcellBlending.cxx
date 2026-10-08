// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/MethodCommandProvider.hh"

#include "FluxReconstructionNavierStokes/NSJacobGradientComputerSubcellBlending.hh"
#include "FluxReconstructionNavierStokes/FluxReconstructionNavierStokes.hh"
#include "Framework/PhysicalModel.hh"
#include "NavierStokes/EulerTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< NSJacobGradientComputerSubcellBlending,
		       FluxReconstructionSolverData,
		       FluxReconstructionNavierStokesModule >
NSJacobGradientComputerSubcellBlendingProvider("ConvRHSJacobNSSubcellBlending");
  
//////////////////////////////////////////////////////////////////////////////
  
NSJacobGradientComputerSubcellBlending::NSJacobGradientComputerSubcellBlending(const std::string& name) :
  ConvRHSJacobFluxReconstructionSubcellBlending(name),
  m_machPData()
{
  // no addConfigOptionsTo(this) here: this class adds no options of its own, and
  // calling it would re-register the base class options (FaceFluxBlending) and
  // throw DuplicateNameException at configure time.
}

//////////////////////////////////////////////////////////////////////////////

CFreal NSJacobGradientComputerSubcellBlending::computeCellMaxMach(const std::vector< State* >& states)
{
  if (m_machPData.size() == 0)
  {
    PhysicalModelStack::getActive()->getImplementor()->getConvectiveTerm()->resizePhysicalData(m_machPData);
  }
  CFreal machMax = 0.;
  for (CFuint iSol = 0; iSol < states.size(); ++iSol)
  {
    m_updateVarSet->computePhysicalData(*states[iSol], m_machPData);
    machMax = std::max(machMax, m_machPData[EulerTerm::V]/m_machPData[EulerTerm::A]);
  }
  return machMax;
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD
