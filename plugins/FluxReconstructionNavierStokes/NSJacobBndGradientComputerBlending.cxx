// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/MethodCommandProvider.hh"

#include "FluxReconstructionNavierStokes/NSJacobBndGradientComputerBlending.hh"
#include "FluxReconstructionNavierStokes/FluxReconstructionNavierStokes.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< NSJacobBndGradientComputerBlending,
		       FluxReconstructionSolverData,
		       FluxReconstructionNavierStokesModule >
NSJacobBndGradientComputerBlendingProvider("ConvBndCorrectionsRHSJacobNSBlending");

MethodCommandProvider< NSJacobBndGradientComputerBlending,
		       FluxReconstructionSolverData,
		       FluxReconstructionNavierStokesModule >
NSJacobBndGradientComputerBlendingTurbProvider("ConvBndCorrectionsRHSJacobTurbBlending");
  
//////////////////////////////////////////////////////////////////////////////
  
NSJacobBndGradientComputerBlending::NSJacobBndGradientComputerBlending(const std::string& name) :
  ConvBndCorrectionsRHSJacobFluxReconstructionBlending(name),
  m_tempGradTerm(),
  m_tempGradTermGhost(),
  m_tempStates(),
  m_tempStatesGhost()
{
}

//////////////////////////////////////////////////////////////////////////////

void NSJacobBndGradientComputerBlending::setup()
{
  CFAUTOTRACE;
  
  // setup parent class
  ConvBndCorrectionsRHSJacobFluxReconstructionBlending::setup();
  
  m_tempGradTerm.resize(m_nbrEqs,m_nbrFaceFlxPnts);
  m_tempGradTermGhost.resize(m_nbrEqs,m_nbrFaceFlxPnts);
  
  m_tempStates.resize(m_nbrFaceFlxPnts);
  m_tempStatesGhost.resize(m_nbrFaceFlxPnts);
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD
