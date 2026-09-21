// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/MethodCommandProvider.hh"

#include "FluxReconstructionNavierStokes/NSJacobBndGradientComputer.hh"
#include "FluxReconstructionNavierStokes/FluxReconstructionNavierStokes.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< NSJacobBndGradientComputer,
		       FluxReconstructionSolverData,
		       FluxReconstructionNavierStokesModule >
NSJacobBndGradientComputerProvider("ConvBndCorrectionsRHSJacobNS");

MethodCommandProvider< NSJacobBndGradientComputer,
		       FluxReconstructionSolverData,
		       FluxReconstructionNavierStokesModule >
NSJacobBndGradientComputerTurbProvider("ConvBndCorrectionsRHSJacobTurb");
  
//////////////////////////////////////////////////////////////////////////////
  
NSJacobBndGradientComputer::NSJacobBndGradientComputer(const std::string& name) :
  ConvBndCorrectionsRHSJacobFluxReconstruction(name),
  m_tempGradTerm(),
  m_tempGradTermGhost(),
  m_tempStates(),
  m_tempStatesGhost()
{
}

//////////////////////////////////////////////////////////////////////////////

void NSJacobBndGradientComputer::setup()
{
  CFAUTOTRACE;
  
  // setup parent class
  ConvBndCorrectionsRHSJacobFluxReconstruction::setup();
  
  m_tempGradTerm.resize(m_nbrEqs,m_nbrFaceFlxPnts);
  m_tempGradTermGhost.resize(m_nbrEqs,m_nbrFaceFlxPnts);
  
  m_tempStates.resize(m_nbrFaceFlxPnts);
  m_tempStatesGhost.resize(m_nbrFaceFlxPnts);
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD
