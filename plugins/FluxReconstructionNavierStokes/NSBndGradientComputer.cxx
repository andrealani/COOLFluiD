// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/MethodCommandProvider.hh"

#include "FluxReconstructionNavierStokes/NSBndGradientComputer.hh"
#include "FluxReconstructionNavierStokes/FluxReconstructionNavierStokes.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< NSBndGradientComputer,
		       FluxReconstructionSolverData,
		       FluxReconstructionNavierStokesModule >
NsBndGradientComputerProvider("ConvBndCorrectionsRHSNS");

MethodCommandProvider< NSBndGradientComputer,
		       FluxReconstructionSolverData,
		       FluxReconstructionNavierStokesModule >
NSBndGradientComputerTurbProvider("ConvBndCorrectionsRHSTurb");
  
//////////////////////////////////////////////////////////////////////////////
  
NSBndGradientComputer::NSBndGradientComputer(const std::string& name) :
  ConvBndCorrectionsRHSFluxReconstruction(name),
  m_tempGradTerm(),
  m_tempGradTermGhost(),
  m_tempStates(),
  m_tempStatesGhost()
{
}

//////////////////////////////////////////////////////////////////////////////

void NSBndGradientComputer::setup()
{
  CFAUTOTRACE;
  
  // setup parent class
  ConvBndCorrectionsRHSFluxReconstruction::setup();
  
  m_tempGradTerm.resize(m_nbrEqs,m_nbrFaceFlxPnts);
  m_tempGradTermGhost.resize(m_nbrEqs,m_nbrFaceFlxPnts);
  
  m_tempStates.resize(m_nbrFaceFlxPnts);
  m_tempStatesGhost.resize(m_nbrFaceFlxPnts);
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD
