// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/MethodStrategyProvider.hh"

#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/NullInterfaceFlux.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider<NullInterfaceFlux,
				  FluxReconstructionSolverData,
				  BaseInterfaceFlux,
				  FluxReconstructionModule >
nullInterfaceFluxStrategyProvider("Null");

//////////////////////////////////////////////////////////////////////////////

NullInterfaceFlux::NullInterfaceFlux(const std::string& name) :
  BaseInterfaceFlux(name)
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

NullInterfaceFlux::~NullInterfaceFlux()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void NullInterfaceFlux::compute()
{
  CFAUTOTRACE;

  CFLog(VERBOSE, "NullInterfaceFlux::compute()\n");
}

//////////////////////////////////////////////////////////////////////////////

void NullInterfaceFlux::setup()
{
  CFAUTOTRACE;
  
  BaseInterfaceFlux::setup();

  CFLog(VERBOSE, "NullInterfaceFlux::setup()\n");
}

//////////////////////////////////////////////////////////////////////////////

void NullInterfaceFlux::unsetup()
{
  CFAUTOTRACE;
  
  BaseInterfaceFlux::unsetup();
  
  CFLog(VERBOSE, "NullInterfaceFlux::unsetup()\n");
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

