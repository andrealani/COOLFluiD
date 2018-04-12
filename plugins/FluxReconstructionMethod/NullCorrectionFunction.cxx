// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/MethodStrategyProvider.hh"

#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/NullCorrectionFunction.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider<NullCorrectionFunction,
				  FluxReconstructionSolverData,
				  BaseCorrectionFunction,
				  FluxReconstructionModule >
nullCorrectionFunctionStrategyProvider("Null");

//////////////////////////////////////////////////////////////////////////////

NullCorrectionFunction::NullCorrectionFunction(const std::string& name) :
  BaseCorrectionFunction(name)
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

NullCorrectionFunction::~NullCorrectionFunction()
{
  CFAUTOTRACE;
}


//////////////////////////////////////////////////////////////////////////////

void NullCorrectionFunction::computeCorrectionFunction(Common::SafePtr< FluxReconstructionElementData > frElemData, std::vector< std::vector< RealVector > >& corrcts)
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void NullCorrectionFunction::computeDivCorrectionFunction(Common::SafePtr< FluxReconstructionElementData > frElemData, std::vector< std::vector< CFreal > >& corrcts)
{
  CFAUTOTRACE;
}
      
//////////////////////////////////////////////////////////////////////////////
      
void NullCorrectionFunction::setup()
{
  CFAUTOTRACE;
  
  // setup parent class
  BaseCorrectionFunction::setup();
}

//////////////////////////////////////////////////////////////////////////////

void NullCorrectionFunction::unsetup()
{
  CFAUTOTRACE;
  
  // setup parent class
  BaseCorrectionFunction::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

